% Plot neural firing rate vs estimated time remaining for each cluster
% For experiment group 'training_running', calculate time remaining estimates
% based on cumulative progress and current speed, then bin by time remaining
% and compute average firing rates for each bin.
%
%   Specify options:
%
%       experiment_groups:      Will generate plots for all trials
%                               and all probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21',
%                                   'training_running'
%                               Should be a cell array of strings with each
%                               entry an experiment group
%
%       save_figs:              true or false, whether to save the figures to pdf
%
%       overwrite:              true or false. If figure pdf's already exist,
%                               whether to overwrite 
%       
%       figure_dir:             cell array of strings specifying which
%                               directory to save pdf's. The directory will
%                               be relative to the directory specified by
%                               path_config.figure_dir (in
%                               `path_config.m`), so that {'one', 'two',
%                               'three'} will save .pdfs to:
%                               <path_config.figure_dir>\one\two\three\      
%
%       time_remaining_bin_size: Size of time remaining bins in seconds (default: 2)
%       time_remaining_cap:      Maximum time remaining to include in analysis (default: 20)
%
% If `save_figs` is true, one pdf will be created for each probe recording,
% and contain:
% 1. Trial-level plots showing speed and time remaining for each trial
% 2. Cluster-level plots showing firing rate vs time remaining for each cluster

% Configuration
experiment_groups       = {'training_running'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'neural_firing_rate_vs_time_remaining', 'training_running'};
time_remaining_bin_size = .1; % seconds
time_remaining_cap      = 10; % seconds
gaussian_kernel_std     = 0.5; % seconds - standard deviation for Gaussian smoothing

% Initialize controller and get probe IDs
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));

for pid = 1:length(probe_ids)
    fprintf('\nProcessing probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    
    data = ctl.load_formatted_data(probe_ids{pid});
    clusters = data.selected_clusters();
    sessions = data.motion_sessions();
    
    fprintf('  Found %d clusters and %d sessions\n', length(clusters), length(sessions));
    
    % Collect all trials with their session index
    all_trials_info = struct('trial', {}, 'session_idx', {}, 'trial_idx', {});
    for s = 1:length(sessions)
        session = sessions{s};
        trials_in_session = session.trials;
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            if any(motion_mask) % Only include trials with motion periods
                all_trials_info(end+1) = struct('trial', trial, 'session_idx', s, 'trial_idx', t); %#ok<AGROW>
            end
        end
    end
    
    fprintf('  Total trials with motion: %d\n', numel(all_trials_info));
    
    % Group trials by individual sessions
    session_trials = {};
    for s = 1:length(sessions)
        session_trials{s} = all_trials_info([all_trials_info.session_idx] == s);
    end
    
    % Define colors for different sessions
    % 120cm sessions: dark blue and cyan
    % 60cm sessions: dark red and orange
    session_colors = {[0 0 0.8], [0 0.8 0.8], [0.8 0 0], [1 0.5 0]}; % dark blue, cyan, dark red, orange
    session_labels = {};
    for s = 1:length(sessions)
        if s <= 2  % First two sessions are 120cm
            session_labels{s} = sprintf('120cm Session %d', s);
        else  % Remaining sessions are 60cm
            session_labels{s} = sprintf('60cm Session %d', s);
        end
    end
    
    % --- 1. Create trial-level plots showing speed and time remaining ---
    fprintf('  Creating trial-level plots...\n');
    for t = 1:length(all_trials_info)
        trial_info = all_trials_info(t);
        trial = trial_info.trial;
        
        % Get motion periods only
        motion_mask = trial.motion_mask();
        if ~any(motion_mask)
            continue;
        end
        
        % Extract motion data
        tvec = trial.probe_t(motion_mask);
        vel = trial.velocity(motion_mask);
        
        % Calculate total expected distance
        dt_vec = [diff(tvec); 0];
        total_distance = sum(vel .* dt_vec);
        
        % Calculate cumulative distance and time remaining at each time point
        cumulative_distance = cumsum(vel .* dt_vec);
        remaining_distance = total_distance - cumulative_distance;
        
        % Calculate time remaining (remaining_distance / current_speed)
        time_remaining = remaining_distance ./ vel;
        time_remaining(vel < 0.1) = NaN; % Avoid division by very small speeds
        
        % Remove unrealistic time remaining outliers
        % Cap time remaining and remove values above the cap
        time_remaining(time_remaining > time_remaining_cap) = NaN;
        
        % Create figure for this trial
        fig_trial = ctl.figs.a4figure('landscape');
        
        % Get the time range for consistent x-axis limits
        time_range = [min(tvec), max(tvec)];
        
        % Top panel: Running speed over time
        subplot(2, 1, 1);
        plot(tvec, vel, 'b-', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('Speed (cm/s)');
        title(sprintf('Trial %d - Speed Profile', trial.trial_id));
        xlim(time_range);
        grid on;
        
        % Bottom panel: Time remaining over time
        subplot(2, 1, 2);
        plot(tvec, time_remaining, 'r-', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('Time Remaining (s)');
        title(sprintf('Trial %d - Estimated Time to Completion', trial.trial_id));
        xlim(time_range);
        grid on;
        
        % Add figure title
        FigureTitle(fig_trial, sprintf('Trial %d - %s', trial.trial_id, probe_ids{pid}));
        ctl.figs.save_fig_to_join();
    end
    
    % --- 2. Create cluster-level plots showing firing rate vs time remaining ---
    fprintf('  Creating cluster-level plots...\n');
    
    % For each cluster, analyze firing rate vs time remaining
    for c = 1:length(clusters)
        cluster = clusters(c);
        
        % Create figure for this cluster with two subplots
        fig_cluster = ctl.figs.a4figure('landscape');
        
        % Store data for both plots
        all_session_data = struct('bin_centers', {}, 'median_firing_rates', {}, ...
                                 'p25_firing_rates', {}, 'p75_firing_rates', {}, ...
                                 'median_firing_rates_norm', {}, 'smoothed_firing_rates', {});
        
        % Process each session separately
        for s = 1:length(sessions)
            session_trials_this_session = session_trials{s};
            
            if isempty(session_trials_this_session)
                continue;
            end
            
            % Collect all time remaining values and corresponding firing rates
            all_time_remaining = [];
            all_firing_rates = [];
            
            for t = 1:length(session_trials_this_session)
                trial_info = session_trials_this_session(t);
                trial = trial_info.trial;
                
                % Get motion periods only
                motion_mask = trial.motion_mask();
                if ~any(motion_mask)
                    continue;
                end
                
                % Extract motion data
                tvec = trial.probe_t(motion_mask);
                vel = trial.velocity(motion_mask);
                
                % Calculate total expected distance
                dt_vec = [diff(tvec); 0];
                total_distance = sum(vel .* dt_vec);
                
                % Calculate cumulative distance and time remaining
                cumulative_distance = cumsum(vel .* dt_vec);
                remaining_distance = total_distance - cumulative_distance;
                time_remaining = remaining_distance ./ vel;
                time_remaining(vel < 0.1) = NaN;
                
                % Remove unrealistic time remaining outliers
                % Cap time remaining and remove values above the cap
                time_remaining(time_remaining > time_remaining_cap) = NaN;
                
                % Get firing rate for this cluster during motion periods
                fr = cluster.fr.get_convolution(tvec);
                
                % Remove NaN values
                valid_mask = ~isnan(time_remaining) & ~isnan(fr);
                if any(valid_mask)
                    all_time_remaining = [all_time_remaining; time_remaining(valid_mask)];
                    all_firing_rates = [all_firing_rates; fr(valid_mask)];
                end
            end
            
            if isempty(all_time_remaining)
                continue;
            end
            
            % Bin time remaining values
            min_time_remaining = min(all_time_remaining);
            max_time_remaining = max(all_time_remaining);
            bin_edges = min_time_remaining:time_remaining_bin_size:max_time_remaining;
            if isempty(bin_edges) || length(bin_edges) < 2
                bin_edges = [min_time_remaining, max_time_remaining];
            end
            
            % Calculate median firing rate and percentiles for each bin
            [~, ~, bin_indices] = histcounts(all_time_remaining, bin_edges);
            bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
            
            median_firing_rates = zeros(1, length(bin_centers));
            p25_firing_rates = zeros(1, length(bin_centers));
            p75_firing_rates = zeros(1, length(bin_centers));
            n_samples_per_bin = zeros(1, length(bin_centers));
            
            for b = 1:length(bin_centers)
                bin_mask = bin_indices == b;
                if any(bin_mask)
                    bin_firing_rates = all_firing_rates(bin_mask);
                    median_firing_rates(b) = median(bin_firing_rates);
                    p25_firing_rates(b) = prctile(bin_firing_rates, 25);
                    p75_firing_rates(b) = prctile(bin_firing_rates, 75);
                    n_samples_per_bin(b) = sum(bin_mask);
                end
            end
            
            % Normalize median firing rates to 0-1 range (independently for each session)
            if max(median_firing_rates) > min(median_firing_rates)
                median_firing_rates_norm = (median_firing_rates - min(median_firing_rates)) / (max(median_firing_rates) - min(median_firing_rates));
            else
                median_firing_rates_norm = zeros(size(median_firing_rates));
            end
            
            % Create Gaussian kernel for smoothing
            kernel_std = gaussian_kernel_std / time_remaining_bin_size; % Convert to bin units
            kernel_size = round(6 * kernel_std); % 6 std for kernel width
            if mod(kernel_size, 2) == 0
                kernel_size = kernel_size + 1; % Ensure odd size
            end
            kernel_center = (kernel_size + 1) / 2;
            kernel = exp(-0.5 * ((1:kernel_size) - kernel_center).^2 / kernel_std^2);
            kernel = kernel / sum(kernel); % Normalize
            
            % Pad the rate data to avoid edge effects using mirror padding
            pad_size = floor(length(kernel) / 2);
            rate_padded = [fliplr(median_firing_rates_norm(1:pad_size)), median_firing_rates_norm, fliplr(median_firing_rates_norm(end-pad_size+1:end))];
            rate_smooth_padded = conv(rate_padded, kernel, 'same');
            smoothed_firing_rates = rate_smooth_padded(pad_size + 1 : end - pad_size);
            
            % Store data for this session
            all_session_data(s).bin_centers = bin_centers;
            all_session_data(s).median_firing_rates = median_firing_rates;
            all_session_data(s).p25_firing_rates = p25_firing_rates;
            all_session_data(s).p75_firing_rates = p75_firing_rates;
            all_session_data(s).median_firing_rates_norm = median_firing_rates_norm;
            all_session_data(s).smoothed_firing_rates = smoothed_firing_rates;
        end
        
        % Plot 1: Original firing rates with percentiles
        subplot(2, 1, 1);
        for s = 1:length(sessions)
            if ~isempty(all_session_data) && s <= length(all_session_data) && ~isempty(all_session_data(s).bin_centers)
                bin_centers = all_session_data(s).bin_centers;
                median_firing_rates = all_session_data(s).median_firing_rates;
                p25_firing_rates = all_session_data(s).p25_firing_rates;
                p75_firing_rates = all_session_data(s).p75_firing_rates;
                
                % Plot this session
                if s == 1
                    plot(bin_centers, median_firing_rates, 'Color', session_colors{s}, 'LineWidth', 2, 'DisplayName', session_labels{s});
                    hold on;
                else
                    plot(bin_centers, median_firing_rates, 'Color', session_colors{s}, 'LineWidth', 2, 'DisplayName', session_labels{s});
                end
                
                % Add 25th/75th percentile shading
                fill([bin_centers, fliplr(bin_centers)], ...
                     [p75_firing_rates, fliplr(p25_firing_rates)], ...
                     session_colors{s}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', sprintf('%s (25th-75th percentile)', session_labels{s}));
            end
        end
        
        hold off;
        xlabel('Time Remaining (s)');
        ylabel('Median Firing Rate (Hz)');
        title(sprintf('Cluster %d - %s', cluster.id, probe_ids{pid}));
        legend('show');
        grid on;
        
        % Plot 2: Normalized firing rates with smoothing
        subplot(2, 1, 2);
        for s = 1:length(sessions)
            if ~isempty(all_session_data) && s <= length(all_session_data) && ~isempty(all_session_data(s).bin_centers)
                bin_centers = all_session_data(s).bin_centers;
                median_firing_rates_norm = all_session_data(s).median_firing_rates_norm;
                smoothed_firing_rates = all_session_data(s).smoothed_firing_rates;
                
                % Plot this session
                if s == 1
                    h1 = plot(bin_centers, median_firing_rates_norm, 'Color', session_colors{s}, 'LineWidth', 1, 'DisplayName', session_labels{s});
                    set(h1, 'Color', [session_colors{s}, 0.5]); % Set transparency
                    hold on;
                    plot(bin_centers, smoothed_firing_rates, 'Color', session_colors{s}, 'LineWidth', 2, 'DisplayName', sprintf('%s (smoothed)', session_labels{s}));
                else
                    h1 = plot(bin_centers, median_firing_rates_norm, 'Color', session_colors{s}, 'LineWidth', 1, 'DisplayName', session_labels{s});
                    set(h1, 'Color', [session_colors{s}, 0.5]); % Set transparency
                    plot(bin_centers, smoothed_firing_rates, 'Color', session_colors{s}, 'LineWidth', 2, 'DisplayName', sprintf('%s (smoothed)', session_labels{s}));
                end
            end
        end
        
        hold off;
        xlabel('Time Remaining (s)');
        ylabel('Normalized Firing Rate (0-1)');
        title(sprintf('Cluster %d - %s (Normalized)', cluster.id, probe_ids{pid}));
        legend('show');
        grid on;
        
        % Add text with sample counts
        sample_text = '';
        for s = 1:length(sessions)
            if ~isempty(session_trials{s})
                sample_text = [sample_text, sprintf('%s: %d, ', session_labels{s}, length(session_trials{s}))];
            end
        end
        if ~isempty(sample_text)
            sample_text = sample_text(1:end-2); % Remove last comma and space
        end
        text(0.02, 0.98, sprintf('Samples: %s', sample_text), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
        
        % Add figure title
        FigureTitle(fig_cluster, sprintf('Cluster %d - %s', cluster.id, probe_ids{pid}));
        ctl.figs.save_fig_to_join();
    end
    
    % Join all plots for this probe into one PDF
    fprintf('  Saving PDF for probe %s...\n', probe_ids{pid});
    fname = sprintf('%s.pdf', probe_ids{pid});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
    
    fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
end

fprintf('\nNeural firing rate vs time remaining analysis completed for %d probe(s)\n', length(probe_ids)); 