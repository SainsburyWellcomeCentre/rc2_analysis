% Plot neural firing rate vs planned trajectory progress for each cluster
% For experiment group 'training_running', fit smooth mathematical models to speed data
% to estimate planned trajectories, then calculate progress based on planned vs actual
% movement and compute average firing rates for each progress bin.
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
%       progress_bin_size:       Size of progress bins in % (default: 5)
%       progress_range:          Range of progress to include in analysis (default: [0, 100])
%       fit_method:              Method for fitting planned trajectory:
%                                   'parabolic' - fit quadratic function
%                                   'gaussian' - fit Gaussian function
%                                   'spline' - fit cubic spline
%
% If `save_figs` is true, one pdf will be created for each probe recording,
% and contain:
% 1. Trial-level plots showing actual vs planned speed and progress for each trial
% 2. Cluster-level plots showing firing rate vs planned progress for each cluster

% Configuration
experiment_groups       = {'training_running'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'neural_firing_rate_vs_planned_trajectory', 'training_running'};
progress_bin_size       = 5; % percentage
progress_range          = [0, 100]; % percentage range
fit_method              = 'parabolic'; % 'parabolic', 'gaussian', or 'spline'

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
    
    % Determine session types based on maximum position in each session
    session_types = cell(length(sessions), 1);
    session_max_positions = zeros(length(sessions), 1);
    
    for s = 1:length(sessions)
        session = sessions{s};
        trials_in_session = session.trials;
        session_max_pos = 0;
        
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            if ~isempty(pos)
                max_pos = max(pos);
                session_max_pos = max(session_max_pos, max_pos);
            end
        end
        
        session_max_positions(s) = session_max_pos;
        if session_max_pos > 90
            session_types{s} = 'long';
        else
            session_types{s} = 'short';
        end
        
        fprintf('  Session %d: max position = %.2f cm, type = %s\n', s, session_max_pos, session_types{s});
    end
    
    % Define colors for different sessions
    % Long sessions (120cm): cold colors (blue, cyan)
    % Short sessions (60cm): warm colors (red, orange)
    cold_colors = {[0 0 0.8], [0 0.8 0.8]}; % dark blue, cyan
    warm_colors = {[0.8 0 0], [1 0.5 0]}; % dark red, orange
    
    session_colors = cell(length(sessions), 1);
    session_labels = {};
    long_count = 0;
    short_count = 0;
    
    for s = 1:length(sessions)
        if strcmp(session_types{s}, 'long')
            long_count = long_count + 1;
            session_colors{s} = cold_colors{mod(long_count-1, length(cold_colors)) + 1};
            session_labels{s} = sprintf('120cm Session %d', s);
        else
            short_count = short_count + 1;
            session_colors{s} = warm_colors{mod(short_count-1, length(warm_colors)) + 1};
            session_labels{s} = sprintf('60cm Session %d', s);
        end
    end
    
    % --- 1. Create trial-level plots showing actual vs planned speed and progress ---
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
        
        % Fit planned trajectory
        [planned_vel, planned_pos, total_planned_distance] = fit_planned_trajectory(tvec, vel, fit_method);
        
        % Calculate progress based on planned trajectory
        progress = calculate_planned_progress(tvec, vel, planned_vel, planned_pos, total_planned_distance);
        
        % Create figure for this trial
        fig_trial = ctl.figs.a4figure('landscape');
        
        % Get the time range for consistent x-axis limits
        time_range = [min(tvec), max(tvec)];
        
        % Top panel: Actual vs planned speed over time
        subplot(3, 1, 1);
        plot(tvec, vel, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Actual Speed');
        hold on;
        plot(tvec, planned_vel, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Planned Speed');
        xlabel('Time (s)');
        ylabel('Speed (cm/s)');
        title(sprintf('Trial %d - Speed Profile (%s fit)', trial.trial_id, fit_method));
        xlim(time_range);
        legend('show');
        grid on;
        hold off;
        
        % Middle panel: Actual vs planned position over time
        subplot(3, 1, 2);
        actual_pos = cumsum(vel .* [diff(tvec); 0]);
        plot(tvec, actual_pos, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Actual Position');
        hold on;
        plot(tvec, planned_pos, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Planned Position');
        xlabel('Time (s)');
        ylabel('Position (cm)');
        title(sprintf('Trial %d - Position Profile', trial.trial_id));
        xlim(time_range);
        legend('show');
        grid on;
        hold off;
        
        % Bottom panel: Progress based on planned trajectory
        subplot(3, 1, 3);
        plot(tvec, progress, 'g-', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('Progress (%)');
        title(sprintf('Trial %d - Planned Trajectory Progress', trial.trial_id));
        xlim(time_range);
        ylim([0, 100]);
        grid on;
        
        % Add figure title
        FigureTitle(fig_trial, sprintf('Trial %d - %s (%s)', trial.trial_id, probe_ids{pid}, fit_method));
        ctl.figs.save_fig_to_join();
    end
    
    % --- 2. Create cluster-level plots showing firing rate vs planned progress ---
    fprintf('  Creating cluster-level plots...\n');
    
    % For each cluster, analyze firing rate vs planned progress
    for c = 1:length(clusters)
        cluster = clusters(c);
        
        % Create figure for this cluster
        fig_cluster = ctl.figs.a4figure('portrait');
        
        % Store data for plotting
        all_session_data = struct('bin_centers', {}, 'median_firing_rates', {}, ...
                                 'p25_firing_rates', {}, 'p75_firing_rates', {});
        
        % Process each session separately
        for s = 1:length(sessions)
            session_trials_this_session = session_trials{s};
            
            if isempty(session_trials_this_session)
                continue;
            end
            
            % Collect all progress values and corresponding firing rates
            all_progress = [];
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
                
                % Fit planned trajectory
                [planned_vel, planned_pos, total_planned_distance] = fit_planned_trajectory(tvec, vel, fit_method);
                
                % Calculate progress based on planned trajectory
                progress = calculate_planned_progress(tvec, vel, planned_vel, planned_pos, total_planned_distance);
                
                % Get firing rate for this cluster during motion periods
                fr = cluster.fr.get_convolution(tvec);
                
                % Remove NaN values
                valid_mask = ~isnan(progress) & ~isnan(fr);
                if any(valid_mask)
                    all_progress = [all_progress; progress(valid_mask)];
                    all_firing_rates = [all_firing_rates; fr(valid_mask)];
                end
            end
            
            if isempty(all_progress)
                continue;
            end
            
            % Bin progress values
            min_progress = max(min(all_progress), progress_range(1));
            max_progress = min(max(all_progress), progress_range(2));
            bin_edges = min_progress:progress_bin_size:max_progress;
            if isempty(bin_edges) || length(bin_edges) < 2
                bin_edges = [min_progress, max_progress];
            end
            
            % Calculate median firing rate and percentiles for each bin
            [~, ~, bin_indices] = histcounts(all_progress, bin_edges);
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
            
            % Store data for this session
            all_session_data(s).bin_centers = bin_centers;
            all_session_data(s).median_firing_rates = median_firing_rates;
            all_session_data(s).p25_firing_rates = p25_firing_rates;
            all_session_data(s).p75_firing_rates = p75_firing_rates;
        end
        
        % Plot firing rates with percentiles
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
        xlabel('Planned Trajectory Progress (%)');
        ylabel('Median Firing Rate (Hz)');
        title(sprintf('Cluster %d - %s (%s fit)', cluster.id, probe_ids{pid}, fit_method));
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
        FigureTitle(fig_cluster, sprintf('Cluster %d - %s (%s)', cluster.id, probe_ids{pid}, fit_method));
        ctl.figs.save_fig_to_join();
    end
    
    % Join all plots for this probe into one PDF
    fprintf('  Saving PDF for probe %s...\n', probe_ids{pid});
    fname = sprintf('%s_%s.pdf', probe_ids{pid}, fit_method);
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
    
    fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
end

fprintf('\nNeural firing rate vs planned trajectory analysis completed for %d probe(s)\n', length(probe_ids));

% Helper function to fit planned trajectory
function [planned_vel, planned_pos, total_planned_distance] = fit_planned_trajectory(tvec, vel, method)
    % Normalize time to [0, 1] for better fitting
    t_norm = (tvec - tvec(1)) / (tvec(end) - tvec(1));
    
    % Constrain fit to periods where running speed is > 2 cm/s
    running_mask = vel > 2.0;
    
    % If we don't have enough running periods, fall back to original method
    if sum(running_mask) < 3
        % Remove outliers for fitting (speeds that are too high or too low)
        valid_mask = vel > 0.1 & vel < prctile(vel, 95);
        t_fit = t_norm(valid_mask);
        vel_fit = vel(valid_mask);
    else
        % Use only running periods for fitting
        t_fit = t_norm(running_mask);
        vel_fit = vel(running_mask);
    end
    
    if length(t_fit) < 3
        % Not enough data for fitting, use simple average
        planned_vel = repmat(mean(vel), size(tvec));
        planned_pos = cumsum(planned_vel) * mean(diff(tvec));
        total_planned_distance = planned_pos(end);
        return;
    end
    
    switch method
        case 'parabolic'
            % Fit quadratic function: v(t) = a*t^2 + b*t + c
            p = polyfit(t_fit, vel_fit, 2);
            planned_vel = polyval(p, t_norm);
            
        case 'gaussian'
            % Fit Gaussian function: v(t) = a*exp(-((t-b)/c)^2)
            % Use initial guess based on data
            [~, max_idx] = max(vel_fit);
            a_guess = max(vel_fit);
            b_guess = t_fit(max_idx);
            c_guess = 0.3;
            
            % Define Gaussian function
            gaussian_fun = @(params, t) params(1) * exp(-((t - params(2)) / params(3)).^2);
            
            % Fit using lsqcurvefit
            try
                params0 = [a_guess, b_guess, c_guess];
                params = lsqcurvefit(gaussian_fun, params0, t_fit, vel_fit);
                planned_vel = gaussian_fun(params, t_norm);
            catch
                % Fallback to parabolic if Gaussian fitting fails
                p = polyfit(t_fit, vel_fit, 2);
                planned_vel = polyval(p, t_norm);
            end
            
        case 'spline'
            % Fit cubic spline
            try
                pp = spline(t_fit, vel_fit);
                planned_vel = ppval(pp, t_norm);
            catch
                % Fallback to parabolic if spline fitting fails
                p = polyfit(t_fit, vel_fit, 2);
                planned_vel = polyval(p, t_norm);
            end
            
        otherwise
            error('Unknown fit method: %s', method);
    end
    
    % Ensure planned velocity is non-negative
    planned_vel = max(planned_vel, 0.1);
    
    % Calculate planned position and total distance
    dt_vec = [diff(tvec); 0];
    planned_pos = cumsum(planned_vel .* dt_vec);
    total_planned_distance = planned_pos(end);
end

% Helper function to calculate progress based on planned trajectory
function progress = calculate_planned_progress(tvec, vel, planned_vel, planned_pos, total_planned_distance)
    % Calculate actual cumulative distance
    dt_vec = [diff(tvec); 0];
    actual_pos = cumsum(vel .* dt_vec);
    
    % Calculate progress as ratio of actual to planned distance
    progress = (actual_pos / total_planned_distance) * 100;
    
    % Ensure progress stays within reasonable bounds
    progress = max(0, min(100, progress));
end 