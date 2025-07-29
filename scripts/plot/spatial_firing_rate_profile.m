% Plot spatial firing rate profile for each cluster using position as x-axis
% For experiment group 'training_running', bin position into 2-cm bins,
% normalize spike counts by occupancy, and smooth with a Gaussian kernel (8-cm s.d.)
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
% If `save_figs` is true, one pdf will be created for each probe recording,
% and contain spatial firing rate profiles for all clusters in that probe.

% Configuration
experiment_groups       = {'training_running'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'spatial_firing_rate', 'training_running'};

% Analysis parameters
bin_size_cm = 2;
gauss_sigma_cm = 8; % standard deviation for smoothing (cm)

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
    
    % Session analysis and trial counting
    for s = 1:length(sessions)
        session = sessions{s};
        fprintf('  Session %d/%d: session_id = %s\n', s, length(sessions), session.session_id);
        trials_in_session = session.trials;
        n_in_motion_trials = 0;
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            if ~isempty(pos)
                n_in_motion_trials = n_in_motion_trials + 1;
                max_pos = max(pos);
                fprintf('    Trial %d: max(position in-motion) = %.2f cm\n', t, max_pos);
            else
                fprintf('    Trial %d: [no in-motion samples]\n', t);
            end
        end
        fprintf('    Total in-motion trials in session: %d\n', n_in_motion_trials);
    end

    % --- Collect all trials with their session index and max in-motion position ---
    all_trials_info = struct('trial', {}, 'session_idx', {}, 'trial_idx', {}, 'max_pos', {});
    for s = 1:length(sessions)
        session = sessions{s};
        trials_in_session = session.trials;
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            if ~isempty(pos)
                max_pos = max(pos);
                all_trials_info(end+1) = struct('trial', trial, 'session_idx', s, 'trial_idx', t, 'max_pos', max_pos); %#ok<AGROW>
            end
        end
    end
    
    % --- Group trials by max position (long/short) ---
    long_trials = all_trials_info([all_trials_info.max_pos] > 90);
    short_trials = all_trials_info([all_trials_info.max_pos] <= 90);
    fprintf('  Total long trials: %d\n', numel(long_trials));
    fprintf('  Total short trials: %d\n', numel(short_trials));

    % --- Prepare for cluster analysis as before, but using long_trials and short_trials ---
    trial_groups = struct('long', long_trials, 'short', short_trials);
    group_names = {'long', 'short'};
    group_labels = {'Long (0-120 cm)', 'Short (60-120 cm)'};
    all_rate_smooth_by_group = struct();
    all_rate_unsmooth_by_group = struct(); % NEW: Store unsmoothed rates
    all_bin_centers_by_group = struct();
    all_avg_velocity_by_group = struct();
    all_occ_by_group = struct();
    cluster_ids = [clusters.id];
    n_clusters = length(clusters);

    % Compute and fill group summary structures for all clusters (needed for later plots)
    fprintf('  Computing spatial firing rate profiles for %d clusters...\n', n_clusters);
    for c = 1:n_clusters
        cluster = clusters(c);
        for g = 1:2
            group = group_names{g};
            trials_struct = trial_groups.(group);
            all_spike_positions = [];
            all_positions = [];
            all_times = [];
            all_velocities = [];
            for k = 1:length(trials_struct)
                trial = trials_struct(k).trial;
                motion_mask = trial.motion_mask();
                pos = trial.position(motion_mask);
                tvec = trial.probe_t(motion_mask);
                vel = trial.velocity(motion_mask);
                all_positions = [all_positions; pos(:)];
                all_times = [all_times; tvec(:)];
                all_velocities = [all_velocities; vel(:)];
                st = cluster.spike_times;
                mask = st >= tvec(1) & st <= tvec(end);
                st = st(mask);
                if ~isempty(st)
                    spike_pos = interp1(tvec, pos, st, 'linear', 'extrap');
                    all_spike_positions = [all_spike_positions; spike_pos];
                end
            end
            if isempty(all_positions)
                continue;
            end
            % Set bin edges explicitly for each group
            if strcmp(group, 'long')
                edges = 0:bin_size_cm:120;
            else
                edges = 0:bin_size_cm:60; % Short trials have positions 0-60 cm
            end
            bin_centers = edges(1:end-1) + bin_size_cm/2;
            
            % For short trials, adjust bin centers for plotting to align with 60-120 cm
            if strcmp(group, 'short')
                plot_bin_centers = bin_centers + 60; % Shift to 60-120 cm for plotting
            else
                plot_bin_centers = bin_centers;
            end
            
            spike_count = histcounts(all_spike_positions, edges);
            dt_vec = [diff(all_times); 0];
            occupancy = zeros(1, length(edges)-1);
            velocity_sum = zeros(1, length(edges)-1);
            for i = 1:length(occupancy)
                in_bin = all_positions >= edges(i) & all_positions < edges(i+1);
                occupancy(i) = sum(dt_vec(in_bin));
                velocity_sum(i) = sum(all_velocities(in_bin) .* dt_vec(in_bin));
            end
            avg_velocity = velocity_sum ./ occupancy;
            avg_velocity(isnan(avg_velocity) | isinf(avg_velocity)) = 0;
            rate = spike_count ./ occupancy;
            rate(isnan(rate) | isinf(rate)) = 0;
            
            % Create proper spatial Gaussian kernel
            kernel_size_cm = 8 * gauss_sigma_cm; % 64 cm total kernel size
            kernel_samples = round(kernel_size_cm / bin_size_cm); 
            x = linspace(-kernel_samples/2, kernel_samples/2, kernel_samples);
            kernel = exp(-x.^2 / (2 * (gauss_sigma_cm/bin_size_cm)^2));
            kernel = kernel / sum(kernel); % normalize
            
            % Pad the rate data to avoid edge effects using mirror padding
            pad_size = floor(length(kernel) / 2);
            rate_padded = [fliplr(rate(1:pad_size)), rate, fliplr(rate(end-pad_size+1:end))];
            rate_smooth_padded = conv(rate_padded, kernel, 'same');
            rate_smooth = rate_smooth_padded(pad_size + 1 : end - pad_size);
            
            % Store for heatmap and group plots
            if ~isfield(all_rate_smooth_by_group, group)
                all_rate_smooth_by_group.(group) = nan(n_clusters, length(rate_smooth));
                all_rate_unsmooth_by_group.(group) = nan(n_clusters, length(rate)); % NEW: Store unsmoothed rates
                all_bin_centers_by_group.(group) = plot_bin_centers; % Use plot_bin_centers for display
                all_avg_velocity_by_group.(group) = nan(n_clusters, length(avg_velocity));
                all_occ_by_group.(group) = nan(n_clusters, length(occupancy));
            end
            all_rate_smooth_by_group.(group)(c, :) = rate_smooth;
            all_rate_unsmooth_by_group.(group)(c, :) = rate; % NEW: Store unsmoothed rates
            all_bin_centers_by_group.(group) = plot_bin_centers; % Use plot_bin_centers for display
            all_avg_velocity_by_group.(group)(c, :) = avg_velocity;
            all_occ_by_group.(group)(c, :) = occupancy;
        end
    end

    % --- Plot average running and occupancy profiles for each group ---
    fprintf('  Creating average running and occupancy profiles...\n');
    fig_profiles = ctl.figs.a4figure('landscape');
    % Average velocity
    subplot(2, 1, 1, 'Parent', fig_profiles);
    hold on;
    for g = 1:2
        group = group_names{g};
        bin_centers = all_bin_centers_by_group.(group);
        avg_vel = nanmean(all_avg_velocity_by_group.(group), 1);
        plot(bin_centers, avg_vel, 'LineWidth', 2, 'DisplayName', group_labels{g});
    end
    hold off;
    xlabel('Position (cm)');
    ylabel('Avg Velocity (cm/s)');
    legend('show');
    % Removed title to avoid overlap with FigureTitle
    grid on;
    % Average occupancy
    subplot(2, 1, 2, 'Parent', fig_profiles);
    hold on;
    for g = 1:2
        group = group_names{g};
        bin_centers = all_bin_centers_by_group.(group);
        avg_occ = nanmean(all_occ_by_group.(group), 1);
        plot(bin_centers, avg_occ, 'LineWidth', 2, 'DisplayName', group_labels{g});
    end
    hold off;
    xlabel('Position (cm)');
    ylabel('Occupancy (s)');
    legend('show');
    % Removed title to avoid overlap with FigureTitle
    grid on;
    
    % Add figure title
    FigureTitle(fig_profiles, sprintf('Average Profiles - %s', probe_ids{pid}));
    ctl.figs.save_fig_to_join();

    % Determine session types based on actual trial data
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

    % For each cluster, create a single figure with both group firing rates overlaid
    fprintf('  Creating individual cluster plots...\n');
    for c = 1:n_clusters
        cluster = clusters(c);
        fig_cluster = ctl.figs.a4figure('landscape');
        
        % Process each session separately
        for s = 1:length(sessions)
            session = sessions{s};
            trials_in_session = session.trials;
            
            % Collect data for this session
            all_spike_positions = [];
            all_positions = [];
            all_times = [];
            all_velocities = [];
            
            for t = 1:length(trials_in_session)
                trial = trials_in_session{t}.to_aligned;
                motion_mask = trial.motion_mask();
                pos = trial.position(motion_mask);
                tvec = trial.probe_t(motion_mask);
                vel = trial.velocity(motion_mask);
                all_positions = [all_positions; pos(:)];
                all_times = [all_times; tvec(:)];
                all_velocities = [all_velocities; vel(:)];
                st = cluster.spike_times;
                mask = st >= tvec(1) & st <= tvec(end);
                st = st(mask);
                if ~isempty(st)
                    spike_pos = interp1(tvec, pos, st, 'linear', 'extrap');
                    all_spike_positions = [all_spike_positions; spike_pos];
                end
            end
            
            if isempty(all_positions)
                continue;
            end
            
            % Set bin edges based on session type (long vs short)
            if strcmp(session_types{s}, 'long')
                edges = 0:bin_size_cm:120;
            else  % Short sessions
                edges = 0:bin_size_cm:60;
            end
            bin_centers = edges(1:end-1) + bin_size_cm/2;
            
            % For short trials, adjust bin centers for plotting to align with 60-120 cm
            if strcmp(session_types{s}, 'short')
                plot_bin_centers = bin_centers + 60; % Shift to 60-120 cm for plotting
            else
                plot_bin_centers = bin_centers;
            end
            
            spike_count = histcounts(all_spike_positions, edges);
            dt_vec = [diff(all_times); 0];
            occupancy = zeros(1, length(edges)-1);
            velocity_sum = zeros(1, length(edges)-1);
            for i = 1:length(occupancy)
                in_bin = all_positions >= edges(i) & all_positions < edges(i+1);
                occupancy(i) = sum(dt_vec(in_bin));
                velocity_sum(i) = sum(all_velocities(in_bin) .* dt_vec(in_bin));
            end
            avg_velocity = velocity_sum ./ occupancy;
            avg_velocity(isnan(avg_velocity) | isinf(avg_velocity)) = 0;
            rate = spike_count ./ occupancy;
            rate(isnan(rate) | isinf(rate)) = 0;
            
            % Create proper spatial Gaussian kernel
            kernel_size_cm = 8 * gauss_sigma_cm; % 64 cm total kernel size
            kernel_samples = round(kernel_size_cm / bin_size_cm); 
            x = linspace(-kernel_samples/2, kernel_samples/2, kernel_samples);
            kernel = exp(-x.^2 / (2 * (gauss_sigma_cm/bin_size_cm)^2));
            kernel = kernel / sum(kernel); % normalize
            
            % Pad the rate data to avoid edge effects using mirror padding
            pad_size = floor(length(kernel) / 2);
            rate_padded = [fliplr(rate(1:pad_size)), rate, fliplr(rate(end-pad_size+1:end))];
            rate_smooth_padded = conv(rate_padded, kernel, 'same');
            rate_smooth = rate_smooth_padded(pad_size + 1 : end - pad_size);
            
            % Plot this session
            hold on;
            % Plot unsmoothed rate with 0.5 alpha (semi-transparent)
            plot(plot_bin_centers, rate, '-', 'LineWidth', 1, 'Color', [session_colors{s}, 0.5], 'DisplayName', sprintf('%s (unsmoothed)', session_labels{s}));
            % Plot smoothed rate (solid line)
            plot(plot_bin_centers, rate_smooth, '-', 'LineWidth', 2, 'Color', session_colors{s}, 'DisplayName', session_labels{s});
        end
        
        hold off;
        xlabel('Position (cm)');
        ylabel('Firing rate (Hz)');
        legend('show');
        % Removed title to avoid overlap with FigureTitle
        grid on;
        
        % Add figure title
        FigureTitle(fig_cluster, sprintf('Cluster %d - %s', cluster.id, probe_ids{pid}));
        ctl.figs.save_fig_to_join();
    end

    % Plot heatmaps for each group in a single figure with four subplots
    fprintf('  Creating heatmap plots...\n');
    fig_heatmaps = ctl.figs.a4figure('landscape');

    % Sort clusters by maximum peak positional firing rate for long trials
    long_rate_mat = all_rate_smooth_by_group.('long');
    peak_positions = zeros(n_clusters, 1);
    for c = 1:n_clusters
        rate_profile = long_rate_mat(c, :);
        bin_centers = all_bin_centers_by_group.('long');
        
        % Find the position of the maximum peak
        [~, max_idx] = max(rate_profile);
        peak_positions(c) = bin_centers(max_idx);
    end
    [~, sort_idx] = sort(peak_positions);
    sorted_cluster_ids = cluster_ids(sort_idx);

    for g = 1:2
        group = group_names{g};
        rate_mat = all_rate_smooth_by_group.(group);
        bin_centers = all_bin_centers_by_group.(group);
        sorted_rate_mat = rate_mat(sort_idx, :);
        
        % Raw heatmap
        subplot(2, 2, g, 'Parent', fig_heatmaps);
        imagesc(bin_centers, 1:n_clusters, sorted_rate_mat);
        colormap('jet'); % Changed to blue-to-yellow colormap
        colorbar;
        xlabel('Position (cm)');
        ylabel('Cluster (sorted by maximum peak position)');
        % Removed title to avoid overlap with FigureTitle
        set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sorted_cluster_ids);
        
        % Set x-axis limits to make short trials visually half width
        if strcmp(group, 'short')
            xlim([60, 120]); % Short trials span 60-120 cm
        else
            xlim([0, 120]); % Long trials span 0-120 cm
        end
        
        % Normalized heatmap
        norm_rate_mat = sorted_rate_mat ./ max(sorted_rate_mat, [], 2);
        norm_rate_mat(isnan(norm_rate_mat)) = 0;
        subplot(2, 2, g+2, 'Parent', fig_heatmaps);
        imagesc(bin_centers, 1:n_clusters, norm_rate_mat);
        colormap('jet'); % Changed to blue-to-yellow colormap
        colorbar;
        xlabel('Position (cm)');
        ylabel('Cluster (sorted by maximum peak position)');
        % Removed title to avoid overlap with FigureTitle
        set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sorted_cluster_ids);
        
        % Set x-axis limits to make short trials visually half width
        if strcmp(group, 'short')
            xlim([60, 120]); % Short trials span 60-120 cm
        else
            xlim([0, 120]); % Long trials span 0-120 cm
        end
    end
    % Removed sgtitle to avoid overlap with FigureTitle
    
    % Add figure title with cleaner formatting
    FigureTitle(fig_heatmaps, sprintf('Heatmaps - %s', probe_ids{pid}));
    ctl.figs.save_fig_to_join();
    
    % Join all plots for this probe into one PDF
    fprintf('  Saving PDF for probe %s...\n', probe_ids{pid});
    fname = sprintf('%s.pdf', probe_ids{pid});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
    
    fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
end

fprintf('\nSpatial firing rate profile analysis completed for %d probe(s)\n', length(probe_ids)); 