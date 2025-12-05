 % Plot spatial firing rate profile for each cluster using position as x-axis
% For experiment group 'training_running', bin position into 2-cm bins,
% normalize spike counts by occupancy, and smooth with a Gaussian kernel (8-cm s.d.)
%
% REFACTORED VERSION using SpatialTuningAnalyzer class
%
% SMOOTHING APPROACH:
%   Spike counts and occupancy are smoothed SEPARATELY before computing rates.
%   This is mathematically cleaner than smoothing rates directly:
%       rate_smooth = smooth(counts) / smooth(occupancy)
%   rather than:
%       rate_smooth = smooth(counts / occupancy)
%
% STATISTICAL TESTING:
%   Spatial tuning significance is tested using a bin-shuffle approach analogous
%   to velocity tuning in ShuffleTuning.m. The rate matrix (n_trials × n_bins) is
%   completely shuffled to break trial-bin associations, and Skaggs information
%   is computed on the shuffled data to create a null distribution.
%
% Requires:
%   - lib/classes/analysis/SpatialTuningAnalyzer.m
%   - lib/spatial_analysis/* helper functions

%%
% Prevent figures from popping up on screen
set(groot, 'DefaultFigureVisible', 'off');

% Configuration
experiment_groups        = {'ambient_light'};
save_figs                = true;
overwrite                = true;
figure_dir               = {'spatial_firing_rate', 'ambient_light'};
plot_single_cluster_fig  = true;
plot_heatmap_cluster_fig = true;
re_run_analysis          = false;  % Set to true to recompute all metrics, false to load cached data

% Parallel processing configuration
use_parallel            = true;   % Set to false to disable parallel processing entirely
max_workers             = 4;      % Maximum number of parallel workers

% Analysis parameters
bin_size_cm = 2;
gauss_sigma_cm = 8; % standard deviation for smoothing (cm)
position_threshold_cm = 90; % threshold for long/short trial classification

% Speed tuning parameters
trial_group_label_for_tuning = 'RT';  % Must match the label used in create_tables.m

% Acceleration tuning parameters
trial_group_label_for_accel_tuning = 'RT';  % Must match the label used in create_tables_acceleration.m

% Statistical analysis parameters
run_statistics = true;      % Set to false to skip statistical analysis
minSpikes = 25;             % minimum spikes per cluster per group
minPeakRate = 1.0;          % Hz minimum peak for field consideration
fieldFrac = 0.7;            % field threshold fraction of rate range
minFieldBins = 5;           % contiguous bins for a field
maxNumFields = 1;           % maximum number of fields
nShuf = 100;                % number of bin-shuffles for significance testing
pThresh = 0.05;             % significance threshold

% Initialize controller and get probe IDs
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));

for pid = 1:length(probe_ids)
    fprintf('\nProcessing probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    
    % Define cache file path
    cache_filename = sprintf('%s_spatial_analysis_cache.mat', probe_ids{pid});
    cache_filepath = fullfile(ctl.figs.curr_dir, cache_filename);
    
    % Load data
    data = ctl.load_formatted_data(probe_ids{pid});
    clusters = data.selected_clusters();
    sessions = data.motion_sessions();
    
    fprintf('  Found %d clusters and %d sessions\n', length(clusters), length(sessions));
    
    % Create analyzer instance
    analyzer = SpatialTuningAnalyzer(clusters, sessions);
    
    % Configure parameters
    analyzer.bin_size_cm = bin_size_cm;
    analyzer.gauss_sigma_cm = gauss_sigma_cm;
    analyzer.position_threshold = position_threshold_cm;
    analyzer.use_parallel = use_parallel;
    analyzer.max_workers = max_workers;
    
    % Configure statistical parameters
    analyzer.stats_params.minSpikes = minSpikes;
    analyzer.stats_params.minPeakRate = minPeakRate;
    analyzer.stats_params.fieldFrac = fieldFrac;
    analyzer.stats_params.minFieldBins = minFieldBins;
    analyzer.stats_params.maxNumFields = maxNumFields;
    analyzer.stats_params.nShuf = nShuf;
    analyzer.stats_params.pThresh = pThresh;
    analyzer.stats_params.run_statistics = run_statistics;
    
    % Try to load cached data if not re-running analysis
    if ~re_run_analysis && exist(cache_filepath, 'file')
        fprintf('  Loading cached analysis data from: %s\n', cache_filename);
        analyzer.load_cache(cache_filepath);
        fprintf('  Loaded cached data for %d clusters\n', length(analyzer.cluster_ids));
        
        % Still need trial_groups for plotting
        analyzer.collect_trials();
        
        % Bin config is reconstructed in load_cache, but verify it exists
        if isempty(analyzer.bin_config)
            analyzer.setup_bins();
        end
    else
        if ~re_run_analysis
            fprintf('  No cache found, computing analysis...\n');
        else
            fprintf('  Re-running analysis (re_run_analysis=true)...\n');
        end
        
        % Run full analysis
        analyzer.analyze_all_clusters();
        
        % Print summary
        analyzer.print_summary();
        
        % Save cache
        fprintf('  Saving analysis data to cache: %s\n', cache_filename);
        analyzer.save_cache(cache_filepath);
        fprintf('  Cache saved successfully\n');
    end
    
    % Extract data for plotting (convert to old format for compatibility)
    cluster_ids = analyzer.cluster_ids;
    n_clusters = length(cluster_ids);
    group_names = analyzer.group_names;
    group_labels = analyzer.group_labels;
    trial_groups = analyzer.trial_groups;
    
    all_bin_centers_by_group = struct();
    all_rate_smooth_by_group = struct();
    all_Q1_rate_smooth_by_group = struct();
    all_Q2_rate_smooth_by_group = struct();
    all_Q3_rate_smooth_by_group = struct();
    all_avg_velocity_by_group = struct();
    all_occ_by_group = struct();
    all_rate_per_trial_by_group = struct();
    all_rate_per_trial_smooth_by_group = struct();
    all_spike_positions_by_group = struct();
    all_global_trial_indices_by_group = struct();
    all_rate_pooled_by_group = struct();
    spatial_tuning_stats = analyzer.tuning_stats;
    
    for g = 1:2
        group = group_names{g};
        all_bin_centers_by_group.(group) = analyzer.bin_config.(group).plot_centers;
        all_rate_smooth_by_group.(group) = analyzer.firing_rates.(group).rate_smooth;
        all_Q1_rate_smooth_by_group.(group) = analyzer.firing_rates.(group).Q1_smooth;
        all_Q2_rate_smooth_by_group.(group) = analyzer.firing_rates.(group).Q2_smooth;
        all_Q3_rate_smooth_by_group.(group) = analyzer.firing_rates.(group).Q3_smooth;
        all_avg_velocity_by_group.(group) = analyzer.firing_rates.(group).avg_velocity;
        all_occ_by_group.(group) = analyzer.firing_rates.(group).occupancy;
        all_rate_per_trial_by_group.(group) = analyzer.firing_rates.(group).rate_per_trial;
        all_rate_per_trial_smooth_by_group.(group) = analyzer.firing_rates.(group).rate_per_trial_smooth;
        all_spike_positions_by_group.(group) = analyzer.firing_rates.(group).spike_positions;
        all_global_trial_indices_by_group.(group) = analyzer.firing_rates.(group).global_trial_indices;
        all_rate_pooled_by_group.(group) = analyzer.firing_rates.(group).rate_pooled;
    end
    
    % --- Plot average running and occupancy profiles for each group ---
    if save_figs
        fprintf('  Creating average running and occupancy profiles...\n');
        fig_profiles = ctl.figs.a4figure('landscape');
        
        % Reconstruct per-trial velocity and occupancy data for plotting
        trial_vel_by_group = struct();
        trial_occ_by_group = struct();
        
        for g = 1:2
            group = group_names{g};
            trials_struct = trial_groups.(group);
            edges = analyzer.bin_config.(group).edges;
            bin_centers = analyzer.bin_config.(group).centers;
            n_bins = analyzer.bin_config.(group).n_bins;
            
            n_trials = length(trials_struct);
            vel_per_trial = nan(n_trials, n_bins);
            occ_per_trial = nan(n_trials, n_bins);
            
            for k = 1:n_trials
                trial = trials_struct(k).trial;
                [~, occ, vel, ~] = compute_trial_firing_rate(trial, clusters(1), edges);
                vel_per_trial(k, :) = vel;
                occ_per_trial(k, :) = occ;
            end
            
            trial_vel_by_group.(group) = vel_per_trial;
            trial_occ_by_group.(group) = occ_per_trial;
        end
        
        % Average velocity plot with individual trials
        subplot(2, 1, 1, 'Parent', fig_profiles);
        hold on;
        for g = 1:2
            group = group_names{g};
            bin_centers = all_bin_centers_by_group.(group);
            
            % Plot individual trials with thin, semi-transparent lines
            vel_trials = trial_vel_by_group.(group);
            for k = 1:size(vel_trials, 1)
                if g == 1  % long trials - blue
                    plot(bin_centers, vel_trials(k, :), 'Color', [0.5 0.5 1], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                else  % short trials - red
                    plot(bin_centers, vel_trials(k, :), 'Color', [1 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                end
            end
            
            % Plot average with thick line
            avg_vel = nanmean(all_avg_velocity_by_group.(group), 1);
            if g == 1  % long trials - blue
                plot(bin_centers, avg_vel, 'Color', [0 0 0.8], 'LineWidth', 2, 'DisplayName', group_labels{g});
            else  % short trials - red
                plot(bin_centers, avg_vel, 'Color', [0.8 0 0], 'LineWidth', 2, 'DisplayName', group_labels{g});
            end
        end
        hold off;
        xlabel('Position (cm)');
        ylabel('Avg Velocity (cm/s)');
        legend('show', 'Location', 'northeastoutside');
        grid on;
        
        % Average occupancy plot with individual trials
        subplot(2, 1, 2, 'Parent', fig_profiles);
        hold on;
        for g = 1:2
            group = group_names{g};
            bin_centers = all_bin_centers_by_group.(group);
            
            % Plot individual trials with thin, semi-transparent lines
            occ_trials = trial_occ_by_group.(group);
            for k = 1:size(occ_trials, 1)
                if g == 1  % long trials - blue
                    plot(bin_centers, occ_trials(k, :), 'Color', [0.5 0.5 1], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                else  % short trials - red
                    plot(bin_centers, occ_trials(k, :), 'Color', [1 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                end
            end
            
            % Plot average with thick line
            avg_occ = nanmean(all_occ_by_group.(group), 1);
            if g == 1  % long trials - blue
                plot(bin_centers, avg_occ, 'Color', [0 0 0.8], 'LineWidth', 2, 'DisplayName', group_labels{g});
            else  % short trials - red
                plot(bin_centers, avg_occ, 'Color', [0.8 0 0], 'LineWidth', 2, 'DisplayName', group_labels{g});
            end
        end
        hold off;
        xlabel('Position (cm)');
        ylabel('Occupancy (s)');
        legend('show', 'Location', 'northeastoutside');
        grid on;
        
        % Add figure title
        FigureTitle(fig_profiles, sprintf('Average Profiles - %s', probe_ids{pid}));
        ctl.figs.save_fig_to_join();
    end

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
    end
    
    % Define colors for different sessions
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

    % For each cluster, create a combined figure with spatial profile and x-normalized comparison
    if plot_single_cluster_fig
        fprintf('  Creating individual cluster plots...\n');
        
        % Define colors for the two groups
        group_colors = struct('long', [0 0 0.8], 'short', [0.8 0 0]); % Blue for long, red for short
        
        % Load speed tuning curves for all clusters (combined long+short trials)
        tuning_curves = cell(n_clusters, 1);
        for c = 1:n_clusters
            try
                tuning_curves{c} = data.load_tuning_curves(cluster_ids(c), trial_group_label_for_tuning);
            catch ME
                fprintf('    Warning: Could not load speed tuning curve for cluster %d\n', cluster_ids(c));
                fprintf('    Error: %s\n', ME.message);
                tuning_curves{c} = [];
            end
        end
        
        % Load acceleration tuning curves for all clusters (combined long+short trials)
        accel_tuning_curves = cell(n_clusters, 1);
        for c = 1:n_clusters
            try
                % load_tuning_curves_acceleration returns a cell array {all, acc, dec}
                % We want the 'all' table (index 1) for combined acceleration tuning
                accel_data = data.load_tuning_curves_acceleration(cluster_ids(c), trial_group_label_for_accel_tuning);
                if ~isempty(accel_data) && length(accel_data) >= 1
                    accel_tuning_curves{c} = accel_data{1};  % Extract 'all' table
                else
                    accel_tuning_curves{c} = [];
                end
            catch ME
                fprintf('    Warning: Could not load acceleration tuning curve for cluster %d\n', cluster_ids(c));
                fprintf('    Error: %s\n', ME.message);
                accel_tuning_curves{c} = [];
            end
        end
        
        for c = 1:n_clusters
            cluster = clusters(c);
            
            % Prepare rate data structure for this cluster
            rate_data = struct();
            for g = 1:2
                group = group_names{g};
                rate_data.(group).Q1_smooth = all_Q1_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).Q2_smooth = all_Q2_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).Q3_smooth = all_Q3_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).rate_per_trial = all_rate_per_trial_by_group.(group){c};
                rate_data.(group).rate_per_trial_smooth = all_rate_per_trial_smooth_by_group.(group){c};
                rate_data.(group).spike_positions = all_spike_positions_by_group.(group){c};
                rate_data.(group).global_trial_indices = all_global_trial_indices_by_group.(group){c};
                rate_data.(group).rate_pooled = all_rate_pooled_by_group.(group)(c, :);
            end
            
            % Get statistics for this cluster (if available)
            cluster_stats = [];
            if ~isempty(spatial_tuning_stats)
                cluster_field = sprintf('cluster_%d', cluster.id);
                if isfield(spatial_tuning_stats, cluster_field)
                    cluster_stats = spatial_tuning_stats.(cluster_field);
                end
            end
            
            % Compute 2D rate maps (position × velocity/acceleration) with equal-occupancy binning
            rate_maps_2d = struct();
            
            % Velocity: use 20 equal-occupancy bins (5% per bin, like tuning curves)
            n_vel_bins = 20;
            sigma_vel = 5;  % Smoothing sigma for velocity (cm/s)
            min_occ = 0.1;  % Minimum occupancy threshold per trial (seconds)
            
            % Compute for long trials
            [rate_maps_2d.vel_long, vel_bin_edges_long, vel_trial_counts_long] = ...
                analyzer.compute_2d_position_velocity_tuning(cluster, 'long', n_vel_bins, sigma_vel, min_occ);
            % Compute for short trials  
            [rate_maps_2d.vel_short, vel_bin_edges_short, vel_trial_counts_short] = ...
                analyzer.compute_2d_position_velocity_tuning(cluster, 'short', n_vel_bins, sigma_vel, min_occ);
            
            % Store bin edges, centers, and trial counts for plotting
            rate_maps_2d.vel_bin_edges_long = vel_bin_edges_long;
            rate_maps_2d.vel_bin_edges_short = vel_bin_edges_short;
            rate_maps_2d.vel_bin_centers_long = (vel_bin_edges_long(1:end-1) + vel_bin_edges_long(2:end)) / 2;
            rate_maps_2d.vel_bin_centers_short = (vel_bin_edges_short(1:end-1) + vel_bin_edges_short(2:end)) / 2;
            rate_maps_2d.vel_trial_counts_long = vel_trial_counts_long;
            rate_maps_2d.vel_trial_counts_short = vel_trial_counts_short;
            rate_maps_2d.vel_n_trials_long = length(trial_groups.long);
            rate_maps_2d.vel_n_trials_short = length(trial_groups.short);
            
            % Acceleration: use 20 equal-occupancy bins (5% per bin, like tuning curves)
            n_accel_bins = 20;
            sigma_accel = 10;  % Smoothing sigma for acceleration (cm/s²) - reduced to preserve structure
            
            % Compute for long trials
            [rate_maps_2d.accel_long, accel_bin_edges_long, accel_trial_counts_long] = ...
                analyzer.compute_2d_position_acceleration_tuning(cluster, 'long', n_accel_bins, sigma_accel, min_occ);
            % Compute for short trials
            [rate_maps_2d.accel_short, accel_bin_edges_short, accel_trial_counts_short] = ...
                analyzer.compute_2d_position_acceleration_tuning(cluster, 'short', n_accel_bins, sigma_accel, min_occ);
            
            % Store bin edges, centers, and trial counts for plotting
            rate_maps_2d.accel_bin_edges_long = accel_bin_edges_long;
            rate_maps_2d.accel_bin_edges_short = accel_bin_edges_short;
            rate_maps_2d.accel_bin_centers_long = (accel_bin_edges_long(1:end-1) + accel_bin_edges_long(2:end)) / 2;
            rate_maps_2d.accel_bin_centers_short = (accel_bin_edges_short(1:end-1) + accel_bin_edges_short(2:end)) / 2;
            rate_maps_2d.accel_trial_counts_long = accel_trial_counts_long;
            rate_maps_2d.accel_trial_counts_short = accel_trial_counts_short;
            rate_maps_2d.accel_n_trials_long = length(trial_groups.long);
            rate_maps_2d.accel_n_trials_short = length(trial_groups.short);
            
            % Store position bins (both centers and edges) for axis
            rate_maps_2d.pos_bins_long = analyzer.bin_config.long.centers;
            rate_maps_2d.pos_bins_short = analyzer.bin_config.short.centers;
            rate_maps_2d.pos_bin_edges_long = analyzer.bin_config.long.edges;
            rate_maps_2d.pos_bin_edges_short = analyzer.bin_config.short.edges;
            
            % Create combined figure using helper function
            fig_cluster = plot_cluster_spatial_profile(cluster.id, all_bin_centers_by_group, ...
                                                       rate_data, group_names, group_labels, ...
                                                       group_colors, probe_ids{pid}, cluster_stats, tuning_curves{c}, accel_tuning_curves{c}, rate_maps_2d);
            
            % Save using the figure management system
            if save_figs
                % Convert to a4figure format for PDF joining
                set(fig_cluster, 'PaperOrientation', 'landscape');
                set(fig_cluster, 'PaperUnits', 'normalized');
                set(fig_cluster, 'PaperPosition', [0 0 1 1]);
                ctl.figs.save_fig_to_join();
            else
                close(fig_cluster);
            end
        end
    end
    
    % Plot heatmaps for each group in a single figure with four subplots
    if plot_heatmap_cluster_fig
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
            rate_mat = all_Q2_rate_smooth_by_group.(group);  % Use median (Q2) instead of mean
            bin_centers = all_bin_centers_by_group.(group);
            sorted_rate_mat = rate_mat(sort_idx, :);
            
            % Raw heatmap
            subplot(2, 2, g, 'Parent', fig_heatmaps);
            imagesc(bin_centers, 1:n_clusters, sorted_rate_mat);
            colormap('jet');
            colorbar;
            xlabel('Position (cm)');
            ylabel('Cluster (sorted by maximum peak position)');
            set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sorted_cluster_ids);
            
            % Set x-axis limits to make short trials visually half width
            if strcmp(group, 'short')
                xlim([60, 120]);
            else
                xlim([0, 120]);
            end
            
            % Normalized heatmap (min-max normalization)
            min_vals = min(sorted_rate_mat, [], 2);
            max_vals = max(sorted_rate_mat, [], 2);
            norm_rate_mat = (sorted_rate_mat - min_vals) ./ (max_vals - min_vals);
            norm_rate_mat(isnan(norm_rate_mat) | isinf(norm_rate_mat)) = 0;
            subplot(2, 2, g+2, 'Parent', fig_heatmaps);
            imagesc(bin_centers, 1:n_clusters, norm_rate_mat);
            colormap('jet');
            colorbar;
            xlabel('Position (cm)');
            ylabel('Cluster (sorted by maximum peak position)');
            set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sorted_cluster_ids);
            
            % Set x-axis limits
            if strcmp(group, 'short')
                xlim([60, 120]);
            else
                xlim([0, 120]);
            end
        end
        
        % Adjust subplot positions BEFORE adding title
        all_axes = findall(fig_heatmaps, 'Type', 'axes');
        for ax_idx = 1:length(all_axes)
            ax = all_axes(ax_idx);
            pos = get(ax, 'Position');
            pos(4) = pos(4) * 0.82;  % Reduce height
            pos(2) = max(0.08, pos(2) - 0.03);  % Lower position
            set(ax, 'Position', pos);
        end
        
        % Add figure title
        FigureTitle(fig_heatmaps, sprintf('Heatmaps - %s', probe_ids{pid}));
        ctl.figs.save_fig_to_join();
        
        % Plot spatially tuned clusters
        if ~isempty(spatial_tuning_stats)
            fprintf('  Creating spatially tuned clusters plot...\n');
            
            % Identify clusters tuned in at least one condition
            spatially_tuned_any = [];
            
            for c = 1:n_clusters
                cluster_field = sprintf('cluster_%d', cluster_ids(c));
                if isfield(spatial_tuning_stats, cluster_field)
                    stats = spatial_tuning_stats.(cluster_field);
                    if (isfield(stats, 'long') && stats.long.isSpatiallyTuned) || ...
                       (isfield(stats, 'short') && stats.short.isSpatiallyTuned)
                        spatially_tuned_any = [spatially_tuned_any; c]; %#ok<AGROW>
                    end
                end
            end
            
            n_tuned = length(spatially_tuned_any);
            fprintf('    Found %d clusters spatially tuned in at least one condition\n', n_tuned);
            
            % Create 2x3 grid plot using helper function
            % Create 2x3 grid plot
            fig_grid = plot_spatial_tuning_grid(spatial_tuning_stats, ...
                all_bin_centers_by_group, all_Q2_rate_smooth_by_group, ...
                cluster_ids, gauss_sigma_cm, probe_ids{pid});
            
            if ~isempty(fig_grid)
                ctl.figs.save_fig_to_join();
            end
            
            % Create 3D visualization
            fig_3d = plot_spatial_tuning_3d(spatial_tuning_stats, ...
                all_bin_centers_by_group, all_Q2_rate_smooth_by_group, ...
                cluster_ids, probe_ids{pid}, ctl);
            
            if ~isempty(fig_3d)
                ctl.figs.save_fig_to_join();
            end
        end
    end
    
    % Join figures if requested
    if save_figs
        fprintf('  Joining figures into single PDF...\n');
        fname = sprintf('%s.pdf', probe_ids{pid});
        ctl.figs.join_figs(fname, overwrite);
        ctl.figs.clear_figs();
    end
    
    fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
    
    % Clean up memory
    clear analyzer data clusters sessions trial_groups
    clear all_rate_smooth_by_group all_bin_centers_by_group all_avg_velocity_by_group
    clear all_Q1_rate_smooth_by_group all_Q2_rate_smooth_by_group all_Q3_rate_smooth_by_group
    clear all_occ_by_group all_rate_per_trial_by_group all_rate_per_trial_smooth_by_group
    clear all_spike_positions_by_group all_global_trial_indices_by_group all_rate_pooled_by_group
    clear spatial_tuning_stats
    
    fprintf('  Memory cleaned up after probe %d/%d\n', pid, length(probe_ids));
end

fprintf('\n=== All probes processed ===\n');

% Restore default figure visibility
set(groot, 'DefaultFigureVisible', 'on');
