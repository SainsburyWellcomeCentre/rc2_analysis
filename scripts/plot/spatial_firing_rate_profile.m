 % Plot spatial firing rate profile for each cluster using position as x-axis
% For experiment group 'training_running', bin position into 2-cm bins,
% normalize spike counts by occupancy, and smooth with a Gaussian kernel (8-cm s.d.)
%
% REFACTORED VERSION using SpatialTuningAnalyzer class and modular helper functions
%
% HELPER FUNCTIONS (in lib/spatial_analysis/):
%   - load_spatial_analysis_cache.m: Load analyzer from cache file
%   - extract_plotting_data_from_analyzer.m: Extract organized plotting data
%   - plot_average_velocity_occupancy_profiles.m: Create velocity/occupancy plots
%   - load_tuning_curves_for_clusters.m: Load velocity/acceleration tuning curves
%   - prepare_cluster_data_for_plotting.m: Prepare rate and 2D map data per cluster
%   - compute_ttg_kinematic_maps.m: Compute 2D TTG × velocity/acceleration maps
%   - create_single_cluster_figures.m: Create individual spatial profile figures for each cluster
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
% Prevent figures from popping up on screen and stealing focus
% Reset any previously set CreateFcn that might cause errors
set(groot, 'DefaultFigureCreateFcn', '');  % Clear any previous CreateFcn
set(groot, 'DefaultFigureVisible', 'off');

% Configuration
experiment_groups        = {'ambient_light'};
save_figs                = true;
overwrite                = true;
figure_dir               = {'spatial_firing_rate', 'ambient_light'};
plot_single_cluster_fig  = true;
plot_heatmap_cluster_fig = true;
re_run_analysis          = false;  % Set to true to recompute all metrics, false to load cached data
plot_velocity_tuning     = true;   % Set to false to skip velocity tuning plots
plot_acceleration_tuning = true;   % Set to false to skip acceleration tuning plots

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
fieldFrac = 0.6;            % field threshold fraction of rate range
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
    
    % Create parameters struct for analysis
    params = struct();
    params.bin_size_cm = bin_size_cm;
    params.gauss_sigma_cm = gauss_sigma_cm;
    params.position_threshold = position_threshold_cm;
    params.use_parallel = use_parallel;
    params.max_workers = max_workers;
    params.stats_params = struct();
    params.stats_params.minSpikes = minSpikes;
    params.stats_params.minPeakRate = minPeakRate;
    params.stats_params.fieldFrac = fieldFrac;
    params.stats_params.minFieldBins = minFieldBins;
    params.stats_params.maxNumFields = maxNumFields;
    params.stats_params.nShuf = nShuf;
    params.stats_params.pThresh = pThresh;
    params.stats_params.run_statistics = run_statistics;
    
    % Decide whether to load from cache or run analysis
    if re_run_analysis || ~exist(cache_filepath, 'file')
        % No cache or force recompute: run full analysis
        if ~exist(cache_filepath, 'file')
            fprintf('  No cache found\n');
        else
            fprintf('  Re-running analysis (re_run_analysis=true)\n');
        end
        analyzer = SpatialTuningAnalyzer(clusters, sessions, params);
        analyzer.run_full_analysis_and_save(cache_filepath);
    else
        % Load from existing cache
        analyzer = load_spatial_analysis_cache(clusters, sessions, cache_filepath);
    end
    
    % Initialize TTG rate storage for heatmap plotting
    all_ttg_rates = struct('long', [], 'short', []);
    ttg_norm_bin_centers = [];
    
    % --- Plot average running and occupancy profiles for each group ---
    if save_figs
        fprintf('  Creating average running and occupancy profiles...\n');
        fig_profiles = plot_average_velocity_occupancy_profiles(analyzer.trial_groups, analyzer.group_names, analyzer.group_labels, ...
            analyzer, clusters, probe_ids{pid}, ctl);
        ctl.figs.save_fig_to_join();
    end
    
    % --- Plot Time-to-Goal analysis figure ---
    if save_figs
        fprintf('  Creating time-to-goal analysis figure...\n');
        fig_ttg_analysis = plot_time_to_goal_analysis(analyzer.trial_groups, analyzer.group_names, analyzer.group_labels, ...
            analyzer.bin_config, probe_ids{pid}, ctl);
        ctl.figs.save_fig_to_join();
    end
    
    % Verify distribution comparison results
    if ~isempty(analyzer.distribution_comparisons)
        n_non_empty = sum(~cellfun(@isempty, analyzer.distribution_comparisons));
        fprintf('  Using distribution comparisons (%d/%d clusters have results)\n', ...
            n_non_empty, length(analyzer.distribution_comparisons));
    else
        fprintf('  WARNING: No distribution comparisons found!\n');
    end

    % Create individual cluster figures
    [all_ttg_rates, ttg_norm_bin_centers, all_ttg_fields] = create_single_cluster_figures(...
        plot_single_cluster_fig, analyzer, clusters, data, ...
        plot_velocity_tuning, plot_acceleration_tuning, ...
        trial_group_label_for_tuning, trial_group_label_for_accel_tuning, ...
        bin_size_cm, gauss_sigma_cm, probe_ids{pid}, save_figs, ctl);
    
    % Plot distribution comparison contingency table
    if save_figs
        fprintf('  Creating distribution comparison contingency table...\n');
        fig_contingency = plot_distribution_comparison_contingency(analyzer.distribution_comparisons, analyzer.cluster_ids, probe_ids{pid});
        
        % Convert to a4figure format for PDF joining
        set(fig_contingency, 'PaperOrientation', 'landscape');
        set(fig_contingency, 'PaperUnits', 'normalized');
        set(fig_contingency, 'PaperPosition', [0 0 1 1]);
        ctl.figs.save_fig_to_join();
    end
    
    % Plot heatmaps for each group in a single figure with four subplots
    if plot_heatmap_cluster_fig
        fprintf('  Creating heatmap plots...\n');
        fig_heatmaps = ctl.figs.a4figure('landscape');
    
        % Sort clusters by maximum peak positional firing rate for long trials
        long_rate_mat = analyzer.firing_rates.long.Q2_smooth;
        peak_positions = zeros(length(analyzer.cluster_ids), 1);
        for c = 1:length(analyzer.cluster_ids)
            rate_profile = long_rate_mat(c, :);
            bin_centers = analyzer.bin_config.long.plot_centers;
            
            % Find the position of the maximum peak
            [~, max_idx] = max(rate_profile);
            peak_positions(c) = bin_centers(max_idx);
        end
        [~, sort_idx] = sort(peak_positions);
        
        % Determine significance for each cluster in each condition
        is_significant_long = false(length(analyzer.cluster_ids), 1);
        is_significant_short = false(length(analyzer.cluster_ids), 1);
        
        if ~isempty(analyzer.tuning_stats)
            for c = 1:length(analyzer.cluster_ids)
                cluster_field = sprintf('cluster_%d', analyzer.cluster_ids(c));
                if isfield(analyzer.tuning_stats, cluster_field)
                    stats = analyzer.tuning_stats.(cluster_field);
                    if isfield(stats, 'long') && stats.long.isSpatiallyTuned
                        is_significant_long(c) = true;
                    end
                    if isfield(stats, 'short') && stats.short.isSpatiallyTuned
                        is_significant_short(c) = true;
                    end
                end
            end
        end
        
        % Sort significance vectors to match sorted clusters
        sorted_is_significant_long = is_significant_long(sort_idx);
        sorted_is_significant_short = is_significant_short(sort_idx);
    
        for g = 1:2
            group = analyzer.group_names{g};
            rate_mat = analyzer.firing_rates.(group).Q2_smooth;  % Use median (Q2) instead of mean
            bin_centers = analyzer.bin_config.(group).plot_centers;
            sorted_rate_mat = rate_mat(sort_idx, :);
            
            % Get significance for this group
            if strcmp(group, 'long')
                is_sig = sorted_is_significant_long;
            else
                is_sig = sorted_is_significant_short;
            end
            
            % Min-max normalization
            min_vals = min(sorted_rate_mat, [], 2);
            max_vals = max(sorted_rate_mat, [], 2);
            norm_rate_mat = (sorted_rate_mat - min_vals) ./ (max_vals - min_vals);
            norm_rate_mat(isnan(norm_rate_mat) | isinf(norm_rate_mat)) = 0;
            
            % Separate significant and non-significant clusters, PRESERVING the global sort order
            % (already sorted by long peak position)
            sig_rate_mat = norm_rate_mat(is_sig, :);
            nonsig_rate_mat = norm_rate_mat(~is_sig, :);
            n_sig = sum(is_sig);
            n_nonsig = sum(~is_sig);
            
            % DO NOT re-sort within significant or non-significant groups
            % The order from the global long-condition sorting is preserved
            
            % Significant clusters (colored) - Row 1
            subplot(2, 2, g, 'Parent', fig_heatmaps);
            if n_sig > 0
                imagesc(bin_centers, 1:n_sig, sig_rate_mat);
                colormap(gca, 'jet');
            else
                % Empty plot if no significant clusters
                axis off;
                text(0.5, 0.5, 'No significant clusters', 'Units', 'normalized', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
            xlabel('Position (cm)');
            if g == 1
                ylabel('Significant Clusters');
            else
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'YTick', []);
            
            % Set x-axis limits
            if strcmp(group, 'short')
                xlim([60, 120]);
            else
                xlim([0, 120]);
            end
            title(analyzer.group_labels{g});
            
            % Non-significant clusters (grayscale) - Row 2
            subplot(2, 2, g+2, 'Parent', fig_heatmaps);
            if n_nonsig > 0
                % Show the actual spatial profile in grayscale
                imagesc(bin_centers, 1:n_nonsig, nonsig_rate_mat);
                colormap(gca, gray(256));
            else
                % Empty plot if no non-significant clusters
                axis off;
                text(0.5, 0.5, 'No non-significant clusters', 'Units', 'normalized', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
            xlabel('Position (cm)');
            if g == 1
                ylabel('Non-Significant Clusters');
            else
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'YTick', []);
            
            % Set x-axis limits
            if strcmp(group, 'short')
                xlim([60, 120]);
            else
                xlim([0, 120]);
            end
        end
        
        % Adjust subplot positions to make right column half width and reduce gap
        all_axes = findall(fig_heatmaps, 'Type', 'axes');
        for ax_idx = 1:length(all_axes)
            ax = all_axes(ax_idx);
            pos = get(ax, 'Position');
            
            % Determine if this is a left or right subplot
            if pos(1) > 0.5  % Right column
                pos(3) = pos(3) * 0.5;  % Make width half
                pos(1) = 0.55;  % Reduce gap between columns (was pos(1) + pos(3) * 0.5)
            else  % Left column
                pos(3) = pos(3) * 1.05;  % Slightly increase left column width
            end
            
            pos(4) = pos(4) * 0.82;  % Reduce height
            pos(2) = max(0.08, pos(2) - 0.03);  % Lower position
            set(ax, 'Position', pos);
        end
        
        % Add figure title
        FigureTitle(fig_heatmaps, sprintf('Heatmaps - %s', probe_ids{pid}));
        ctl.figs.save_fig_to_join();
        
        % --- Time-to-Goal Heatmaps (Sorted by Long Peak Position) ---
        if save_figs
            fprintf('  Creating Time-to-Goal Heatmaps...\n');
            fig_ttg_heatmaps = plot_time_to_goal_heatmaps(all_ttg_rates, all_ttg_fields, ...
                ttg_norm_bin_centers, analyzer.cluster_ids, probe_ids{pid}, ctl);
            ctl.figs.save_fig_to_join();
        end
        
        % Plot spatially tuned clusters
        if ~isempty(analyzer.tuning_stats)
            fprintf('  Creating spatially tuned clusters plot...\n');
            
            % Identify clusters tuned in at least one condition
            spatially_tuned_any = [];
            
            for c = 1:length(analyzer.cluster_ids)
                cluster_field = sprintf('cluster_%d', analyzer.cluster_ids(c));
                if isfield(analyzer.tuning_stats, cluster_field)
                    stats = analyzer.tuning_stats.(cluster_field);
                    if (isfield(stats, 'long') && stats.long.isSpatiallyTuned) || ...
                       (isfield(stats, 'short') && stats.short.isSpatiallyTuned)
                        spatially_tuned_any = [spatially_tuned_any; c]; %#ok<AGROW>
                    end
                end
            end
            
            n_tuned = length(spatially_tuned_any);
            fprintf('    Found %d clusters spatially tuned in at least one condition\n', n_tuned);
            
            % Create 3D visualization and get fit results
            [fig_3d, fit_results] = plot_spatial_tuning_3d(analyzer.tuning_stats, ...
                analyzer.bin_config, analyzer.firing_rates, ...
                analyzer.cluster_ids, probe_ids{pid}, ctl);
            
            if ~isempty(fig_3d)
                ctl.figs.save_fig_to_join();
            end
            

        end
    end

    
    % Save TTG data to a separate cache file for pooled analysis
    if ~isempty(all_ttg_rates.long) || ~isempty(all_ttg_rates.short)
        ttg_cache_filename = sprintf('%s_ttg_cache.mat', probe_ids{pid});
        ttg_cache_filepath = fullfile(ctl.figs.curr_dir, ttg_cache_filename);
        fprintf('  Saving TTG data to cache: %s\n', ttg_cache_filename);
        save(ttg_cache_filepath, 'all_ttg_rates', 'ttg_norm_bin_centers', 'analyzer.cluster_ids', '-v7.3');
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
    clear analyzer data clusters sessions
    clear all_ttg_rates ttg_norm_bin_centers
    
    fprintf('  Memory cleaned up after probe %d/%d\n', pid, length(probe_ids));
end

fprintf('\n=== All probes processed ===\n');

% Restore default figure visibility
set(groot, 'DefaultFigureVisible', 'on');
