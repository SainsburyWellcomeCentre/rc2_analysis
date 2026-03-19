function [all_ttg_rates, ttg_norm_bin_centers, all_ttg_fields] = create_single_cluster_figures(...
    plot_single_cluster_fig, analyzer, clusters, data, ...
    plot_velocity_tuning, plot_acceleration_tuning, ...
    trial_group_label_for_tuning, trial_group_label_for_accel_tuning, ...
    bin_size_cm, gauss_sigma_cm, probe_id, save_figs, ctl, re_run_analysis)
% CREATE_SINGLE_CLUSTER_FIGURES Create individual spatial profile figures for each cluster
%
% Creates a combined 6-row figure for each cluster showing spatial firing rate profiles,
% velocity/acceleration tuning curves (combined and split by trial type), 2D rate maps,
% and x-normalized comparisons between trial groups.
%
% INPUTS:
%   plot_single_cluster_fig           - Flag to enable/disable plotting
%   analyzer                          - SpatialTuningAnalyzer object with analysis results
%   clusters                          - Array of cluster objects
%   data                              - Data object for loading tuning curves
%   plot_velocity_tuning              - Flag to include velocity tuning
%   plot_acceleration_tuning          - Flag to include acceleration tuning
%   trial_group_label_for_tuning      - Label for velocity tuning trials
%   trial_group_label_for_accel_tuning - Label for acceleration tuning trials
%   bin_size_cm                       - Spatial bin size in cm
%   gauss_sigma_cm                    - Gaussian smoothing sigma in cm
%   probe_id                          - String identifier for the probe
%   save_figs                         - Flag to save figures
%   ctl                               - RC2Analysis controller object
%
% OUTPUTS:
%   all_ttg_rates                     - Struct with TTG rates for long/short trials
%   ttg_norm_bin_centers              - TTG bin centers for plotting
%   all_ttg_fields                    - Cell array of TTG field information per cluster
%
% See also: plot_cluster_spatial_profile, prepare_cluster_data_for_plotting,
%           load_tuning_curves_for_clusters

    % Initialize outputs
    all_ttg_rates = struct('long', [], 'short', []);
    ttg_norm_bin_centers = [];
    all_ttg_fields = cell(length(analyzer.cluster_ids), 1);
    
    % Early return if plotting is disabled
    if ~plot_single_cluster_fig
        return;
    end
    
    fprintf('  Creating individual cluster plots (%d clusters)...\n', length(analyzer.cluster_ids));
    
    % Define colors for the two groups
    group_colors = struct('long', [0 0 0.8], 'short', [0.8 0 0]); % Blue for long, red for short
    
    % Load tuning curves for all clusters using helper function
    [tuning_curves, accel_tuning_curves] = load_tuning_curves_for_clusters(...
        data, analyzer.cluster_ids, plot_velocity_tuning, plot_acceleration_tuning, ...
        trial_group_label_for_tuning, trial_group_label_for_accel_tuning);
    
    % Compute split tuning curves with caching
    split_cache_filename = sprintf('%s_split_tuning_cache.mat', probe_id);
    split_cache_filepath = fullfile(ctl.figs.curr_dir, split_cache_filename);
    
    if exist(split_cache_filepath, 'file') && ~re_run_analysis
        fprintf('  Loading split tuning curves from cache...\n');
        load(split_cache_filepath, 'tuning_curves_long', 'tuning_curves_short', ...
             'accel_tuning_curves_long', 'accel_tuning_curves_short');
    else
        fprintf('  Computing split tuning curves (this may take a while)...\n');
        [tuning_curves_long, tuning_curves_short, accel_tuning_curves_long, accel_tuning_curves_short] = ...
            compute_split_tuning_curves(data, clusters, analyzer.trial_groups, ...
                                         plot_velocity_tuning, plot_acceleration_tuning, ...
                                         tuning_curves, accel_tuning_curves);
        % Save to cache
        fprintf('  Saving split tuning curves to cache...\n');
        save(split_cache_filepath, 'tuning_curves_long', 'tuning_curves_short', ...
             'accel_tuning_curves_long', 'accel_tuning_curves_short', '-v7.3');
    end
    
    for c = 1:length(analyzer.cluster_ids)
        % Show progress every 10 clusters to avoid cursor stealing
        if mod(c, 10) == 1 || c == length(analyzer.cluster_ids)
            fprintf('    Processing cluster %d/%d...\n', c, length(analyzer.cluster_ids));
        end
        
        cluster = clusters(c);
        
        % Prepare all cluster data for plotting using helper function
        [rate_data, rate_maps_2d, cluster_stats, ttg_data, ttg_stats] = ...
            prepare_cluster_data_for_plotting(cluster, c, analyzer, ...
                                               analyzer.trial_groups, bin_size_cm, gauss_sigma_cm, ...
                                               plot_velocity_tuning, plot_acceleration_tuning);
        
        % Get TTG bin centers for later use
        cluster_field = sprintf('cluster_%d', cluster.id);
        ttg_norm_bin_centers = analyzer.ttg_analysis.(cluster_field).ttg_norm_bin_centers;
        
        % Store median TTG rates for later heatmap plotting
        for g = 1:2
            group = analyzer.group_names{g};
            Q2_smooth = quantile(ttg_data.(group).rate_per_trial_smooth, 0.50, 1);
            all_ttg_rates.(group)(c, :) = Q2_smooth;
        end
        
        % Extract bin centers for plotting
        bin_centers_for_plot = struct();
        bin_centers_for_plot.long = analyzer.bin_config.long.plot_centers;
        bin_centers_for_plot.short = analyzer.bin_config.short.plot_centers;

        % Create combined figure using helper function
        try
            [fig_cluster, ttg_fields] = plot_cluster_spatial_profile(cluster.id, bin_centers_for_plot, ...
                                                       rate_data, analyzer.group_names, analyzer.group_labels, ...
                                                       group_colors, probe_id, cluster_stats, tuning_curves{c}, accel_tuning_curves{c}, ...
                                                       rate_maps_2d, analyzer.distribution_comparisons{c}, ttg_data, ttg_stats, ...
                                                       tuning_curves_long{c}, tuning_curves_short{c}, ...
                                                       accel_tuning_curves_long{c}, accel_tuning_curves_short{c});
            
            % Store TTG field information for this cluster
            all_ttg_fields{c} = ttg_fields;
        catch ME
            warning('Error creating figure for cluster %d: %s', cluster.id, ME.message);
            fprintf('Error occurred at: %s\n', ME.stack(1).name);
            fprintf('Line: %d\n', ME.stack(1).line);
            if length(ME.stack) > 1
                fprintf('Called from: %s (line %d)\n', ME.stack(2).name, ME.stack(2).line);
            end
            continue;
        end
        
        % Verify figure was created successfully
        if isempty(fig_cluster) || ~isvalid(fig_cluster)
            warning('Failed to create figure for cluster %d, skipping...', cluster.id);
            continue;
        end
        
        % Keep figure invisible to prevent focus stealing
        set(fig_cluster, 'Visible', 'off');
        
        % Save using the figure management system
        if save_figs
            % Convert to a4figure format for PDF joining (portrait for 6-row layout)
            set(fig_cluster, 'PaperOrientation', 'portrait');
            set(fig_cluster, 'PaperUnits', 'normalized');
            set(fig_cluster, 'PaperPosition', [0 0 1 1]);
            set(0, 'CurrentFigure', fig_cluster);  % Set as current figure for save_fig_to_join
            ctl.figs.save_fig_to_join();
        else
            close(fig_cluster);
        end
    end
    % Save tuning summary text file
    if save_figs
        fprintf('  Saving tuning summary...\n');
        save_tuning_summary(tuning_curves, accel_tuning_curves, ...
                    tuning_curves_long, tuning_curves_short, ...
                    accel_tuning_curves_long, accel_tuning_curves_short, ...
                    analyzer.cluster_ids, ctl.figs.curr_dir, probe_id);
    end
end