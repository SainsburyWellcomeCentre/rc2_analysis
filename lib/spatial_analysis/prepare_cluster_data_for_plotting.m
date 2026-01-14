function [rate_data, rate_maps_2d, cluster_stats, ttg_data, ttg_stats] = ...
    prepare_cluster_data_for_plotting(cluster, c, analyzer, ...
                                       trial_groups, bin_size_cm, gauss_sigma_cm, ...
                                       plot_velocity, plot_acceleration)
% PREPARE_CLUSTER_DATA_FOR_PLOTTING Prepare all data needed for cluster plotting
%
% SYNTAX:
%   [rate_data, rate_maps_2d, cluster_stats, ttg_data, ttg_stats] = ...
%       prepare_cluster_data_for_plotting(cluster, c, analyzer, ...
%                                          trial_groups, bin_size_cm, gauss_sigma_cm, ...
%                                          plot_velocity, plot_acceleration)
%
% INPUTS:
%   cluster             - Cluster object
%   c                   - Cluster index
%   analyzer            - SpatialTuningAnalyzer instance
%   trial_groups        - Struct with trial arrays
%   bin_size_cm         - Bin size in cm
%   gauss_sigma_cm      - Gaussian smoothing sigma in cm
%   plot_velocity       - Boolean flag
%   plot_acceleration   - Boolean flag
%
% OUTPUTS:
%   rate_data           - Struct with rate data for this cluster
%   rate_maps_2d        - Struct with 2D rate maps
%   cluster_stats       - Statistics for this cluster
%   ttg_data            - Time-to-goal data
%   ttg_stats           - Time-to-goal statistics

    group_names = analyzer.group_names;
    
    % --- Get downsampled long data from cache ---
    downsampled = analyzer.downsampled_long_rates{c};
    rate_long_ds_smooth = downsampled.rate_per_trial_smooth;
    Q1_long_ds = downsampled.Q1_smooth;
    Q2_long_ds = downsampled.Q2_smooth;
    Q3_long_ds = downsampled.Q3_smooth;
    bin_centers_long_ds = downsampled.bin_centers;
    edges_long_ds = downsampled.edges;
    
    % --- Prepare rate data structure ---
    rate_data = struct();
    for g = 1:2
        group = group_names{g};
        rate_data.(group).Q1_smooth = analyzer.firing_rates.(group).Q1_smooth(c, :);
        rate_data.(group).Q2_smooth = analyzer.firing_rates.(group).Q2_smooth(c, :);
        rate_data.(group).Q3_smooth = analyzer.firing_rates.(group).Q3_smooth(c, :);
        rate_data.(group).rate_per_trial = analyzer.firing_rates.(group).rate_per_trial{c};
        rate_data.(group).rate_per_trial_smooth = analyzer.firing_rates.(group).rate_per_trial_smooth{c};
        rate_data.(group).spike_positions = analyzer.firing_rates.(group).spike_positions{c};
        rate_data.(group).global_trial_indices = analyzer.firing_rates.(group).global_trial_indices{c};
    end
    
    % Add downsampled long data
    rate_data.long_downsampled = struct();
    rate_data.long_downsampled.Q1_smooth = Q1_long_ds;
    rate_data.long_downsampled.Q2_smooth = Q2_long_ds;
    rate_data.long_downsampled.Q3_smooth = Q3_long_ds;
    rate_data.long_downsampled.rate_per_trial_smooth = rate_long_ds_smooth;
    rate_data.long_downsampled.bin_centers = bin_centers_long_ds;
    rate_data.long_downsampled.edges = edges_long_ds;
    
    % --- Get cluster statistics ---
    cluster_stats = [];
    if ~isempty(analyzer.tuning_stats)
        cluster_field = sprintf('cluster_%d', cluster.id);
        if isfield(analyzer.tuning_stats, cluster_field)
            cluster_stats = analyzer.tuning_stats.(cluster_field);
        end
    end
    
    % --- Get pre-computed 2D rate maps from cache ---
    rate_maps_2d = analyzer.rate_maps_2d{c};
    
    % --- Get TTG data ---
    cluster_field = sprintf('cluster_%d', cluster.id);
    ttg_data = analyzer.ttg_analysis.(cluster_field).ttg_data;
    ttg_stats = analyzer.ttg_analysis.(cluster_field).ttg_stats;
end
