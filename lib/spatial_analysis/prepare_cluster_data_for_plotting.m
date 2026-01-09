function [rate_data, rate_maps_2d, cluster_stats, ttg_data, ttg_stats] = ...
    prepare_cluster_data_for_plotting(cluster, c, analyzer, plot_data, ...
                                       trial_groups, bin_size_cm, gauss_sigma_cm, ...
                                       plot_velocity, plot_acceleration)
% PREPARE_CLUSTER_DATA_FOR_PLOTTING Prepare all data needed for cluster plotting
%
% SYNTAX:
%   [rate_data, rate_maps_2d, cluster_stats, ttg_data, ttg_stats] = ...
%       prepare_cluster_data_for_plotting(cluster, c, analyzer, plot_data, ...
%                                          trial_groups, bin_size_cm, gauss_sigma_cm, ...
%                                          plot_velocity, plot_acceleration)
%
% INPUTS:
%   cluster             - Cluster object
%   c                   - Cluster index
%   analyzer            - SpatialTuningAnalyzer instance
%   plot_data           - Struct from extract_plotting_data_from_analyzer
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

    group_names = plot_data.group_names;
    
    % --- Compute downsampled long data ---
    [rate_long_ds_smooth, Q1_long_ds, Q2_long_ds, Q3_long_ds, bin_centers_long_ds, edges_long_ds] = ...
        compute_downsampled_long_rates(trial_groups, cluster, bin_size_cm, gauss_sigma_cm);
    
    % --- Prepare rate data structure ---
    rate_data = struct();
    for g = 1:2
        group = group_names{g};
        rate_data.(group).Q1_smooth = plot_data.Q1_rate_smooth.(group)(c, :);
        rate_data.(group).Q2_smooth = plot_data.Q2_rate_smooth.(group)(c, :);
        rate_data.(group).Q3_smooth = plot_data.Q3_rate_smooth.(group)(c, :);
        rate_data.(group).rate_per_trial = plot_data.rate_per_trial.(group){c};
        rate_data.(group).rate_per_trial_smooth = plot_data.rate_per_trial_smooth.(group){c};
        rate_data.(group).spike_positions = plot_data.spike_positions.(group){c};
        rate_data.(group).global_trial_indices = plot_data.global_trial_indices.(group){c};
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
    if ~isempty(plot_data.spatial_tuning_stats)
        cluster_field = sprintf('cluster_%d', cluster.id);
        if isfield(plot_data.spatial_tuning_stats, cluster_field)
            cluster_stats = plot_data.spatial_tuning_stats.(cluster_field);
        end
    end
    
    % --- Compute 2D rate maps ---
    rate_maps_2d = struct();
    
    % Compute velocity and acceleration maps in parallel
    if plot_velocity
        n_vel_bins = 20;
        sigma_vel = 5;
        
        % Pre-allocate outputs
        vel_maps = cell(1, 2);
        vel_edges = cell(1, 2);
        vel_counts = cell(1, 2);
        group_names_local = {'long', 'short'};
        
        % Parallel computation for both groups
        parfor g = 1:2
            [vel_maps{g}, vel_edges{g}, vel_counts{g}] = ...
                analyzer.compute_2d_position_velocity_tuning(cluster, group_names_local{g}, n_vel_bins, sigma_vel);
        end
        
        % Unpack results
        rate_maps_2d.vel_long = vel_maps{1};
        rate_maps_2d.vel_short = vel_maps{2};
        vel_bin_edges_long = vel_edges{1};
        vel_bin_edges_short = vel_edges{2};
        vel_trial_counts_long = vel_counts{1};
        vel_trial_counts_short = vel_counts{2};
        
        % Store metadata
        rate_maps_2d.vel_bin_edges_long = vel_bin_edges_long;
        rate_maps_2d.vel_bin_edges_short = vel_bin_edges_short;
        rate_maps_2d.vel_bin_centers_long = (vel_bin_edges_long(1:end-1) + vel_bin_edges_long(2:end)) / 2;
        rate_maps_2d.vel_bin_centers_short = (vel_bin_edges_short(1:end-1) + vel_bin_edges_short(2:end)) / 2;
        rate_maps_2d.vel_trial_counts_long = vel_trial_counts_long;
        rate_maps_2d.vel_trial_counts_short = vel_trial_counts_short;
        rate_maps_2d.vel_n_trials_long = length(trial_groups.long);
        rate_maps_2d.vel_n_trials_short = length(trial_groups.short);
    end
    
    if plot_acceleration
        n_accel_bins = 20;
        sigma_accel = 10;
        
        % Pre-allocate outputs
        accel_maps = cell(1, 2);
        accel_edges = cell(1, 2);
        accel_counts = cell(1, 2);
        group_names_local = {'long', 'short'};
        
        % Parallel computation for both groups
        parfor g = 1:2
            [accel_maps{g}, accel_edges{g}, accel_counts{g}] = ...
                analyzer.compute_2d_position_acceleration_tuning(cluster, group_names_local{g}, n_accel_bins, sigma_accel);
        end
        
        % Unpack results
        rate_maps_2d.accel_long = accel_maps{1};
        rate_maps_2d.accel_short = accel_maps{2};
        accel_bin_edges_long = accel_edges{1};
        accel_bin_edges_short = accel_edges{2};
        accel_trial_counts_long = accel_counts{1};
        accel_trial_counts_short = accel_counts{2};
        
        % Store metadata
        rate_maps_2d.accel_bin_edges_long = accel_bin_edges_long;
        rate_maps_2d.accel_bin_edges_short = accel_bin_edges_short;
        rate_maps_2d.accel_bin_centers_long = (accel_bin_edges_long(1:end-1) + accel_bin_edges_long(2:end)) / 2;
        rate_maps_2d.accel_bin_centers_short = (accel_bin_edges_short(1:end-1) + accel_bin_edges_short(2:end)) / 2;
        rate_maps_2d.accel_trial_counts_long = accel_trial_counts_long;
        rate_maps_2d.accel_trial_counts_short = accel_trial_counts_short;
        rate_maps_2d.accel_n_trials_long = length(trial_groups.long);
        rate_maps_2d.accel_n_trials_short = length(trial_groups.short);
    end
    
    % Store position bins
    rate_maps_2d.pos_bins_long = analyzer.bin_config.long.centers;
    rate_maps_2d.pos_bins_short = analyzer.bin_config.short.centers;
    rate_maps_2d.pos_bin_edges_long = analyzer.bin_config.long.edges;
    rate_maps_2d.pos_bin_edges_short = analyzer.bin_config.short.edges;
    
    % --- Get TTG data ---
    cluster_field = sprintf('cluster_%d', cluster.id);
    ttg_data = analyzer.ttg_analysis.(cluster_field).ttg_data;
    ttg_stats = analyzer.ttg_analysis.(cluster_field).ttg_stats;
    
    % --- Compute TTG × kinematic maps ---
    rate_maps_2d = compute_ttg_kinematic_maps(ttg_data, trial_groups, rate_maps_2d, ...
                                               plot_velocity, plot_acceleration);
end
