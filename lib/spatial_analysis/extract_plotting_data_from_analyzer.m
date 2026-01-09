function plot_data = extract_plotting_data_from_analyzer(analyzer)
% EXTRACT_PLOTTING_DATA_FROM_ANALYZER Extract data from analyzer for plotting
%
% SYNTAX:
%   plot_data = extract_plotting_data_from_analyzer(analyzer)
%
% INPUTS:
%   analyzer   - SpatialTuningAnalyzer instance with completed analysis
%
% OUTPUTS:
%   plot_data  - Struct containing all data needed for plotting, organized by group
%
% DESCRIPTION:
%   Extracts and organizes data from the analyzer into a convenient format
%   for plotting functions. This consolidates repetitive data extraction code.

    plot_data = struct();
    
    % Basic info
    plot_data.cluster_ids = analyzer.cluster_ids;
    plot_data.n_clusters = length(analyzer.cluster_ids);
    plot_data.group_names = analyzer.group_names;
    plot_data.group_labels = analyzer.group_labels;
    plot_data.trial_groups = analyzer.trial_groups;
    
    % Statistics
    plot_data.spatial_tuning_stats = analyzer.tuning_stats;
    plot_data.dist_comparison_results = analyzer.distribution_comparisons;
    
    % Initialize per-group data structures
    for g = 1:2
        group = plot_data.group_names{g};
        
        % Bin configuration
        plot_data.bin_centers.(group) = analyzer.bin_config.(group).plot_centers;
        plot_data.edges.(group) = analyzer.bin_config.(group).edges;
        
        % Firing rate data
        plot_data.Q1_rate_smooth.(group) = analyzer.firing_rates.(group).Q1_smooth;
        plot_data.Q2_rate_smooth.(group) = analyzer.firing_rates.(group).Q2_smooth;
        plot_data.Q3_rate_smooth.(group) = analyzer.firing_rates.(group).Q3_smooth;
        plot_data.avg_velocity.(group) = analyzer.firing_rates.(group).avg_velocity;
        plot_data.occupancy.(group) = analyzer.firing_rates.(group).occupancy;
        plot_data.rate_per_trial.(group) = analyzer.firing_rates.(group).rate_per_trial;
        plot_data.rate_per_trial_smooth.(group) = analyzer.firing_rates.(group).rate_per_trial_smooth;
        plot_data.spike_positions.(group) = analyzer.firing_rates.(group).spike_positions;
        plot_data.global_trial_indices.(group) = analyzer.firing_rates.(group).global_trial_indices;
    end
end
