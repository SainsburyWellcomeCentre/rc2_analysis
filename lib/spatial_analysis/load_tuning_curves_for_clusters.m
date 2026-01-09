function [tuning_curves, accel_tuning_curves] = load_tuning_curves_for_clusters(...
    data, cluster_ids, plot_velocity, plot_acceleration, ...
    trial_group_label_velocity, trial_group_label_accel)
% LOAD_TUNING_CURVES_FOR_CLUSTERS Load velocity and acceleration tuning curves
%
% SYNTAX:
%   [tuning_curves, accel_tuning_curves] = load_tuning_curves_for_clusters(...
%       data, cluster_ids, plot_velocity, plot_acceleration, ...
%       trial_group_label_velocity, trial_group_label_accel)
%
% INPUTS:
%   data                       - Data object with load methods
%   cluster_ids                - Array of cluster IDs
%   plot_velocity              - Boolean; load velocity curves if true
%   plot_acceleration          - Boolean; load acceleration curves if true
%   trial_group_label_velocity - Label for velocity tuning group
%   trial_group_label_accel    - Label for acceleration tuning group
%
% OUTPUTS:
%   tuning_curves              - Cell array of velocity tuning curves
%   accel_tuning_curves        - Cell array of acceleration tuning curves
%
% DESCRIPTION:
%   Loads tuning curves for all clusters. Returns empty arrays for curves
%   if corresponding plot flag is false. Handles errors gracefully with warnings.

    n_clusters = length(cluster_ids);
    
    % Initialize outputs
    tuning_curves = cell(n_clusters, 1);
    accel_tuning_curves = cell(n_clusters, 1);
    
    % Load velocity tuning curves
    if plot_velocity
        for c = 1:n_clusters
            try
                tuning_curves{c} = data.load_tuning_curves(...
                    cluster_ids(c), trial_group_label_velocity);
            catch ME
                fprintf('    Warning: Could not load velocity tuning curve for cluster %d\n', ...
                        cluster_ids(c));
                fprintf('    Error: %s\n', ME.message);
                tuning_curves{c} = [];
            end
        end
    end
    
    % Load acceleration tuning curves
    if plot_acceleration
        for c = 1:n_clusters
            try
                % load_tuning_curves_acceleration returns {all, acc, dec}
                % We want the 'all' table (index 1) for combined acceleration
                accel_data = data.load_tuning_curves_acceleration(...
                    cluster_ids(c), trial_group_label_accel);
                if ~isempty(accel_data) && length(accel_data) >= 1
                    accel_tuning_curves{c} = accel_data{1};  % Extract 'all' table
                else
                    accel_tuning_curves{c} = [];
                end
            catch ME
                fprintf('    Warning: Could not load acceleration tuning curve for cluster %d\n', ...
                        cluster_ids(c));
                fprintf('    Error: %s\n', ME.message);
                accel_tuning_curves{c} = [];
            end
        end
    end
end
