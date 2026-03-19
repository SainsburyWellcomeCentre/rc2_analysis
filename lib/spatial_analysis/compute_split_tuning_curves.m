function [tuning_curves_long, tuning_curves_short, accel_tuning_curves_long, accel_tuning_curves_short] = ...
    compute_split_tuning_curves(data, clusters, trial_groups, plot_velocity, plot_acceleration, tuning_curves_combined, accel_tuning_curves_combined)
% COMPUTE_SPLIT_TUNING_CURVES Compute velocity and acceleration tuning curves split by trial type
%
% SYNTAX:
%   [tuning_curves_long, tuning_curves_short, accel_tuning_curves_long, accel_tuning_curves_short] = ...
%       compute_split_tuning_curves(data, clusters, trial_groups, plot_velocity, plot_acceleration, tuning_curves_combined, accel_tuning_curves_combined)
%
% INPUTS:
%   data                        - Data object with probe_id
%   clusters                    - Array of cluster objects
%   trial_groups                - Struct with 'long' and 'short' trial arrays
%   plot_velocity               - Boolean; compute velocity curves if true
%   plot_acceleration           - Boolean; compute acceleration curves if true
%   tuning_curves_combined      - Cell array of combined velocity tuning curves (to extract best model)
%   accel_tuning_curves_combined - Cell array of combined acceleration tuning curves (to extract best model)
%
% OUTPUTS:
%   tuning_curves_long          - Cell array of velocity tuning curves for long trials
%   tuning_curves_short         - Cell array of velocity tuning curves for short trials
%   accel_tuning_curves_long    - Cell array of acceleration tuning curves for long trials
%   accel_tuning_curves_short   - Cell array of acceleration tuning curves for short trials
%
% DESCRIPTION:
%   Computes tuning curves separately for long and short trials WITHOUT shuffling.
%   Reuses the best model from the combined data for efficiency.

    n_clusters = length(clusters);
    
    % Initialize outputs
    tuning_curves_long = cell(n_clusters, 1);
    tuning_curves_short = cell(n_clusters, 1);
    accel_tuning_curves_long = cell(n_clusters, 1);
    accel_tuning_curves_short = cell(n_clusters, 1);
    
    % Get probe_id from the data object
    probe_id = data.probe_id;
    
    % Extract trial objects from trial group structs
    % trial_groups.long and trial_groups.short are struct arrays with .trial field
    long_trial_objects = arrayfun(@(x) x.trial, trial_groups.long, 'UniformOutput', false);
    short_trial_objects = arrayfun(@(x) x.trial, trial_groups.short, 'UniformOutput', false);
    
    % Process velocity tuning curves
    if plot_velocity
        % Create TuningTable instances for each trial group
        tt_long = TuningTable(probe_id);
        tt_long.add_trials(long_trial_objects, repmat({'long'}, length(long_trial_objects), 1));
        
        tt_short = TuningTable(probe_id);
        tt_short.add_trials(short_trial_objects, repmat({'short'}, length(short_trial_objects), 1));
        
        % Compute tuning curves for each cluster (NO SHUFFLING)
        for c = 1:n_clusters
            try
                % Bin the data for long trials
                [tuning_long, ~] = tt_long.vtc.fr_curve_count(clusters(c));
                
                % Create output struct with bin info and tuning data
                tuning_curves_long{c} = struct();
                tuning_curves_long{c}.bin_edges = tt_long.velocity_bins.bin_edges;
                tuning_curves_long{c}.bin_centers = tt_long.velocity_bins.bin_centers;
                tuning_curves_long{c}.tuning = tuning_long;  % n_bins x n_trials matrix

                % Add model selection and shuffling
                mst = ModelSelectionTuning(tuning_long, tt_long.velocity_bins.bin_centers);
                tuning_curves_long{c}.shuffled = mst.get_summary();
                
            catch ME
                fprintf('    Warning: Could not compute velocity tuning (long) for cluster %d\n', ...
                        clusters(c).id);
                fprintf('    Error: %s\n', ME.message);
                tuning_curves_long{c} = [];
            end
            
            try
                % Bin the data for short trials
                [tuning_short, ~] = tt_short.vtc.fr_curve_count(clusters(c));
                
                % Create output struct with bin info and tuning data
                tuning_curves_short{c} = struct();
                tuning_curves_short{c}.bin_edges = tt_short.velocity_bins.bin_edges;
                tuning_curves_short{c}.bin_centers = tt_short.velocity_bins.bin_centers;
                tuning_curves_short{c}.tuning = tuning_short;  % n_bins x n_trials matrix

                % Add model selection and shuffling
                mst = ModelSelectionTuning(tuning_short, tt_short.velocity_bins.bin_centers);
                tuning_curves_short{c}.shuffled = mst.get_summary();
                
            catch ME
                fprintf('    Warning: Could not compute velocity tuning (short) for cluster %d\n', ...
                        clusters(c).id);
                fprintf('    Error: %s\n', ME.message);
                tuning_curves_short{c} = [];
            end
        end
    end
    
    % Process acceleration tuning curves
    if plot_acceleration
        % Create TuningTableAcc instances for each trial group
        tta_long = TuningTableAcc(probe_id);
        tta_long.add_trials(long_trial_objects, repmat({'long'}, length(long_trial_objects), 1));
        
        tta_short = TuningTableAcc(probe_id);
        tta_short.add_trials(short_trial_objects, repmat({'short'}, length(short_trial_objects), 1));
        
        % Compute acceleration tuning curves for each cluster (NO SHUFFLING)
        for c = 1:n_clusters
            try
                % Bin the data for long trials
                [tuning_long, ~] = tta_long.atc.fr_curve_count(clusters(c));
                
                % Create output struct with bin info and tuning data
                accel_tuning_curves_long{c} = struct();
                accel_tuning_curves_long{c}.bin_edges = tta_long.acceleration_bins.bin_edges;
                accel_tuning_curves_long{c}.bin_centers = tta_long.acceleration_bins.bin_centers;
                accel_tuning_curves_long{c}.tuning = tuning_long;  % n_bins x n_trials matrix

                % Add model selection and shuffling
                mst = ModelSelectionTuning(tuning_long, tta_long.acceleration_bins.bin_centers);
                accel_tuning_curves_long{c}.shuffled = mst.get_summary();
                
            catch ME
                fprintf('    Warning: Could not compute acceleration tuning (long) for cluster %d\n', ...
                        clusters(c).id);
                fprintf('    Error: %s\n', ME.message);
                accel_tuning_curves_long{c} = [];
            end
            
            try
                % Bin the data for short trials
                [tuning_short, ~] = tta_short.atc.fr_curve_count(clusters(c));
                
                % Create output struct with bin info and tuning data
                accel_tuning_curves_short{c} = struct();
                accel_tuning_curves_short{c}.bin_edges = tta_short.acceleration_bins.bin_edges;
                accel_tuning_curves_short{c}.bin_centers = tta_short.acceleration_bins.bin_centers;
                accel_tuning_curves_short{c}.tuning = tuning_short;  % n_bins x n_trials matrix

                % Add model selection and shuffling
                mst = ModelSelectionTuning(tuning_short, tta_short.acceleration_bins.bin_centers);
                accel_tuning_curves_short{c}.shuffled = mst.get_summary();
                
            catch ME
                fprintf('    Warning: Could not compute acceleration tuning (short) for cluster %d\n', ...
                        clusters(c).id);
                fprintf('    Error: %s\n', ME.message);
                accel_tuning_curves_short{c} = [];
            end
        end
    end
end
