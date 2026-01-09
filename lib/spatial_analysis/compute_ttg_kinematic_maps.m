function rate_maps_2d = compute_ttg_kinematic_maps(ttg_data, trial_groups, rate_maps_2d, ...
                                                    plot_velocity, plot_acceleration)
% COMPUTE_TTG_KINEMATIC_MAPS Compute 2D maps of TTG × velocity/acceleration
%
% SYNTAX:
%   rate_maps_2d = compute_ttg_kinematic_maps(ttg_data, trial_groups, rate_maps_2d, ...
%                                              plot_velocity, plot_acceleration)
%
% INPUTS:
%   ttg_data          - Struct with TTG analysis results for 'long' and 'short'
%   trial_groups      - Struct with trial arrays for 'long' and 'short'
%   rate_maps_2d      - Struct to append results to (must contain bin edges)
%   plot_velocity     - Boolean; compute velocity maps if true
%   plot_acceleration - Boolean; compute acceleration maps if true
%
% OUTPUTS:
%   rate_maps_2d      - Updated struct with TTG × kinematic rate maps
%
% DESCRIPTION:
%   Computes 2D rate maps showing how firing rate varies with normalized
%   time-to-goal and either velocity or acceleration. Maps are computed by
%   assigning each time point to bins in both dimensions and using the
%   smoothed TTG firing rate for that trial.

    group_names = {'long', 'short'};
    
    % --- Compute TTG × Velocity maps ---
    if plot_velocity
        for g = 1:2
            group = group_names{g};
            
            % Get velocity bin edges (pre-computed)
            if strcmp(group, 'long')
                vel_edges = rate_maps_2d.vel_bin_edges_long;
            else
                vel_edges = rate_maps_2d.vel_bin_edges_short;
            end
            n_vel_bins = length(vel_edges) - 1;
            
            % Get TTG data
            ttg_edges = ttg_data.(group).bin_edges_norm;
            n_ttg_bins = length(ttg_edges) - 1;
            rate_per_trial_smooth = ttg_data.(group).rate_per_trial_smooth;
            
            trials = trial_groups.(group);
            n_trials = length(trials);
            
            % Compute map
            [rate_map, trial_counts, ttg_bin_occupancy] = compute_ttg_kinematic_map_for_group(...
                trials, ttg_edges, n_ttg_bins, vel_edges, n_vel_bins, ...
                rate_per_trial_smooth, 'velocity');
            
            % Store results
            rate_maps_2d.(['ttg_vel_' group]) = rate_map;
            rate_maps_2d.(['ttg_vel_trial_counts_' group]) = trial_counts;
            rate_maps_2d.(['ttg_vel_n_trials_' group]) = n_trials;
            rate_maps_2d.(['ttg_vel_ttg_occupancy_' group]) = ttg_bin_occupancy;
            rate_maps_2d.(['ttg_bin_edges_' group]) = ttg_edges;
        end
    end
    
    % --- Compute TTG × Acceleration maps ---
    if plot_acceleration
        for g = 1:2
            group = group_names{g};
            
            % Get acceleration bin edges (pre-computed)
            if strcmp(group, 'long')
                accel_edges = rate_maps_2d.accel_bin_edges_long;
            else
                accel_edges = rate_maps_2d.accel_bin_edges_short;
            end
            n_accel_bins = length(accel_edges) - 1;
            
            % Get TTG data
            ttg_edges = ttg_data.(group).bin_edges_norm;
            n_ttg_bins = length(ttg_edges) - 1;
            rate_per_trial_smooth = ttg_data.(group).rate_per_trial_smooth;
            
            trials = trial_groups.(group);
            n_trials = length(trials);
            
            % Compute map
            [rate_map, trial_counts, ttg_bin_occupancy] = compute_ttg_kinematic_map_for_group(...
                trials, ttg_edges, n_ttg_bins, accel_edges, n_accel_bins, ...
                rate_per_trial_smooth, 'acceleration');
            
            % Store results
            rate_maps_2d.(['ttg_accel_' group]) = rate_map;
            rate_maps_2d.(['ttg_accel_trial_counts_' group]) = trial_counts;
            rate_maps_2d.(['ttg_accel_n_trials_' group]) = n_trials;
            rate_maps_2d.(['ttg_accel_ttg_occupancy_' group]) = ttg_bin_occupancy;
            rate_maps_2d.(['ttg_bin_edges_' group]) = ttg_edges;
        end
    end
end


function [rate_map, trial_counts, ttg_bin_occupancy] = compute_ttg_kinematic_map_for_group(...
    trials, ttg_edges, n_ttg_bins, kinematic_edges, n_kinematic_bins, ...
    rate_per_trial_smooth, kinematic_type)
% COMPUTE_TTG_KINEMATIC_MAP_FOR_GROUP Compute 2D rate map for one group
%
% INPUTS:
%   trials              - Array of trial structs
%   ttg_edges           - Edges for TTG bins
%   n_ttg_bins          - Number of TTG bins
%   kinematic_edges     - Edges for velocity or acceleration bins
%   n_kinematic_bins    - Number of kinematic bins
%   rate_per_trial_smooth - Smoothed firing rates per trial (n_trials × n_ttg_bins)
%   kinematic_type      - 'velocity' or 'acceleration'
%
% OUTPUTS:
%   rate_map            - 2D rate map (n_ttg_bins × n_kinematic_bins)
%   trial_counts        - Number of trials contributing to each bin
%   ttg_bin_occupancy   - Number of trials visiting each TTG bin (1 × n_ttg_bins)

    n_trials = length(trials);
    rate_maps_per_trial = nan(n_trials, n_ttg_bins, n_kinematic_bins);
    bin_sampled_per_trial = zeros(n_trials, n_ttg_bins, n_kinematic_bins);
    ttg_bin_occupancy_per_trial = zeros(n_trials, n_ttg_bins);  % Track TTG bin occupancy separately
    
    for t = 1:n_trials
        trial = trials(t).trial;
        motion_mask = trial.motion_mask();
        
        % Get data during motion
        motion_times = trial.probe_t(motion_mask);
        if strcmp(kinematic_type, 'velocity')
            kinematic_data = trial.velocity(motion_mask);
        else
            kinematic_data = trial.acceleration(motion_mask);
        end
        
        n_motion_samples = length(motion_times);
        if n_motion_samples < 2
            continue;
        end
        
        % Create continuous time from number of samples (removes stationary gaps)
        % This matches the approach in compute_time_to_goal_tuning.m
        sampling_rate = 10000;  % Hz
        trial_duration = n_motion_samples / sampling_rate;  % Total motion duration in seconds
        continuous_time = linspace(0, trial_duration, n_motion_samples)';
        
        % Compute normalized TTG using continuous time
        ttg = trial_duration - continuous_time;
        ttg_norm = (ttg / trial_duration) * 100;
        
        % Initialize map for this trial
        rate_map_trial = nan(n_ttg_bins, n_kinematic_bins);
        
        % For each TTG bin
        for b = 1:n_ttg_bins
            in_ttg_bin = ttg_norm >= ttg_edges(b) & ttg_norm < ttg_edges(b+1);
            
            if any(in_ttg_bin)
                % Mark this TTG bin as occupied in this trial
                ttg_bin_occupancy_per_trial(t, b) = 1;
                
                kinematic_in_bin = kinematic_data(in_ttg_bin);
                
                % Assign to kinematic bins
                for k = 1:n_kinematic_bins
                    in_kinematic_bin = kinematic_in_bin >= kinematic_edges(k) & ...
                                       kinematic_in_bin < kinematic_edges(k+1);
                    
                    if any(in_kinematic_bin)
                        bin_sampled_per_trial(t, b, k) = 1;
                        rate_map_trial(b, k) = rate_per_trial_smooth(t, b);
                    end
                end
            end
        end
        rate_maps_per_trial(t, :, :) = rate_map_trial;
    end
    
    % Compute median and counts
    rate_map = squeeze(nanmedian(rate_maps_per_trial, 1));
    trial_counts = squeeze(sum(bin_sampled_per_trial, 1));
    
    % Compute TTG bin occupancy (how many trials visited each TTG bin)
    ttg_bin_occupancy = sum(ttg_bin_occupancy_per_trial, 1);  % Vector of length n_ttg_bins
end
