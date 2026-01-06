function [ttg_data, ttg_norm_bin_centers] = compute_time_to_goal_tuning(cluster, trial_groups, group_names, all_global_trial_indices_by_group, cluster_idx)
% COMPUTE_TIME_TO_GOAL_TUNING Compute time-to-goal tuning for a cluster
%
%   [ttg_data, ttg_norm_bin_centers] = compute_time_to_goal_tuning(cluster, ...
%       trial_groups, group_names, all_global_trial_indices_by_group, cluster_idx)
%
% Inputs:
%   cluster - Cluster object with spike_times property
%   trial_groups - Struct with 'long' and 'short' fields containing trial arrays
%   group_names - Cell array of group names {'long', 'short'}
%   all_global_trial_indices_by_group - Struct with global trial indices per group
%   cluster_idx - Index of this cluster in the all_* arrays
%
% Outputs:
%   ttg_data - Struct with fields for each group containing:
%              - spike_ttg: Cell array of absolute TTG values per trial
%              - spike_ttg_norm: Cell array of normalized TTG values (0-100%) per trial
%              - bin_edges_norm: Bin edges for normalized TTG
%              - bin_centers_norm: Bin centers for normalized TTG
%              - trial_durations: Vector of trial durations
%              - rate_per_trial_smooth: Smoothed firing rates per trial
%   ttg_norm_bin_centers - Vector of normalized TTG bin centers (0-100%)

ttg_data = struct();

% First pass: find max TTG across all trials
max_ttg_observed = 0;
for g = 1:2
    group = group_names{g};
    trials = trial_groups.(group);
    n_trials = length(trials);
    
    for t = 1:n_trials
        trial = trials(t).trial;
        trial_duration = trial.probe_t(end) - trial.probe_t(1);
        if trial_duration > max_ttg_observed
            max_ttg_observed = trial_duration;
        end
    end
end

% Use the full TTG range (no occupancy threshold)
ttg_max = ceil(max_ttg_observed);

% Process each group
for g = 1:2
    group = group_names{g};
    trials = trial_groups.(group);
    n_trials = length(trials);
    
    spike_ttg_all = cell(n_trials, 1);  % Time-to-goal for each spike (absolute)
    spike_ttg_norm_all = cell(n_trials, 1);  % Time-to-goal normalized (0-100%)
    trial_durations = zeros(n_trials, 1);  % Duration of each trial in seconds
    global_indices = all_global_trial_indices_by_group.(group){cluster_idx};
    
    for t = 1:n_trials
        trial = trials(t).trial;  % Extract the trial object from the struct
        
        % Get motion mask and motion samples
        motion_mask = trial.motion_mask();
        motion_times = trial.probe_t(motion_mask);
        n_motion_samples = length(motion_times);
        
        if n_motion_samples < 2
            continue;  % Skip trials with insufficient motion data
        end
        
        % Get spike times for this cluster and trial
        all_spike_times = cluster.spike_times;
        
        % Create continuous time from number of samples (removes stationary gaps)
        n_samples = length(motion_times);
        sampling_rate = 10000;  % Hz
        trial_duration = n_samples / sampling_rate;  % Total motion duration in seconds
        
        % Create continuous time vector from 0 to trial_duration
        continuous_time = linspace(0, trial_duration, n_samples)';
        
        trial_durations(t) = trial_duration;  % Store duration
        
        % Filter spike times within motion time window
        mask = all_spike_times >= motion_times(1) & all_spike_times <= motion_times(end);
        spike_times = all_spike_times(mask);
        
        if isempty(spike_times)
            spike_ttg_norm_all{t} = [];
            continue;
        end
        
        % For each spike, find which motion sample it corresponds to
        % and compute TTG in the continuous time space
        spike_ttg_norm_per_spike = zeros(length(spike_times), 1);
        
        for s = 1:length(spike_times)
            % Find the closest motion sample to this spike
            [~, closest_idx] = min(abs(motion_times - spike_times(s)));
            % Get the continuous time for this sample
            spike_continuous_time = continuous_time(closest_idx);
            % Time to goal in continuous time
            ttg_time = trial_duration - spike_continuous_time;
            % Normalize by motion duration to get percentage
            spike_ttg_norm_per_spike(s) = (ttg_time / trial_duration) * 100;
        end
        
        % Store normalized TTG values (0-100%)
        spike_ttg_norm_all{t} = spike_ttg_norm_per_spike;
    end
    
    % Store for this group - NORMALIZED binning (0-100%)
    n_norm_bins = 50;  % 50 bins = 2% per bin
    ttg_norm_bin_edges = linspace(0, 100, n_norm_bins + 1);
    ttg_norm_bin_centers = (ttg_norm_bin_edges(1:end-1) + ttg_norm_bin_edges(2:end)) / 2;
    ttg_data.(group).spike_ttg = spike_ttg_all;
    ttg_data.(group).spike_ttg_norm = spike_ttg_norm_all;
    ttg_data.(group).bin_edges_norm = ttg_norm_bin_edges;
    ttg_data.(group).bin_centers_norm = ttg_norm_bin_centers;
    ttg_data.(group).trial_durations = trial_durations;
    ttg_data.(group).global_trial_indices = global_indices;
    ttg_data.(group).ttg_max = ttg_max;
    
    % --- Compute TTG median rates for heatmap plotting ---
    % Use time-based binning with constant occupancy per bin
    spike_ttg_norm = ttg_data.(group).spike_ttg_norm;
    bin_edges = ttg_data.(group).bin_edges_norm;
    n_bins = length(bin_edges) - 1;
    n_trials = length(spike_ttg_norm);
    bin_width_percent = bin_edges(2) - bin_edges(1);
    
    spike_counts_per_trial = nan(n_trials, n_bins);
    occupancy_per_trial = nan(n_trials, n_bins);
    
    for t = 1:n_trials
        trial_duration = trial_durations(t);
        
        if trial_duration > 0
            % Constant occupancy per bin based on motion duration and bin width
            occ_per_bin = trial_duration * (bin_width_percent / 100);
            
            % Compute spike counts per bin
            if ~isempty(spike_ttg_norm{t})
                spike_counts_per_trial(t, :) = histcounts(spike_ttg_norm{t}, bin_edges);
            else
                spike_counts_per_trial(t, :) = 0;
            end
            
            % Set constant occupancy for all bins
            occupancy_per_trial(t, :) = occ_per_bin;
        end
    end
    
    % Smooth spike counts and occupancy separately
    sigma_percent = 8;
    bin_width_percent = bin_edges(2) - bin_edges(1);
    sigma_bins = sigma_percent / bin_width_percent;
    spike_counts_smooth = zeros(size(spike_counts_per_trial));
    occupancy_smooth = zeros(size(occupancy_per_trial));
    for t = 1:n_trials
        spike_counts_smooth(t, :) = imgaussfilt(spike_counts_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
        occupancy_smooth(t, :) = imgaussfilt(occupancy_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
    end
    rate_per_trial_smooth = spike_counts_smooth ./ occupancy_smooth;
    rate_per_trial_smooth(occupancy_smooth == 0) = nan;
    
    % Store smoothed rates for 2D map computation
    ttg_data.(group).rate_per_trial_smooth = rate_per_trial_smooth;
end

% Return normalized bin centers (same for both groups)
ttg_norm_bin_centers = ttg_norm_bin_centers;

end
