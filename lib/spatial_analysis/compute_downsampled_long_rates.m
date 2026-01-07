function [rate_per_trial_smooth, Q1_smooth, Q2_smooth, Q3_smooth, bin_centers, edges] = compute_downsampled_long_rates(trial_groups, cluster, bin_size_cm, gauss_sigma_cm)
% COMPUTE_DOWNSAMPLED_LONG_RATES Compute spatial firing rates with downsampled bins
%
%   [rate_per_trial_smooth, Q1, Q2, Q3, bin_centers, edges] = ...
%       compute_downsampled_long_rates(trial_groups, cluster, bin_size_cm, gauss_sigma_cm)
%
%   For fair comparison with short trials in relative position analysis, this function
%   processes long trials with HALF the spatial resolution (every other bin removed)
%   BEFORE computing firing rates and smoothing.
%
%   Inputs:
%       trial_groups   - Struct with 'long' field containing trial data
%       cluster        - Cluster object with spike_times
%       bin_size_cm    - Bin size (cm), but will use 2x this for downsampled version
%       gauss_sigma_cm - Gaussian smoothing sigma (cm)
%
%   Outputs:
%       rate_per_trial_smooth - Smoothed firing rates (n_trials × n_bins)
%       Q1_smooth      - First quartile across trials
%       Q2_smooth      - Median across trials
%       Q3_smooth      - Third quartile across trials
%       bin_centers    - Bin centers for downsampled bins
%       edges          - Bin edges for downsampled bins

    % Create downsampled edges (double the bin size)
    downsampled_bin_size = bin_size_cm * 2;
    edges = 0:downsampled_bin_size:120;
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    n_bins = length(bin_centers);
    
    trials_struct = trial_groups.long;
    n_trials = length(trials_struct);
    
    % Compute firing rates per trial with downsampled bins
    spike_counts_per_trial = nan(n_trials, n_bins);
    occupancy_per_trial = nan(n_trials, n_bins);
    
    for t = 1:n_trials
        trial = trials_struct(t).trial;
        [~, occ, ~, counts] = compute_trial_firing_rate(trial, cluster, edges);
        spike_counts_per_trial(t, :) = counts;
        occupancy_per_trial(t, :) = occ;
    end
    
    % Smooth spike counts and occupancy separately (same approach as original)
    sigma_bins = gauss_sigma_cm / downsampled_bin_size;
    spike_counts_smooth = zeros(size(spike_counts_per_trial));
    occupancy_smooth = zeros(size(occupancy_per_trial));
    
    for t = 1:n_trials
        spike_counts_smooth(t, :) = imgaussfilt(spike_counts_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
        occupancy_smooth(t, :) = imgaussfilt(occupancy_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
    end
    
    % Compute smoothed rates
    rate_per_trial_smooth = spike_counts_smooth ./ occupancy_smooth;
    rate_per_trial_smooth(occupancy_smooth == 0) = nan;
    
    % Compute quartiles across trials
    Q1_smooth = quantile(rate_per_trial_smooth, 0.25, 1);
    Q2_smooth = quantile(rate_per_trial_smooth, 0.50, 1);
    Q3_smooth = quantile(rate_per_trial_smooth, 0.75, 1);
end
