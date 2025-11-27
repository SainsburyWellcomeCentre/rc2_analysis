function results = compute_quartiles_and_smooth(rate_per_trial, bin_size_cm, gauss_sigma_cm, verbose)
% COMPUTE_QUARTILES_AND_SMOOTH Compute quartiles across trials and apply smoothing
%
%   results = compute_quartiles_and_smooth(rate_per_trial, bin_size_cm, gauss_sigma_cm)
%   results = compute_quartiles_and_smooth(rate_per_trial, bin_size_cm, gauss_sigma_cm, verbose)
%
%   Computes Q1 (25th percentile), Q2 (median), Q3 (75th percentile), and mean
%   firing rates across trials, then applies Gaussian smoothing to each.
%
%   Inputs:
%       rate_per_trial - Matrix of firing rates (n_trials x n_bins)
%       bin_size_cm    - Size of each spatial bin (cm)
%       gauss_sigma_cm - Standard deviation of Gaussian kernel (cm)
%       verbose        - (optional) If true, print statistics (default: false)
%
%   Outputs:
%       results - Structure containing:
%           .Q1_unsmooth  - 25th percentile (unsmoothed)
%           .Q1_smooth    - 25th percentile (smoothed)
%           .Q2_unsmooth  - 50th percentile / median (unsmoothed)
%           .Q2_smooth    - 50th percentile / median (smoothed)
%           .Q3_unsmooth  - 75th percentile (unsmoothed)
%           .Q3_smooth    - 75th percentile (smoothed)
%           .mean_unsmooth - Mean across trials (unsmoothed)
%           .mean_smooth   - Mean across trials (smoothed)
%           .n_trials      - Number of trials used
%
%   Example:
%       results = compute_quartiles_and_smooth(rate_per_trial, 2, 8, true);

    if nargin < 4
        verbose = false;
    end
    
    % Number of trials
    n_trials = size(rate_per_trial, 1);
    
    % Compute quartiles across trials (dimension 1)
    Q1_rate = prctile(rate_per_trial, 25, 1);
    Q2_rate = prctile(rate_per_trial, 50, 1);  % Median
    Q3_rate = prctile(rate_per_trial, 75, 1);
    mean_rate = nanmean(rate_per_trial, 1);
    
    % Print unsmoothed statistics if verbose
    if verbose
        fprintf('      Unsmoothed statistics (n=%d trials):\n', n_trials);
        fprintf('        Q1 (25%%): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(Q1_rate), min(Q1_rate), max(Q1_rate));
        fprintf('        Q2 (50%%, median): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(Q2_rate), min(Q2_rate), max(Q2_rate));
        fprintf('        Q3 (75%%): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(Q3_rate), min(Q3_rate), max(Q3_rate));
        fprintf('        Mean: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(mean_rate), min(mean_rate), max(mean_rate));
    end
    
    % Replace NaN with 0 for smoothing
    Q1_for_smooth = Q1_rate; Q1_for_smooth(isnan(Q1_for_smooth)) = 0;
    Q2_for_smooth = Q2_rate; Q2_for_smooth(isnan(Q2_for_smooth)) = 0;
    Q3_for_smooth = Q3_rate; Q3_for_smooth(isnan(Q3_for_smooth)) = 0;
    mean_for_smooth = mean_rate; mean_for_smooth(isnan(mean_for_smooth)) = 0;
    
    % Apply smoothing to each
    Q1_smooth = smooth_spatial_rate(Q1_for_smooth, bin_size_cm, gauss_sigma_cm);
    Q2_smooth = smooth_spatial_rate(Q2_for_smooth, bin_size_cm, gauss_sigma_cm);
    Q3_smooth = smooth_spatial_rate(Q3_for_smooth, bin_size_cm, gauss_sigma_cm);
    mean_smooth = smooth_spatial_rate(mean_for_smooth, bin_size_cm, gauss_sigma_cm);
    
    % Print smoothed statistics if verbose
    if verbose
        fprintf('      After smoothing (sigma=%.1f cm):\n', gauss_sigma_cm);
        fprintf('        Q1 (25%%): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(Q1_smooth), min(Q1_smooth), max(Q1_smooth));
        fprintf('        Q2 (50%%, median): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(Q2_smooth), min(Q2_smooth), max(Q2_smooth));
        fprintf('        Q3 (75%%): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(Q3_smooth), min(Q3_smooth), max(Q3_smooth));
        fprintf('        Mean: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                nanmean(mean_smooth), min(mean_smooth), max(mean_smooth));
    end
    
    % Package results
    results = struct();
    results.Q1_unsmooth = Q1_rate;
    results.Q1_smooth = Q1_smooth;
    results.Q2_unsmooth = Q2_rate;
    results.Q2_smooth = Q2_smooth;
    results.Q3_unsmooth = Q3_rate;
    results.Q3_smooth = Q3_smooth;
    results.mean_unsmooth = mean_rate;
    results.mean_smooth = mean_smooth;
    results.n_trials = n_trials;
end
