function infoNull = compute_shuffled_spatial_info(rate_per_trial_smooth, occ_smooth, nShuf)
% COMPUTE_SHUFFLED_SPATIAL_INFO Bin-shuffle test for spatial tuning significance
%
%   infoNull = compute_shuffled_spatial_info(rate_per_trial_smooth, occ_smooth, nShuf)
%
%   Performs a shuffle test analogous to velocity tuning in ShuffleTuning.m.
%   Completely shuffles all rate matrix values to break trial-bin associations,
%   then recomputes Skaggs information for each shuffle to create a null distribution.
%
%   This approach tests whether the observed spatial structure is significant by
%   comparing the observed Skaggs information to a null distribution where spatial
%   relationships are destroyed.
%
%   Inputs:
%       rate_per_trial_smooth - (n_trials × n_bins) smoothed firing rates per trial
%       occ_smooth            - (1 × n_bins) average occupancy across trials
%       nShuf                 - number of shuffles to perform
%
%   Output:
%       infoNull - (nShuf × 1) vector of Skaggs information values for shuffled data
%
%   The shuffle procedure:
%       1. Take the rate matrix (n_trials × n_bins)
%       2. Shuffle ALL values with replacement using randi (like ShuffleTuning.m)
%       3. Compute median rate across shuffled trials
%       4. Calculate Skaggs information on shuffled median rate
%       5. Repeat nShuf times to build null distribution
%
%   See also: ShuffleTuning, skaggs_info, spatial_firing_rate_profile

    infoNull = zeros(nShuf, 1);
    
    % Shuffle with replacement (following ShuffleTuning.m line 121)
    for s = 1:nShuf
        % Completely shuffle all rate values - breaks trial-bin associations
        I = randi(numel(rate_per_trial_smooth), size(rate_per_trial_smooth));
        rate_shuffled = rate_per_trial_smooth(I);
        
        % Compute median of shuffled rates
        median_rate_shuf = nanmedian(rate_shuffled, 1);
        
        % Compute Skaggs info on shuffled data
        infoNull(s) = skaggs_info(median_rate_shuf, occ_smooth);
    end
end
