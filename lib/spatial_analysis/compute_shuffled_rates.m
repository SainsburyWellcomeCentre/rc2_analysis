function rate_per_trial_shuf = compute_shuffled_rates(trial_data_cache, shift_offsets, edges, bin_size_cm, gauss_sigma_cm, n_bins)
% COMPUTE_SHUFFLED_RATES Compute shuffled firing rates for all trials
%
%   rate_per_trial_shuf = compute_shuffled_rates(trial_data_cache, shift_offsets, ...
%                                                edges, bin_size_cm, gauss_sigma_cm, n_bins)
%
%   Efficiently computes shuffled firing rates for circular shuffle test by using
%   pre-cached trial data to avoid redundant loading and processing.
%
%   Inputs:
%       trial_data_cache - Cell array of structs containing pre-computed trial data:
%                           .pos - position vector
%                           .tvec - time vector
%                           .st_in_trial - spike times within trial
%                           .duration - trial duration
%                           .occ_smooth - smoothed occupancy
%       shift_offsets    - Array of random shift fractions [0,1] for each trial
%       edges            - Spatial bin edges (cm)
%       bin_size_cm      - Size of each spatial bin (cm)
%       gauss_sigma_cm   - Standard deviation of Gaussian smoothing kernel (cm)
%       n_bins           - Number of spatial bins
%
%   Output:
%       rate_per_trial_shuf - Shuffled firing rates (n_trials × n_bins)
%
%   This function is optimized for use in parfor loops for parallel shuffle testing.

    n_trials = length(trial_data_cache);
    rate_per_trial_shuf = nan(n_trials, n_bins);
    
    for k = 1:n_trials
        trial_cache = trial_data_cache{k};
        
        % Skip empty trials
        if isempty(trial_cache)
            continue;
        end
        
        % Skip trials with no spikes
        if isempty(trial_cache.st_in_trial)
            rate_per_trial_shuf(k, :) = 0;
            continue;
        end
        
        % Circular shift: add random offset and wrap
        shift_offset = shift_offsets(k) * trial_cache.duration;
        st_shifted = trial_cache.st_in_trial + shift_offset;
        
        % Wrap times that exceed trial duration
        tvec_start = trial_cache.tvec(1);
        tvec_end = trial_cache.tvec(end);
        st_shifted(st_shifted > tvec_end) = st_shifted(st_shifted > tvec_end) - trial_cache.duration;
        
        % Interpolate positions for shifted spikes
        spike_pos_shuf = interp1(trial_cache.tvec, trial_cache.pos, st_shifted, 'linear', 'extrap');
        
        % Compute rate for this shuffled trial
        spike_count_shuf = histcounts(spike_pos_shuf, edges);
        count_smooth_shuf = smooth_spatial_rate(spike_count_shuf, bin_size_cm, gauss_sigma_cm);
        rate_per_trial_shuf(k, :) = count_smooth_shuf ./ trial_cache.occ_smooth;
        rate_per_trial_shuf(k, isnan(rate_per_trial_shuf(k, :)) | isinf(rate_per_trial_shuf(k, :))) = 0;
    end
end
