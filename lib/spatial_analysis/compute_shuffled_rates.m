function rate_per_trial_shuf = compute_shuffled_rates(trial_data_cache, shift_offsets, edges, bin_size_cm, gauss_sigma_cm, n_bins)
% COMPUTE_SHUFFLED_RATES [DEPRECATED] Compute shuffled firing rates using circular shift
%
%   *** DEPRECATED: This function is no longer used. ***
%   Use compute_shuffled_spatial_info.m instead, which implements a bin-shuffle
%   approach analogous to velocity tuning in ShuffleTuning.m
%
%   rate_per_trial_shuf = compute_shuffled_rates(trial_data_cache, shift_offsets, ...
%                                                edges, bin_size_cm, gauss_sigma_cm, n_bins)
%
%   This function previously implemented circular shuffle testing by shifting spike
%   times within each trial. This approach has been replaced by a simpler bin-shuffle
%   method that completely shuffles the rate matrix to break trial-bin associations.
%
%   See also: compute_shuffled_spatial_info, ShuffleTuning

    warning('compute_shuffled_rates is deprecated. Use compute_shuffled_spatial_info instead.');
    
    % Original implementation preserved for backward compatibility

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
