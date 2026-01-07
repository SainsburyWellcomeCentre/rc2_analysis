function ttg_stats = compute_ttg_tuning_statistics(ttg_data, group_names, stats_params)
% COMPUTE_TTG_TUNING_STATISTICS Test significance of TTG tuning using bin-shuffle
%
%   ttg_stats = compute_ttg_tuning_statistics(ttg_data, group_names, stats_params)
%
%   Tests whether time-to-goal (TTG) tuning is significant by shuffling the
%   rate matrix to break trial-bin associations and comparing observed TTG
%   information against a null distribution.
%
%   Inputs:
%       ttg_data      - Struct with TTG data per group (from compute_time_to_goal_tuning)
%       group_names   - Cell array {'long', 'short'}
%       stats_params  - Struct with fields:
%                       .minSpikes - minimum spikes required
%                       .nShuf     - number of shuffles
%                       .pThresh   - significance threshold
%
%   Output:
%       ttg_stats - Struct with fields for each group:
%                   .info        - observed TTG information (bits/spike)
%                   .infoNull    - shuffled info values (1 x nShuf)
%                   .pVal        - p-value from shuffle test
%                   .isTTGTuned  - boolean significance flag
%                   .nSpikes     - total spike count
%                   .meanRate    - mean firing rate (Hz)

ttg_stats = struct();

for g = 1:length(group_names)
    group = group_names{g};
    
    % Initialize empty stats
    ttg_stats.(group) = struct(...
        'info', nan, ...
        'infoNull', [], ...
        'pVal', nan, ...
        'isTTGTuned', false, ...
        'nSpikes', 0, ...
        'meanRate', nan);
    
    if ~isfield(ttg_data, group)
        continue;
    end
    
    % Get data
    rate_per_trial_smooth = ttg_data.(group).rate_per_trial_smooth;
    trial_durations = ttg_data.(group).trial_durations;
    bin_edges = ttg_data.(group).bin_edges_norm;
    
    [n_trials, n_bins] = size(rate_per_trial_smooth);
    
    if n_trials == 0 || n_bins == 0
        continue;
    end
    
    % Compute occupancy per trial (constant across bins for normalized TTG)
    bin_width_percent = bin_edges(2) - bin_edges(1);
    occupancy_per_trial = nan(n_trials, n_bins);
    
    for t = 1:n_trials
        if trial_durations(t) > 0
            occ_per_bin = trial_durations(t) * (bin_width_percent / 100);
            occupancy_per_trial(t, :) = occ_per_bin;
        end
    end
    
    % Compute spike counts from rates
    spike_counts_per_trial = rate_per_trial_smooth .* occupancy_per_trial;
    
    % Total spikes and mean rate
    total_spikes = nansum(spike_counts_per_trial(:));
    total_time = nansum(occupancy_per_trial(:));
    mean_rate = total_spikes / total_time;
    
    ttg_stats.(group).nSpikes = total_spikes;
    ttg_stats.(group).meanRate = mean_rate;
    
    % Check minimum spikes criterion
    if total_spikes < stats_params.minSpikes
        continue;
    end
    
    % Compute observed TTG information
    % Aggregate across trials: sum spike counts and occupancy per bin
    total_spike_counts = nansum(spike_counts_per_trial, 1);
    total_occupancy = nansum(occupancy_per_trial, 1);
    
    % Rate map: average rate per bin across trials
    rate_map = total_spike_counts ./ total_occupancy;
    rate_map(total_occupancy == 0) = 0;
    
    % Compute observed info
    info_obs = ttg_info(rate_map, total_occupancy);
    ttg_stats.(group).info = info_obs;
    
    % Shuffle test: randomly permute rate matrix to break trial-bin associations
    info_null = nan(stats_params.nShuf, 1);
    
    for shuf = 1:stats_params.nShuf
        % Shuffle: randomly permute all elements of the rate matrix
        rate_shuffled = rate_per_trial_smooth(:);
        rate_shuffled = rate_shuffled(randperm(length(rate_shuffled)));
        rate_shuffled = reshape(rate_shuffled, n_trials, n_bins);
        
        % Recompute spike counts with shuffled rates
        spike_counts_shuffled = rate_shuffled .* occupancy_per_trial;
        
        % Aggregate shuffled data
        total_spike_counts_shuf = nansum(spike_counts_shuffled, 1);
        rate_map_shuf = total_spike_counts_shuf ./ total_occupancy;
        rate_map_shuf(total_occupancy == 0) = 0;
        
        % Compute shuffled info
        info_null(shuf) = ttg_info(rate_map_shuf, total_occupancy);
    end
    
    ttg_stats.(group).infoNull = info_null;
    
    % Compute p-value (one-tailed: observed > null)
    % Using +1 correction to avoid p=0 (same formula as spatial tuning)
    p_val = (sum(info_null >= info_obs) + 1) / (stats_params.nShuf + 1);
    ttg_stats.(group).pVal = p_val;
    
    % Determine significance
    ttg_stats.(group).isTTGTuned = (p_val < stats_params.pThresh);
end

end
