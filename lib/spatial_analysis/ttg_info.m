function info = ttg_info(rateMap, occupancy)
% TTG_INFO Compute time-to-goal information content (analogous to Skaggs et al. 1993)
%
%   info = ttg_info(rateMap, occupancy)
%
%   Computes the mutual information between neural firing and time-to-goal,
%   quantifying how much information the spike train conveys about the time
%   remaining until the goal. Higher values indicate more specific tuning to
%   time-to-goal.
%
%   Note: Unlike spatial information where occupancy varies across bins, TTG
%   information with normalized time bins has approximately constant occupancy
%   (p_i ≈ 1/n_bins). This makes the measure primarily reflect firing rate
%   modulation strength rather than occupancy correction.
%
%   Inputs:
%       rateMap   - Firing rate per TTG bin (Hz), 1D array
%       occupancy - Time spent in each TTG bin (s), 1D array
%
%   Output:
%       info - TTG information in bits/spike
%
%   Formula: I = Σ p_i × (λ_i / λ) × log2(λ_i / λ)
%   where:
%       p_i  = probability of occupying bin i (normalized occupancy)
%       λ_i  = firing rate in bin i
%       λ    = overall mean firing rate (total spikes / total time)
%
%   Reference (spatial version):
%       Skaggs, W. E., McNaughton, B. L., Gothard, K. M., & Markus, E. J. (1993).
%       An information-theoretic approach to deciphering the hippocampal code.
%       In Advances in neural information processing systems (pp. 1030-1037).
%
%   Example:
%       rate = [0.5, 2.3, 5.1, 2.8, 0.3];  % Hz
%       occ = [0.8, 1.2, 1.5, 1.0, 0.5];   % seconds
%       info = ttg_info(rate, occ);        % bits/spike

    % Normalize occupancy to probability
    p_i = occupancy / sum(occupancy);
    
    % Overall mean firing rate (total spikes / total time)
    % Since rateMap = spikes/occupancy, we can recover total spikes
    lam = sum(rateMap .* occupancy) / sum(occupancy);
    
    % Avoid log(0) and division by zero
    mask = (rateMap > 0) & (p_i > 0) & (lam > 0);
    
    if ~any(mask) || lam == 0
        info = 0;
        return;
    end
    
    % TTG information
    lam_i = rateMap;
    info = nansum( p_i(mask) .* (lam_i(mask) ./ lam) .* log2(lam_i(mask) ./ lam) );
end
