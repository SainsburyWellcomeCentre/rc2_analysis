function info = skaggs_info(rateMap, occupancy)
% SKAGGS_INFO Compute spatial information content (Skaggs et al. 1993)
%
%   info = skaggs_info(rateMap, occupancy)
%
%   Computes the spatial information content of a firing rate map, which
%   quantifies how much information the spike train conveys about the
%   animal's position. Higher values indicate more spatially specific firing.
%
%   Inputs:
%       rateMap   - Firing rate per spatial bin (Hz), 1D array
%       occupancy - Time spent in each spatial bin (s), 1D array
%
%   Output:
%       info - Spatial information in bits/spike
%
%   Formula: I = Σ p_i × (λ_i / λ) × log2(λ_i / λ)
%   where:
%       p_i  = probability of occupying bin i (normalized occupancy)
%       λ_i  = firing rate in bin i
%       λ    = overall mean firing rate (total spikes / total time)
%
%   Reference:
%       Skaggs, W. E., McNaughton, B. L., Gothard, K. M., & Markus, E. J. (1993).
%       An information-theoretic approach to deciphering the hippocampal code.
%       In Advances in neural information processing systems (pp. 1030-1037).
%
%   Example:
%       rate = [0.5, 2.3, 5.1, 2.8, 0.3];  % Hz
%       occ = [0.8, 1.2, 1.5, 1.0, 0.5];   % seconds
%       info = skaggs_info(rate, occ);     % bits/spike

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
    
    % Skaggs information
    lam_i = rateMap;
    info = nansum( p_i(mask) .* (lam_i(mask) ./ lam) .* log2(lam_i(mask) ./ lam) );
end
