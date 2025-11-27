function [rate, occupancy, avg_velocity] = compute_trial_firing_rate(trial, cluster, edges)
% COMPUTE_TRIAL_FIRING_RATE Compute spatial firing rate for a single trial
%
%   [rate, occupancy, avg_velocity] = compute_trial_firing_rate(trial, cluster, edges)
%
%   Computes the firing rate per spatial bin for a single trial by dividing
%   spike counts by occupancy time in each bin.
%
%   Inputs:
%       trial   - Trial object with motion_mask, position, probe_t, velocity methods
%       cluster - Cluster object with spike_times property
%       edges   - Bin edges for spatial binning (cm), e.g., 0:2:120
%
%   Outputs:
%       rate         - Firing rate per spatial bin (Hz), NaN for unoccupied bins
%       occupancy    - Time spent in each bin (s)
%       avg_velocity - Average velocity in each bin (cm/s)
%
%   Example:
%       edges = 0:2:120;  % 2-cm bins from 0 to 120 cm
%       [rate, occ, vel] = compute_trial_firing_rate(trial, cluster, edges);

    n_bins = length(edges) - 1;
    
    % Get motion data from trial
    motion_mask = trial.motion_mask();
    pos = trial.position(motion_mask);
    tvec = trial.probe_t(motion_mask);
    vel = trial.velocity(motion_mask);
    
    % Handle empty trials
    if isempty(pos)
        rate = nan(1, n_bins);
        occupancy = nan(1, n_bins);
        avg_velocity = nan(1, n_bins);
        return;
    end
    
    % Get spike times within trial time window
    st = cluster.spike_times;
    mask = st >= tvec(1) & st <= tvec(end);
    st = st(mask);
    
    % Interpolate spike positions from spike times
    if ~isempty(st)
        spike_pos = interp1(tvec, pos, st, 'linear', 'extrap');
    else
        spike_pos = [];
    end
    
    % Compute spike count per bin
    spike_count = histcounts(spike_pos, edges);
    
    % Compute occupancy and velocity per bin
    dt_vec = [diff(tvec); 0];
    occupancy = zeros(1, n_bins);
    velocity_sum = zeros(1, n_bins);
    
    for i = 1:n_bins
        in_bin = pos >= edges(i) & pos < edges(i+1);
        occupancy(i) = sum(dt_vec(in_bin));
        velocity_sum(i) = sum(vel(in_bin) .* dt_vec(in_bin));
    end
    
    % Compute firing rate (spikes / occupancy time)
    rate = spike_count ./ occupancy;
    rate(isnan(rate) | isinf(rate)) = NaN;  % Mark unoccupied bins as NaN
    
    % Compute time-weighted average velocity
    avg_velocity = velocity_sum ./ occupancy;
    avg_velocity(isnan(avg_velocity) | isinf(avg_velocity)) = NaN;
end
