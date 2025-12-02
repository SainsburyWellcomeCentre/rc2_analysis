function fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, group_names, group_labels, group_colors, probe_id)
% PLOT_CLUSTER_SPATIAL_PROFILE Create a combined figure for a single cluster
%
%   fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, ...
%                                      group_names, group_labels, group_colors, probe_id)
%
%   Creates a 3-panel figure (1 row) showing:
%       - Panel 1: Combined position-based raster for both long and short trials
%       - Panel 2: Smoothed traces with median and IQR shading, plus pooled rate (dashed)
%       - Panel 3: X-normalized comparison with percentage x-axis, plus pooled rate (dashed)
%
%   Inputs:
%       cluster_id         - Numeric ID of the cluster
%       bin_centers_by_group - Struct with fields 'long' and 'short' containing bin centers
%       rate_data          - Struct with fields for each group containing:
%                              .Q1_smooth, .Q2_smooth, .Q3_smooth
%                              .rate_per_trial (raw rates, n_trials x n_bins)
%                              .rate_per_trial_smooth (smoothed rates, n_trials x n_bins)
%                              .spike_positions (cell array of spike positions per trial)
%                              .rate_pooled (pooled trial rate: smooth(sum(counts))/smooth(sum(occ)))
%       group_names        - Cell array {'long', 'short'}
%       group_labels       - Cell array with descriptive labels for legend
%       group_colors       - Struct with RGB colors for each group
%       probe_id           - String identifier for the probe (for title)
%
%   Outputs:
%       fig - Figure handle

    fig = figure('Position', [50, 50, 1500, 400]);
    
    % Compute max firing rate as 120% of the third quartile across all data
    all_Q3_values = [];
    for g = 1:length(group_names)
        group = group_names{g};
        if isfield(rate_data, group) && isfield(rate_data.(group), 'Q3_smooth')
            all_Q3_values = [all_Q3_values, rate_data.(group).Q3_smooth];
        end
    end
    
    % Handle edge cases for max firing rate calculation
    if isempty(all_Q3_values) || all(isnan(all_Q3_values))
        max_firing_rate = 10;  % Default fallback value
    else
        max_Q3 = max(all_Q3_values);
        if isnan(max_Q3) || max_Q3 <= 0
            max_firing_rate = 10;  % Default fallback value
        else
            max_firing_rate = max_Q3 * 1.2;  % Set max to 120% of third quartile
        end
    end
    
    %% Panel 1 (col 1): Combined position-based raster for both long and short trials
    subplot(1, 3, 1);
    hold on;
    
    % Collect all trials with their global indices and group info
    all_trials_data = [];
    
    % Collect long trials
    if isfield(rate_data, 'long') && isfield(rate_data.long, 'spike_positions') && ...
       isfield(rate_data.long, 'global_trial_indices')
        spike_positions = rate_data.long.spike_positions;
        global_indices = rate_data.long.global_trial_indices;
        n_long_trials = length(spike_positions);
        color = group_colors.long;
        
        for t = 1:n_long_trials
            all_trials_data = [all_trials_data; struct(...
                'global_idx', global_indices(t), ...
                'spike_positions', {spike_positions{t}}, ...
                'color', color, ...
                'group', 'long')]; %#ok<AGROW>
        end
    end
    
    % Collect short trials
    if isfield(rate_data, 'short') && isfield(rate_data.short, 'spike_positions') && ...
       isfield(rate_data.short, 'global_trial_indices')
        spike_positions = rate_data.short.spike_positions;
        global_indices = rate_data.short.global_trial_indices;
        n_short_trials = length(spike_positions);
        color = group_colors.short;
        
        for t = 1:n_short_trials
            all_trials_data = [all_trials_data; struct(...
                'global_idx', global_indices(t), ...
                'spike_positions', {spike_positions{t}}, ...
                'color', color, ...
                'group', 'short')]; %#ok<AGROW>
        end
    end
    
    % Sort by global trial index to show actual trial order
    if ~isempty(all_trials_data)
        [~, sort_idx] = sort([all_trials_data.global_idx]);
        all_trials_data = all_trials_data(sort_idx);
        
        % Plot trials in original order
        for t = 1:length(all_trials_data)
            trial_data = all_trials_data(t);
            if ~isempty(trial_data.spike_positions)
                n_spikes = length(trial_data.spike_positions);
                scatter(trial_data.spike_positions, t*ones(1, n_spikes), 8, ...
                       trial_data.color, 'filled', 'MarkerFaceAlpha', 0.6);
            end
        end
        
        set(gca, 'YDir', 'reverse', 'YLim', [0, length(all_trials_data)+1]);
        if length(all_trials_data) > 1
            set(gca, 'YTick', [1, length(all_trials_data)]);
        end
    end
    
    hold off;
    xlabel('Position (cm)');
    ylabel('Trial # (chronological)');
    title('Position Raster (Combined)');
    grid on;
    xlim([0, 120]);
    set(gca, 'XTick', 0:20:120);
    
    %% Panel 2 (col 2): Smoothed traces with median and IQR shading
    subplot(1, 3, 2);
    hold on;
    
    for g = 1:length(group_names)
        group = group_names{g};
        
        if ~isfield(rate_data, group)
            continue;
        end
        
        bin_centers = bin_centers_by_group.(group);
        rate_per_trial_smooth = rate_data.(group).rate_per_trial_smooth;
        Q1_smooth = rate_data.(group).Q1_smooth;
        Q2_smooth = rate_data.(group).Q2_smooth;
        Q3_smooth = rate_data.(group).Q3_smooth;
        
        color = group_colors.(group);
        
        % Plot Q1-Q3 shaded area (interquartile range)
        fill_x = [bin_centers, fliplr(bin_centers)];
        fill_y = [Q1_smooth, fliplr(Q3_smooth)];
        fill(fill_x, fill_y, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
        
        % Plot individual smoothed trials as thin lines
        if ~isempty(rate_per_trial_smooth)
            n_trials = size(rate_per_trial_smooth, 1);
            for t = 1:n_trials
                plot(bin_centers, rate_per_trial_smooth(t, :), '-', 'LineWidth', 0.5, 'Color', [color, 0.3], ...
                     'HandleVisibility', 'off');
            end
        end
        
        % Plot median (Q2) as thick solid line
        plot(bin_centers, Q2_smooth, '-', 'LineWidth', 2.5, 'Color', color, ...
             'HandleVisibility', 'off');
        
        % Plot pooled rate as thick dashed line if available
        if isfield(rate_data.(group), 'rate_pooled')
            rate_pooled = rate_data.(group).rate_pooled;
            plot(bin_centers, rate_pooled, '--', 'LineWidth', 2.5, 'Color', color, ...
                 'HandleVisibility', 'off');
        end
    end
    
    hold off;
    xlabel('Position (cm)');
    ylabel('Firing rate (Hz)');
    title('Smoothed Rates with Median & IQR');
    grid on;
    xlim([0, 120]);
    ylim([0, max_firing_rate]);
    set(gca, 'XTick', 0:20:120);
    
    %% Panel 3 (col 3): X-normalized comparison
    subplot(1, 3, 3);
    hold on;
    
    % Common normalized grid (0 to 1, displayed as 0% to 100%)
    x_norm = linspace(0, 1, 100);
    x_percent = x_norm * 100;  % For display as percentage
    
    for g = 1:length(group_names)
        group = group_names{g};
        
        if ~isfield(rate_data, group)
            continue;
        end
        
        pos = bin_centers_by_group.(group)(:);
        Q1_smooth = rate_data.(group).Q1_smooth(:);
        Q2_smooth = rate_data.(group).Q2_smooth(:);
        Q3_smooth = rate_data.(group).Q3_smooth(:);
        
        % Skip if all NaN
        if all(isnan(Q2_smooth))
            continue;
        end
        
        % Normalize position (start=0, end=1)
        pos_norm = (pos - min(pos)) / (max(pos) - min(pos));
        
        % Keep only finite points
        idx = isfinite(pos_norm) & isfinite(Q2_smooth);
        
        if sum(idx) < 2
            continue;
        end
        
        pos_n = pos_norm(idx);
        Q1_n = Q1_smooth(idx);
        Q2_n = Q2_smooth(idx);
        Q3_n = Q3_smooth(idx);
        
        % Interpolate to common grid
        Q1_i = interp1(pos_n, Q1_n, x_norm, 'linear', 'extrap');
        Q2_i = interp1(pos_n, Q2_n, x_norm, 'linear', 'extrap');
        Q3_i = interp1(pos_n, Q3_n, x_norm, 'linear', 'extrap');
        
        color = group_colors.(group);
        
        % Plot Q1-Q3 shaded area
        fill_x = [x_percent, fliplr(x_percent)];
        fill_y = [Q1_i, fliplr(Q3_i)];
        fill(fill_x, fill_y, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
        
        % Plot median as solid line
        plot(x_percent, Q2_i, '-', 'LineWidth', 2, 'Color', color, ...
             'HandleVisibility', 'off');
        
        % Plot pooled rate as dashed line if available
        if isfield(rate_data.(group), 'rate_pooled')
            rate_pooled = rate_data.(group).rate_pooled(:);
            idx_pooled = isfinite(pos_norm) & isfinite(rate_pooled);
            
            if sum(idx_pooled) >= 2
                pos_n_pooled = pos_norm(idx_pooled);
                rate_pooled_n = rate_pooled(idx_pooled);
                rate_pooled_i = interp1(pos_n_pooled, rate_pooled_n, x_norm, 'linear', 'extrap');
                plot(x_percent, rate_pooled_i, '--', 'LineWidth', 2, 'Color', color, ...
                     'HandleVisibility', 'off');
            end
        end
    end
    
    % Compute correlation between long and short medians
    if isfield(rate_data, 'long') && isfield(rate_data, 'short')
        pos_long = bin_centers_by_group.long(:);
        Q2_long = rate_data.long.Q2_smooth(:);
        pos_short = bin_centers_by_group.short(:);
        Q2_short = rate_data.short.Q2_smooth(:);
        
        pos_long_norm = (pos_long - min(pos_long)) / (max(pos_long) - min(pos_long));
        pos_short_norm = (pos_short - min(pos_short)) / (max(pos_short) - min(pos_short));
        
        idxL = isfinite(pos_long_norm) & isfinite(Q2_long);
        idxS = isfinite(pos_short_norm) & isfinite(Q2_short);
        
        if sum(idxL) >= 2 && sum(idxS) >= 2
            Q2_long_i = interp1(pos_long_norm(idxL), Q2_long(idxL), x_norm, 'linear', 'extrap');
            Q2_short_i = interp1(pos_short_norm(idxS), Q2_short(idxS), x_norm, 'linear', 'extrap');
            r = corr(Q2_long_i(:), Q2_short_i(:), 'rows', 'complete');
            title(sprintf('Normalized Comparison (r = %.2f)', r));
        else
            title('Normalized Comparison');
        end
    else
        title('Normalized Comparison');
    end
    
    hold off;
    xlabel('Track position (%)');
    ylabel('Firing rate (Hz)');
    grid on;
    xlim([0, 100]);
    ylim([0, max_firing_rate]);
    set(gca, 'XTick', 0:10:100);
    
    % Overall figure title
    sgtitle(sprintf('Cluster %d - %s', cluster_id, probe_id), 'FontWeight', 'bold');
end
