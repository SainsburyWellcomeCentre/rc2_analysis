function fig = plot_time_to_goal_heatmaps(all_ttg_rates, all_ttg_fields, ttg_norm_bin_centers, cluster_ids, probe_id, ctl)
% PLOT_TIME_TO_GOAL_HEATMAPS Create TTG heatmaps sorted by peak position
%
%   fig = plot_time_to_goal_heatmaps(all_ttg_rates, all_ttg_fields, ...
%         ttg_norm_bin_centers, cluster_ids, probe_id, ctl)
%
% Inputs:
%   all_ttg_rates - Struct with 'long' and 'short' fields containing
%                   n_clusters × n_bins matrices of median firing rates
%   all_ttg_fields - Cell array (n_clusters × 1) containing field detection results
%   ttg_norm_bin_centers - Vector of bin centers (0-100%)
%   cluster_ids - Vector of cluster IDs
%   probe_id - String identifier for the probe
%   ctl - Control structure with figure management
%
% Outputs:
%   fig - Figure handle
%
% The figure contains 1 row × 2 columns:
%   Col 1: Long trials TTG heatmap
%   Col 2: Short trials TTG heatmap
% Clusters are sorted by the peak position (max median rate) in the long condition

n_clusters = length(cluster_ids);

% Create figure
fig = ctl.figs.a4figure('landscape');

% Collect peak locations for sorting (using median rates, same as position heatmaps)
peak_locations = zeros(n_clusters, 1);

for c = 1:n_clusters
    % Find the bin with maximum median firing rate in long trials
    [~, max_idx] = max(all_ttg_rates.long(c, :));
    peak_locations(c) = ttg_norm_bin_centers(max_idx);
end

% Sort clusters by peak position
[~, sort_perm] = sort(peak_locations);
sorted_indices = sort_perm;

% Prepare data for plotting
ttg_bin_centers = ttg_norm_bin_centers; % 0-100%

% Plot Long (Col 1) and Short (Col 2)
for col = 1:2
    if col == 1
        group_name = 'long';
        title_str = 'Long Trials';
    else
        group_name = 'short';
        title_str = 'Short Trials';
    end
    
    subplot(1, 2, col, 'Parent', fig);
    
    if n_clusters > 0
        % Extract rates
        rates = all_ttg_rates.(group_name)(sorted_indices, :);
        
        % Normalize each row (0-1)
        min_vals = min(rates, [], 2);
        max_vals = max(rates, [], 2);
        range_vals = max_vals - min_vals;
        range_vals(range_vals == 0) = 1; % Avoid div by zero
        rates_norm = (rates - min_vals) ./ range_vals;
        
        imagesc(ttg_bin_centers, 1:n_clusters, rates_norm);
        colormap(gca, 'jet');
    else
        axis off;
        text(0.5, 0.5, 'No clusters', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center');
    end
    
    xlabel('Time to Goal (%)');
    if col == 1
        ylabel(sprintf('Clusters (n=%d)', n_clusters));
    else
        set(gca, 'YTickLabel', []);
    end
    
    title(title_str);
    set(gca, 'XDir', 'reverse'); % Time to goal: 100% -> 0%
    xlim([0, 100]);
end

% Adjust subplot positions to make them larger
all_axes = findall(fig, 'Type', 'axes');
for ax_idx = 1:length(all_axes)
    ax = all_axes(ax_idx);
    pos = get(ax, 'Position');
    % Increase width and height
    pos(3) = pos(3) * 1.15;  % Width
    pos(4) = pos(4) * 1.2;   % Height
    % Adjust positioning
    if pos(1) < 0.5  % Left subplot
        pos(1) = pos(1) - 0.02;
    else  % Right subplot
        pos(1) = pos(1) + 0.02;
    end
    pos(2) = pos(2) - 0.05;
    set(ax, 'Position', pos);
end

FigureTitle(fig, sprintf('Time-to-Goal Heatmaps - %s', probe_id));

end
