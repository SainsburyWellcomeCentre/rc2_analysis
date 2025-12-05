% 3D visualization of spatial firing rate profiles
% Creates 3D plots showing:
%   - x-axis: position (cm) or percentage of track
%   - y-axis: cluster ID (categorical)
%   - z-axis: firing rate (Hz)
%
% This script is called by spatial_firing_rate_profile.m and uses variables
% from that script's workspace.
%
% Required variables from parent workspace:
%   - spatial_tuning_stats: structure with spatial tuning statistics
%   - cluster_ids: array of cluster IDs
%   - all_Q2_rate_smooth_by_group: median firing rates by group
%   - all_bin_centers_by_group: bin centers by group
%   - pastel_colors: colors for each cluster
%   - probe_ids{pid}: current probe ID
%   - ctl: RC2Analysis controller object
%   - n_clusters: number of clusters

if ~exist('spatial_tuning_stats', 'var')
    fprintf('  Skipping 3D visualization: spatial_tuning_stats not available\n');
    return;
end

fprintf('  Creating 3D spatial tuning visualization...\n');

% Identify clusters that are spatially tuned in AT LEAST ONE condition
spatially_tuned_any = [];

for c = 1:n_clusters
    cluster_field = sprintf('cluster_%d', cluster_ids(c));
    if isfield(spatial_tuning_stats, cluster_field)
        stats = spatial_tuning_stats.(cluster_field);
        % Must be tuned in at least one group
        if (isfield(stats, 'long') && stats.long.isSpatiallyTuned) || ...
           (isfield(stats, 'short') && stats.short.isSpatiallyTuned)
            spatially_tuned_any(end+1) = c;
        end
    end
end

n_tuned = length(spatially_tuned_any);
fprintf('    Found %d clusters spatially tuned in at least one condition\n', n_tuned);

if n_tuned == 0
    fprintf('    No spatially tuned clusters found. Skipping 3D plot.\n');
    return;
end

% Sort clusters by peak position in long trials (or short if not tuned in long)
peak_pos = zeros(n_tuned, 1);
for i = 1:n_tuned
    c = spatially_tuned_any(i);
    cluster_field = sprintf('cluster_%d', cluster_ids(c));
    stats = spatial_tuning_stats.(cluster_field);
    if isfield(stats, 'long') && stats.long.isSpatiallyTuned
        peak_pos(i) = stats.long.peakPosition;
    elseif isfield(stats, 'short') && stats.short.isSpatiallyTuned
        peak_pos(i) = stats.short.peakPosition;
    end
end
[~, sort_idx] = sort(peak_pos);
spatially_tuned_sorted = spatially_tuned_any(sort_idx);

% Generate pastel colors for all clusters (ensure unique colors)
pastel_colors = generate_pastel_colors(n_tuned);

% Define Gaussian fit function
gauss_fit = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2 * p(3)^2));
options = optimset('Display', 'off');

% Create figure with 2x2 subplots (portrait orientation)
fig_3d = ctl.figs.a4figure('portrait');

% Reduce resolution for faster rendering
n_gauss_points = 50;  % Reduced from 200 to speed up rendering

% --- Subplot 1: Long trials (unnormalized) ---
subplot(2, 2, 1, 'Parent', fig_3d);
hold on;
for i = 1:n_tuned
    c = spatially_tuned_sorted(i);
    cluster_id = cluster_ids(c);
    cluster_field = sprintf('cluster_%d', cluster_id);
    color = pastel_colors(i, :);
    
    bin_centers_long = all_bin_centers_by_group.long;
    median_rate_long = all_Q2_rate_smooth_by_group.long(c, :);
    
    % Create y-values (cluster ID repeated for each position)
    y_vals = ones(size(bin_centers_long)) * i;  % Use sorted index for y-axis
    
    % Plot median data points
    plot3(bin_centers_long, y_vals, median_rate_long, 'o', 'Color', color, ...
          'MarkerSize', 4, 'MarkerFaceColor', color);
    
    % Fit and plot Gaussian with shaded area if tuned in long
    if isfield(spatial_tuning_stats.(cluster_field), 'long') && ...
       spatial_tuning_stats.(cluster_field).long.isSpatiallyTuned
        peakRate = max(median_rate_long);
        peakIdx = find(median_rate_long == peakRate, 1);
        peakPos = bin_centers_long(peakIdx);
        sigma_guess = 10;
        p0 = [peakRate, peakPos, sigma_guess];
        try
            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_long, median_rate_long, [], [], options);
            x_fine = linspace(min(bin_centers_long), max(bin_centers_long), n_gauss_points);
            y_fine = gauss_fit(p_fit, x_fine);
            
            % Create shaded surface using surf() instead of fill3 for better rendering
            y_offset = 0.3;  % Width of the ribbon
            [X_mesh, Y_mesh] = meshgrid(x_fine, [i - y_offset/2, i + y_offset/2]);
            Z_mesh = repmat(y_fine, 2, 1);
            surf(X_mesh, Y_mesh, Z_mesh, 'FaceColor', color, 'EdgeColor', 'none', ...
                 'FaceAlpha', 0.6, 'FaceLighting', 'gouraud');
        catch
            % Fit failed, skip Gaussian
        end
    end
end
xlabel('Position (cm)');
ylabel('Cluster (sorted)');
zlabel('Firing Rate (Hz)');
title('Long Trials');
grid on;
xlim([0, 120]);
ylim([0.5, n_tuned + 0.5]);
view([45, 35]);  % More top-down view (azimuth=45°, elevation=35°)
% Set y-ticks to show cluster IDs
yticks(1:n_tuned);
yticklabels(arrayfun(@(i) sprintf('%d', cluster_ids(spatially_tuned_sorted(i))), 1:n_tuned, 'UniformOutput', false));
hold off;

% --- Subplot 2: Long trials (normalized) ---
subplot(2, 2, 2, 'Parent', fig_3d);
hold on;
for i = 1:n_tuned
    c = spatially_tuned_sorted(i);
    cluster_id = cluster_ids(c);
    cluster_field = sprintf('cluster_%d', cluster_id);
    color = pastel_colors(i, :);
    
    bin_centers_long = all_bin_centers_by_group.long;
    median_rate_long = all_Q2_rate_smooth_by_group.long(c, :);
    
    % Normalize to [0, 1]
    min_rate = min(median_rate_long);
    max_rate = max(median_rate_long);
    if max_rate > min_rate
        median_rate_norm = (median_rate_long - min_rate) / (max_rate - min_rate);
    else
        median_rate_norm = zeros(size(median_rate_long));
    end
    
    % Create y-values (cluster ID repeated for each position)
    y_vals = ones(size(bin_centers_long)) * i;
    
    % Plot normalized data points
    plot3(bin_centers_long, y_vals, median_rate_norm, 'o', 'Color', color, ...
          'MarkerSize', 4, 'MarkerFaceColor', color);
    
    % Fit and plot normalized Gaussian with shaded area if tuned in long
    if isfield(spatial_tuning_stats.(cluster_field), 'long') && ...
       spatial_tuning_stats.(cluster_field).long.isSpatiallyTuned
        peakIdx = find(median_rate_norm == max(median_rate_norm), 1);
        peakPos = bin_centers_long(peakIdx);
        sigma_guess = 10;
        p0 = [1.0, peakPos, sigma_guess];
        try
            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_long, median_rate_norm, [], [], options);
            x_fine = linspace(min(bin_centers_long), max(bin_centers_long), n_gauss_points);
            y_fine = gauss_fit(p_fit, x_fine);
            
            % Create shaded surface using surf() instead of fill3 for better rendering
            y_offset = 0.3;
            [X_mesh, Y_mesh] = meshgrid(x_fine, [i - y_offset/2, i + y_offset/2]);
            Z_mesh = repmat(y_fine, 2, 1);
            surf(X_mesh, Y_mesh, Z_mesh, 'FaceColor', color, 'EdgeColor', 'none', ...
                 'FaceAlpha', 0.6, 'FaceLighting', 'gouraud');
        catch
            % Fit failed, skip Gaussian
        end
    end
end
xlabel('Position (cm)');
ylabel('Cluster (sorted)');
zlabel('Normalized Firing Rate');
title('Long Trials (normalized)');
grid on;
xlim([0, 120]);
ylim([0.5, n_tuned + 0.5]);
zlim([0, 1]);
view([45, 35]);
yticks(1:n_tuned);
yticklabels(arrayfun(@(i) sprintf('%d', cluster_ids(spatially_tuned_sorted(i))), 1:n_tuned, 'UniformOutput', false));
hold off;

% --- Subplot 3: Short trials (unnormalized, 60-120 cm) ---
subplot(2, 2, 3, 'Parent', fig_3d);
hold on;
for i = 1:n_tuned
    c = spatially_tuned_sorted(i);
    cluster_id = cluster_ids(c);
    cluster_field = sprintf('cluster_%d', cluster_id);
    color = pastel_colors(i, :);
    
    bin_centers_short = all_bin_centers_by_group.short;
    median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
    
    % Create y-values (cluster ID repeated for each position)
    y_vals = ones(size(bin_centers_short)) * i;
    
    % Plot data points (using actual cm positions 60-120)
    plot3(bin_centers_short, y_vals, median_rate_short, 'o', 'Color', color, ...
          'MarkerSize', 4, 'MarkerFaceColor', color);
    
    % Fit and plot Gaussian with shaded area if tuned in short
    if isfield(spatial_tuning_stats.(cluster_field), 'short') && ...
       spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
        bin_centers_cm = bin_centers_short;
        peakRate = max(median_rate_short);
        peakIdx = find(median_rate_short == peakRate, 1);
        peakPos_cm = bin_centers_cm(peakIdx);
        sigma_guess = 5;
        p0 = [peakRate, peakPos_cm, sigma_guess];
        try
            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_cm, median_rate_short, [], [], options);
            x_fine = linspace(min(bin_centers_cm), max(bin_centers_cm), n_gauss_points);
            y_fine = gauss_fit(p_fit, x_fine);
            
            % Create shaded surface using surf() instead of fill3 for better rendering
            y_offset = 0.3;
            [X_mesh, Y_mesh] = meshgrid(x_fine, [i - y_offset/2, i + y_offset/2]);
            Z_mesh = repmat(y_fine, 2, 1);
            surf(X_mesh, Y_mesh, Z_mesh, 'FaceColor', color, 'EdgeColor', 'none', ...
                 'FaceAlpha', 0.6, 'FaceLighting', 'gouraud');
        catch
            % Fit failed, skip Gaussian
        end
    end
end
xlabel('Position (cm)');
ylabel('Cluster (sorted)');
zlabel('Firing Rate (Hz)');
title('Short Trials');
grid on;
xlim([60, 120]);
ylim([0.5, n_tuned + 0.5]);
view([45, 35]);
yticks(1:n_tuned);
yticklabels(arrayfun(@(i) sprintf('%d', cluster_ids(spatially_tuned_sorted(i))), 1:n_tuned, 'UniformOutput', false));
hold off;

% --- Subplot 4: Short trials (normalized, 60-120 cm) ---
subplot(2, 2, 4, 'Parent', fig_3d);
hold on;
for i = 1:n_tuned
    c = spatially_tuned_sorted(i);
    cluster_id = cluster_ids(c);
    cluster_field = sprintf('cluster_%d', cluster_id);
    color = pastel_colors(i, :);
    
    bin_centers_short = all_bin_centers_by_group.short;
    median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
    
    % Normalize to [0, 1]
    min_rate = min(median_rate_short);
    max_rate = max(median_rate_short);
    if max_rate > min_rate
        median_rate_norm = (median_rate_short - min_rate) / (max_rate - min_rate);
    else
        median_rate_norm = zeros(size(median_rate_short));
    end
    
    % Create y-values (cluster ID repeated for each position)
    y_vals = ones(size(bin_centers_short)) * i;
    
    % Plot normalized data points (using actual cm positions 60-120)
    plot3(bin_centers_short, y_vals, median_rate_norm, 'o', 'Color', color, ...
          'MarkerSize', 4, 'MarkerFaceColor', color);
    
    % Fit and plot normalized Gaussian with shaded area if tuned in short
    if isfield(spatial_tuning_stats.(cluster_field), 'short') && ...
       spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
        bin_centers_cm = bin_centers_short;
        peakIdx = find(median_rate_norm == max(median_rate_norm), 1);
        peakPos_cm = bin_centers_cm(peakIdx);
        sigma_guess = 5;
        p0 = [1.0, peakPos_cm, sigma_guess];
        try
            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_cm, median_rate_norm, [], [], options);
            x_fine = linspace(min(bin_centers_cm), max(bin_centers_cm), n_gauss_points);
            y_fine = gauss_fit(p_fit, x_fine);
            
            % Create shaded surface using surf() instead of fill3 for better rendering
            y_offset = 0.3;
            [X_mesh, Y_mesh] = meshgrid(x_fine, [i - y_offset/2, i + y_offset/2]);
            Z_mesh = repmat(y_fine, 2, 1);
            surf(X_mesh, Y_mesh, Z_mesh, 'FaceColor', color, 'EdgeColor', 'none', ...
                 'FaceAlpha', 0.6, 'FaceLighting', 'gouraud');
        catch
            % Fit failed, skip Gaussian
        end
    end
end
xlabel('Position (cm)');
ylabel('Cluster (sorted)');
zlabel('Normalized Firing Rate');
title('Short Trials (normalized)');
grid on;
xlim([60, 120]);
ylim([0.5, n_tuned + 0.5]);
zlim([0, 1]);
view([45, 35]);
yticks(1:n_tuned);
yticklabels(arrayfun(@(i) sprintf('%d', cluster_ids(spatially_tuned_sorted(i))), 1:n_tuned, 'UniformOutput', false));
hold off;

% Add figure title with adjusted position
fig_title = FigureTitle(fig_3d, sprintf('3D Spatial Tuning - %s', probe_ids{pid}));
fig_title.y_position = 0.99;

% Adjust subplot positions to create more space at top
all_axes = findall(fig_3d, 'Type', 'axes');
for ax_idx = 1:length(all_axes)
    ax = all_axes(ax_idx);
    if ~strcmp(get(ax, 'Visible'), 'off')  % Skip the hidden axis used by FigureTitle
        pos = get(ax, 'Position');
        pos(4) = pos(4) * 0.90;  % Reduce height by 10% to create space at top
        set(ax, 'Position', pos);
    end
end

% Force rendering before saving to avoid hanging
drawnow;

% Save figure - this integrates it into the PDF joining system
ctl.figs.save_fig_to_join();

fprintf('  3D spatial tuning visualization completed\n');
