function fig_grid = plot_spatial_tuning_grid(spatial_tuning_stats, all_bin_centers_by_group, all_Q2_rate_smooth_by_group, cluster_ids, gauss_sigma_cm, probe_id)
% PLOT_SPATIAL_TUNING_GRID Create 2x3 grid plot with Gaussian fits for spatially tuned clusters
%
%   fig_grid = plot_spatial_tuning_grid(spatial_tuning_stats, all_bin_centers_by_group, ...
%       all_Q2_rate_smooth_by_group, cluster_ids, gauss_sigma_cm, probe_id)
%
%   Creates a 2x3 grid visualization showing spatially tuned clusters with
%   Gaussian fits overlaid on the median firing rates. Includes:
%   - Row 1: Long trials (0-120 cm) - unnormalized and normalized
%   - Row 2: Short trials (60-120 cm) - unnormalized and normalized
%   - Row 3: Short trials (0-100% of track) - unnormalized and normalized
%
%   Inputs:
%       spatial_tuning_stats         - Struct with tuning stats per cluster
%       all_bin_centers_by_group     - Struct with 'long' and 'short' bin centers
%       all_Q2_rate_smooth_by_group  - Struct with 'long' and 'short' median rates
%       cluster_ids                  - Array of cluster IDs
%       gauss_sigma_cm               - Gaussian sigma for initial fit guess
%       probe_id                     - String identifier for the probe (for title)
%
%   Returns:
%       fig_grid - Figure handle (empty if no spatially tuned clusters found)
%
%   Example:
%       fig = plot_spatial_tuning_grid(stats, bin_centers, median_rates, ids, 8, 'probe_001');

    % Identify clusters that are spatially tuned in AT LEAST ONE condition
    n_clusters = length(cluster_ids);
    spatially_tuned_any = [];
    
    for c = 1:n_clusters
        cluster_field = sprintf('cluster_%d', cluster_ids(c));
        if isfield(spatial_tuning_stats, cluster_field)
            stats = spatial_tuning_stats.(cluster_field);
            if (isfield(stats, 'long') && stats.long.isSpatiallyTuned) || ...
               (isfield(stats, 'short') && stats.short.isSpatiallyTuned)
                spatially_tuned_any = [spatially_tuned_any; c]; %#ok<AGROW>
            end
        end
    end
    
    n_tuned = length(spatially_tuned_any);
    fprintf('    Found %d clusters spatially tuned in at least one condition\n', n_tuned);
    
    if n_tuned == 0
        fprintf('    No spatially tuned clusters found. Skipping plot.\n');
        fig_grid = [];
        return;
    end
    
    % Sort by peak position
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
    
    % Generate pastel colors
    pastel_colors = generate_pastel_colors(n_tuned);
    
    % Create figure
    fprintf('  Creating 2x3 grid visualization of spatially tuned clusters...\n');
    fig_grid = figure('Position', [100, 100, 1200, 900], 'PaperPositionMode', 'auto', 'Visible', 'off');
    
    % Define Gaussian fit function
    gauss_fit = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2 * p(3)^2));
    options = optimset('Display', 'off');
    
    %% Row 1, Column 1: Long trials (unnormalized)
    subplot('Position', [0.08, 0.72, 0.38, 0.20], 'Parent', fig_grid);
    hold on;
    for i = 1:n_tuned
        c = spatially_tuned_sorted(i);
        cluster_id = cluster_ids(c);
        cluster_field = sprintf('cluster_%d', cluster_id);
        color = pastel_colors(i, :);
        
        bin_centers_long = all_bin_centers_by_group.long;
        median_rate_long = all_Q2_rate_smooth_by_group.long(c, :);
        
        h = plot(bin_centers_long, median_rate_long, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('C%d', cluster_id));
        h.Color(4) = 0.6;
        
        if isfield(spatial_tuning_stats.(cluster_field), 'long') && spatial_tuning_stats.(cluster_field).long.isSpatiallyTuned
            peak_pos = spatial_tuning_stats.(cluster_field).long.peakPosition;
            peak_rate = spatial_tuning_stats.(cluster_field).long.peakRate;
            p0 = [peak_rate, peak_pos, gauss_sigma_cm];
            try
                p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_long, median_rate_long, [0, 0, 1], [inf, 120, 50], options);
                x_fit = linspace(0, 120, 300);
                y_fit = gauss_fit(p_fit, x_fit);
                fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            catch
            end
        end
    end
    xlabel('Position (cm)'); ylabel('Firing Rate (Hz)'); title('Long Trials (0-120 cm)');
    grid on; xlim([0, 120]); hold off;
    
    %% Row 1, Column 2: Long trials (normalized)
    subplot('Position', [0.54, 0.72, 0.38, 0.20], 'Parent', fig_grid);
    hold on;
    for i = 1:n_tuned
        c = spatially_tuned_sorted(i);
        cluster_id = cluster_ids(c);
        cluster_field = sprintf('cluster_%d', cluster_id);
        color = pastel_colors(i, :);
        
        bin_centers_long = all_bin_centers_by_group.long;
        median_rate_long = all_Q2_rate_smooth_by_group.long(c, :);
        min_rate = min(median_rate_long); max_rate = max(median_rate_long);
        if max_rate > min_rate
            median_rate_norm = (median_rate_long - min_rate) / (max_rate - min_rate);
        else
            median_rate_norm = zeros(size(median_rate_long));
        end
        
        h = plot(bin_centers_long, median_rate_norm, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3);
        h.Color(4) = 0.6;
        
        if isfield(spatial_tuning_stats.(cluster_field), 'long') && spatial_tuning_stats.(cluster_field).long.isSpatiallyTuned
            peak_pos = spatial_tuning_stats.(cluster_field).long.peakPosition;
            p0 = [1.0, peak_pos, gauss_sigma_cm];
            try
                p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_long, median_rate_norm, [0, 0, 1], [2, 120, 50], options);
                x_fit = linspace(0, 120, 300);
                y_fit = gauss_fit(p_fit, x_fit);
                fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            catch
            end
        end
    end
    xlabel('Position (cm)'); ylabel('Normalized Firing Rate'); title('Long Trials (normalized)');
    grid on; xlim([0, 120]); ylim([0, 1]); hold off;
    
    %% Row 2, Column 2: Short trials (unnormalized, 60-120 cm)
    subplot('Position', [0.23, 0.42, 0.23, 0.20], 'Parent', fig_grid);
    hold on;
    for i = 1:n_tuned
        c = spatially_tuned_sorted(i);
        cluster_id = cluster_ids(c);
        cluster_field = sprintf('cluster_%d', cluster_id);
        color = pastel_colors(i, :);
        
        bin_centers_short = all_bin_centers_by_group.short;
        median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
        
        h = plot(bin_centers_short, median_rate_short, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3);
        h.Color(4) = 0.6;
        
        if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
            peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
            peak_rate = spatial_tuning_stats.(cluster_field).short.peakRate;
            p0 = [peak_rate, peak_pos, gauss_sigma_cm];
            try
                p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_short, [0, 60, 1], [inf, 120, 50], options);
                x_fit = linspace(60, 120, 300);
                y_fit = gauss_fit(p_fit, x_fit);
                fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            catch
            end
        end
    end
    xlabel('Position (cm)'); ylabel('Firing Rate (Hz)'); title('Short Trials (60-120 cm)');
    grid on; xlim([60, 120]); hold off;
    
    %% Row 2, Column 4: Short trials (normalized, 60-120 cm)
    subplot('Position', [0.69, 0.42, 0.23, 0.20], 'Parent', fig_grid);
    hold on;
    for i = 1:n_tuned
        c = spatially_tuned_sorted(i);
        cluster_id = cluster_ids(c);
        cluster_field = sprintf('cluster_%d', cluster_id);
        color = pastel_colors(i, :);
        
        bin_centers_short = all_bin_centers_by_group.short;
        median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
        min_rate = min(median_rate_short); max_rate = max(median_rate_short);
        if max_rate > min_rate
            median_rate_norm = (median_rate_short - min_rate) / (max_rate - min_rate);
        else
            median_rate_norm = zeros(size(median_rate_short));
        end
        
        h = plot(bin_centers_short, median_rate_norm, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3);
        h.Color(4) = 0.6;
        
        if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
            peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
            p0 = [1.0, peak_pos, gauss_sigma_cm];
            try
                p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_norm, [0, 60, 1], [2, 120, 50], options);
                x_fit = linspace(60, 120, 300);
                y_fit = gauss_fit(p_fit, x_fit);
                fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            catch
            end
        end
    end
    xlabel('Position (cm)'); ylabel('Normalized Firing Rate'); title('Short Trials (normalized, 60-120 cm)');
    grid on; xlim([60, 120]); ylim([0, 1]); hold off;
    
    %% Row 3, Column 1: Short trials with percentage x-axis
    subplot('Position', [0.08, 0.12, 0.38, 0.20], 'Parent', fig_grid);
    hold on;
    for i = 1:n_tuned
        c = spatially_tuned_sorted(i);
        cluster_id = cluster_ids(c);
        cluster_field = sprintf('cluster_%d', cluster_id);
        color = pastel_colors(i, :);
        
        bin_centers_short = all_bin_centers_by_group.short;
        median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
        bin_centers_percent = (bin_centers_short - 60) / 60 * 100;
        
        h = plot(bin_centers_percent, median_rate_short, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3);
        h.Color(4) = 0.6;
        
        if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
            peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
            peak_rate = spatial_tuning_stats.(cluster_field).short.peakRate;
            p0 = [peak_rate, peak_pos, gauss_sigma_cm];
            try
                p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_short, [0, 60, 1], [inf, 120, 50], options);
                x_fit_cm = linspace(60, 120, 300);
                x_fit_percent = (x_fit_cm - 60) / 60 * 100;
                y_fit = gauss_fit(p_fit, x_fit_cm);
                fill([x_fit_percent, fliplr(x_fit_percent)], [y_fit, zeros(size(y_fit))], color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            catch
            end
        end
    end
    xlabel('Track Position (%)'); ylabel('Firing Rate (Hz)'); title('Short Trials (0-100% of track)');
    grid on; xlim([0, 100]); hold off;
    
    %% Row 3, Column 2: Short trials normalized with percentage x-axis
    subplot('Position', [0.54, 0.12, 0.38, 0.20], 'Parent', fig_grid);
    hold on;
    for i = 1:n_tuned
        c = spatially_tuned_sorted(i);
        cluster_id = cluster_ids(c);
        cluster_field = sprintf('cluster_%d', cluster_id);
        color = pastel_colors(i, :);
        
        bin_centers_short = all_bin_centers_by_group.short;
        median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
        bin_centers_percent = (bin_centers_short - 60) / 60 * 100;
        min_rate = min(median_rate_short); max_rate = max(median_rate_short);
        if max_rate > min_rate
            median_rate_norm = (median_rate_short - min_rate) / (max_rate - min_rate);
        else
            median_rate_norm = zeros(size(median_rate_short));
        end
        
        h = plot(bin_centers_percent, median_rate_norm, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3);
        h.Color(4) = 0.6;
        
        if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
            peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
            p0 = [1.0, peak_pos, gauss_sigma_cm];
            try
                p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_norm, [0, 60, 1], [2, 120, 50], options);
                x_fit_cm = linspace(60, 120, 300);
                x_fit_percent = (x_fit_cm - 60) / 60 * 100;
                y_fit = gauss_fit(p_fit, x_fit_cm);
                fill([x_fit_percent, fliplr(x_fit_percent)], [y_fit, zeros(size(y_fit))], color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            catch
            end
        end
    end
    xlabel('Track Position (%)'); ylabel('Normalized Firing Rate'); title('Short Trials (normalized, 0-100% of track)');
    grid on; xlim([0, 100]); ylim([0, 1]); hold off;
    
    %% Adjust subplot positions
    all_axes = findall(fig_grid, 'Type', 'axes');
    for ax_idx = 1:length(all_axes)
        ax = all_axes(ax_idx);
        pos = get(ax, 'Position');
        pos(4) = pos(4) * 0.88;
        pos(2) = pos(2) * 0.85;
        set(ax, 'Position', pos);
    end
    
    % Add figure title
    FigureTitle(fig_grid, sprintf('Spatial Tuning Grid - %s', probe_id));
end
