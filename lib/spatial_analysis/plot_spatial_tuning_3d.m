function fig_3d = plot_spatial_tuning_3d(spatial_tuning_stats, all_bin_centers_by_group, all_Q2_rate_smooth_by_group, cluster_ids, probe_id, ctl)
% PLOT_SPATIAL_TUNING_3D Create 3D visualization of spatial tuning
%
%   fig_3d = plot_spatial_tuning_3d(spatial_tuning_stats, all_bin_centers_by_group, ...
%       all_Q2_rate_smooth_by_group, cluster_ids, probe_id, ctl)
%
%   Creates a 2x2 grid of 3D plots showing spatially tuned clusters with
%   Gaussian fits overlaid on the median firing rates:
%   - Top row: Long trials (unnormalized and normalized)
%   - Bottom row: Short trials (unnormalized and normalized)
%
%   Inputs:
%       spatial_tuning_stats         - Struct with tuning stats per cluster
%       all_bin_centers_by_group     - Struct with 'long' and 'short' bin centers
%       all_Q2_rate_smooth_by_group  - Struct with 'long' and 'short' median rates
%       cluster_ids                  - Array of cluster IDs
%       probe_id                     - String identifier for the probe (for title)
%       ctl                          - RC2Analysis controller object
%
%   Returns:
%       fig_3d - Figure handle (empty if no spatially tuned clusters found)

    % Identify clusters that are spatially tuned in at least one condition
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
        fprintf('    No spatially tuned clusters found. Skipping 3D plot.\n');
        fig_3d = [];
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
    fprintf('  Creating 3D spatial tuning visualization...\n');
    fig_3d = ctl.figs.a4figure('portrait');
    
    % Gaussian fit function and options
    gauss_fit = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2 * p(3)^2));
    options = optimset('Display', 'off');
    n_gauss_points = 50;
    
    % Helper function to plot 3D subplot
    function plot_3d_subplot(subplot_pos, bin_centers, rates, trial_type, normalize, x_lim, title_str)
        subplot(2, 2, subplot_pos, 'Parent', fig_3d);
        hold on;
        
        for i = 1:n_tuned
            c = spatially_tuned_sorted(i);
            cluster_id = cluster_ids(c);
            cluster_field = sprintf('cluster_%d', cluster_id);
            color = pastel_colors(i, :);
            
            rate_data = rates(c, :);
            
            % Normalize if requested
            if normalize
                min_rate = min(rate_data);
                max_rate = max(rate_data);
                if max_rate > min_rate
                    rate_data = (rate_data - min_rate) / (max_rate - min_rate);
                else
                    rate_data = zeros(size(rate_data));
                end
            end
            
            % Create y-values (sorted index for each cluster)
            y_vals = ones(size(bin_centers)) * i;
            
            % Plot data points
            plot3(bin_centers, y_vals, rate_data, 'o', 'Color', color, ...
                  'MarkerSize', 4, 'MarkerFaceColor', color);
            
            % Fit and plot Gaussian if tuned
            if isfield(spatial_tuning_stats.(cluster_field), trial_type) && ...
               spatial_tuning_stats.(cluster_field).(trial_type).isSpatiallyTuned
                
                [~, peakIdx] = max(rate_data);
                peakPos = bin_centers(peakIdx);
                if strcmp(trial_type, 'long')
                    sigma_guess = 10;
                else
                    sigma_guess = 5;
                end
                if normalize
                    p0 = [1.0, peakPos, sigma_guess];
                else
                    p0 = [max(rate_data), peakPos, sigma_guess];
                end
                
                try
                    p_fit = lsqcurvefit(gauss_fit, p0, bin_centers, rate_data, [], [], options);
                    x_fine = linspace(min(bin_centers), max(bin_centers), n_gauss_points);
                    y_fine = gauss_fit(p_fit, x_fine);
                    
                    % Create shaded surface
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
        if normalize
            zlabel('Normalized Firing Rate');
        else
            zlabel('Firing Rate (Hz)');
        end
        title(title_str);
        grid on;
        xlim(x_lim);
        ylim([0.5, n_tuned + 0.5]);
        if normalize
            zlim([0, 1]);
        end
        view([45, 35]);
        yticks(1:n_tuned);
        yticklabels(arrayfun(@(i) sprintf('%d', cluster_ids(spatially_tuned_sorted(i))), ...
                            1:n_tuned, 'UniformOutput', false));
        hold off;
    end
    
    % Create 4 subplots
    plot_3d_subplot(1, all_bin_centers_by_group.long, all_Q2_rate_smooth_by_group.long, ...
                    'long', false, [0, 120], 'Long Trials');
    plot_3d_subplot(2, all_bin_centers_by_group.long, all_Q2_rate_smooth_by_group.long, ...
                    'long', true, [0, 120], 'Long Trials (normalized)');
    plot_3d_subplot(3, all_bin_centers_by_group.short, all_Q2_rate_smooth_by_group.short, ...
                    'short', false, [60, 120], 'Short Trials');
    plot_3d_subplot(4, all_bin_centers_by_group.short, all_Q2_rate_smooth_by_group.short, ...
                    'short', true, [60, 120], 'Short Trials (normalized)');
    
    % Add figure title
    fig_title = FigureTitle(fig_3d, sprintf('3D Spatial Tuning - %s', probe_id));
    fig_title.y_position = 0.99;
    
    % Adjust subplot positions
    all_axes = findall(fig_3d, 'Type', 'axes');
    for ax_idx = 1:length(all_axes)
        ax = all_axes(ax_idx);
        if ~strcmp(get(ax, 'Visible'), 'off')
            pos = get(ax, 'Position');
            pos(4) = pos(4) * 0.90;
            set(ax, 'Position', pos);
        end
    end
    
    drawnow;
    fprintf('  3D visualization completed\n');
end
