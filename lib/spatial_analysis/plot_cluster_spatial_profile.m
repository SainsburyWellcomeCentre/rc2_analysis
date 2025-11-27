function fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, group_names, group_labels, group_colors, probe_id)
% PLOT_CLUSTER_SPATIAL_PROFILE Create a combined figure for a single cluster
%
%   fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, ...
%                                      group_names, group_labels, group_colors, probe_id)
%
%   Creates a 2-panel figure showing:
%     - Left: Spatial firing rate profile with IQR shading (absolute position)
%     - Right: X-normalized comparison between long and short tracks
%
%   Inputs:
%       cluster_id         - Numeric ID of the cluster
%       bin_centers_by_group - Struct with fields 'long' and 'short' containing bin centers
%       rate_data          - Struct with fields for each group containing:
%                              .mean_smooth, .Q1_smooth, .Q2_smooth, .Q3_smooth
%       group_names        - Cell array {'long', 'short'}
%       group_labels       - Cell array with descriptive labels for legend
%       group_colors       - Struct with RGB colors for each group
%       probe_id           - String identifier for the probe (for title)
%
%   Outputs:
%       fig - Figure handle
%
%   Example:
%       rate_data.long.mean_smooth = mean_rate_long;
%       rate_data.long.Q1_smooth = Q1_long;
%       % ... etc
%       fig = plot_cluster_spatial_profile(42, bin_centers, rate_data, ...
%                                          {'long','short'}, {'Long','Short'}, colors, 'probe1');

    fig = figure('Position', [100, 100, 1200, 500]);
    
    %% Panel 1: Spatial firing rate with IQR shading
    subplot(1, 2, 1);
    hold on;
    
    for g = 1:length(group_names)
        group = group_names{g};
        
        if ~isfield(rate_data, group)
            continue;
        end
        
        bin_centers = bin_centers_by_group.(group);
        rate_smooth = rate_data.(group).mean_smooth;
        Q1_smooth = rate_data.(group).Q1_smooth;
        Q2_smooth = rate_data.(group).Q2_smooth;
        Q3_smooth = rate_data.(group).Q3_smooth;
        
        % Skip if all NaN
        if all(isnan(rate_smooth))
            continue;
        end
        
        % Plot Q1-Q3 shaded area (interquartile range)
        fill_x = [bin_centers, fliplr(bin_centers)];
        fill_y = [Q1_smooth, fliplr(Q3_smooth)];
        fill(fill_x, fill_y, group_colors.(group), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
        
        % Plot median (Q2) as main line
        plot(bin_centers, Q2_smooth, '-', 'LineWidth', 2, 'Color', group_colors.(group), ...
             'DisplayName', sprintf('%s median', group_labels{g}));
        
        % Plot mean as dashed line for comparison
        plot(bin_centers, rate_smooth, '--', 'LineWidth', 1.5, 'Color', group_colors.(group), ...
             'DisplayName', sprintf('%s mean', group_labels{g}));
    end
    
    hold off;
    xlabel('Position (cm)');
    ylabel('Firing rate (Hz)');
    legend('show', 'Location', 'best');
    title('Spatial Firing Rate Profile');
    grid on;
    
    %% Panel 2: X-normalized comparison
    subplot(1, 2, 2);
    
    % Get data for both groups
    pos_long = bin_centers_by_group.long(:);
    rate_long = rate_data.long.mean_smooth(:);
    pos_short = bin_centers_by_group.short(:);
    rate_short = rate_data.short.mean_smooth(:);
    
    % Check if we have valid data
    has_valid_data = ~all(isnan(rate_long)) && ~all(isnan(rate_short));
    
    if has_valid_data
        % Normalize position (start=0, end=1)
        pos_long_norm = (pos_long - min(pos_long)) / (max(pos_long) - min(pos_long));
        pos_short_norm = (pos_short - min(pos_short)) / (max(pos_short) - min(pos_short));
        
        % Keep only finite points
        idxL = isfinite(pos_long_norm) & isfinite(rate_long);
        idxS = isfinite(pos_short_norm) & isfinite(rate_short);
        
        pos_long_n = pos_long_norm(idxL);
        rate_long_n = rate_long(idxL);
        pos_short_n = pos_short_norm(idxS);
        rate_short_n = rate_short(idxS);
        
        if numel(pos_long_n) >= 2 && numel(pos_short_n) >= 2
            % Common normalized grid
            x_norm = linspace(0, 1, 100);
            
            % Interpolate both curves
            rate_long_i = interp1(pos_long_n, rate_long_n, x_norm, 'linear', 'extrap');
            rate_short_i = interp1(pos_short_n, rate_short_n, x_norm, 'linear', 'extrap');
            
            % Compute correlation
            r = corr(rate_long_i(:), rate_short_i(:), 'rows', 'complete');
            
            % Plot
            hold on;
            plot(x_norm, rate_long_i, '-', 'LineWidth', 2, 'Color', group_colors.long, ...
                 'DisplayName', group_labels{1});
            plot(x_norm, rate_short_i, '-', 'LineWidth', 2, 'Color', group_colors.short, ...
                 'DisplayName', group_labels{2});
            hold off;
            
            xlabel('Normalized track position (start \rightarrow end)');
            ylabel('Firing rate (Hz)');
            legend('show', 'Location', 'best');
            title(sprintf('X-Normalized Comparison (r = %.2f)', r));
            grid on;
        else
            text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', ...
                 'Units', 'normalized', 'FontSize', 12);
            title('X-Normalized Comparison');
        end
    else
        text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', ...
             'Units', 'normalized', 'FontSize', 12);
        title('X-Normalized Comparison');
    end
    
    % Overall figure title
    sgtitle(sprintf('Cluster %d - %s', cluster_id, probe_id), 'FontWeight', 'bold');
end
