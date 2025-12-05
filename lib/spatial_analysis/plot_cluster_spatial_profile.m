function fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, group_names, group_labels, group_colors, probe_id, stats, tuning_data, accel_tuning_data, rate_maps_2d)
% PLOT_CLUSTER_SPATIAL_PROFILE Create a combined figure for a single cluster
%
%   fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, ...
%                                      group_names, group_labels, group_colors, probe_id, stats, tuning_data, accel_tuning_data, rate_maps_2d)
%
%   Creates a 2-row, 4-column figure showing:
%       Row 1:
%           - Panel 1: Combined position-based raster for both long and short trials
%           - Panel 2: Smoothed traces with median and IQR shading, plus pooled rate (dashed)
%           - Panel 3: Shuffle distribution histograms (long trials top, short trials bottom)
%           - Panel 4: X-normalized comparison with percentage x-axis, plus pooled rate (dashed)
%       Row 2:
%           - Panel 5: Speed tuning curve (combined long+short trials)
%           - Panel 6: Speed tuning shuffle histogram
%           - Panel 7: Acceleration tuning curve (combined long+short trials)
%           - Panel 8: Acceleration tuning shuffle histogram
%       Row 3:
%           - Panel 9: Velocity × Position contour (long trials)
%           - Panel 10: Velocity × Position contour (short trials)
%           - Panel 11: Acceleration × Position contour (long trials)
%           - Panel 12: Acceleration × Position contour (short trials)
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
%       stats              - (Optional) Struct with statistical test results per group
%       tuning_data        - (Optional) Struct with speed tuning curve data (from data.load_tuning_curves)
%       accel_tuning_data  - (Optional) Struct with acceleration tuning curve data (from data.load_tuning_curves_acceleration)
%       rate_maps_2d       - (Optional) Struct with 2D rate maps for position×velocity and position×acceleration
%
%   Outputs:
%       fig - Figure handle

    % Handle missing arguments
    if nargin < 8
        stats = [];
    end
    if nargin < 9
        tuning_data = [];
    end
    if nargin < 10
        accel_tuning_data = [];
    end
    if nargin < 11
        rate_maps_2d = [];
    end

    fig = figure('Position', [50, 50, 2400, 1400], 'Visible', 'off');
    
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
    subplot(3, 4, 1);
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
                % Plot vertical bars instead of dots
                for s = 1:length(trial_data.spike_positions)
                    % Adjust position for short trials (shift by 60 cm)
                    pos = trial_data.spike_positions(s);
                    if strcmp(trial_data.group, 'short')
                        pos = pos + 60;
                    end
                    plot([pos, pos], ...
                         [t-0.4, t+0.4], 'Color', trial_data.color, 'LineWidth', 1);
                end
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
    grid on;
    xlim([0, 120]);
    set(gca, 'XTick', 0:20:120);
    
    %% Panel 2 (col 2): Smoothed traces with median and IQR shading
    subplot(3, 4, 2);
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
    
    % Add asterisks on peaks and field bars for spatially tuned groups
    if ~isempty(stats)
        for g = 1:length(group_names)
            group = group_names{g};
            if isfield(stats, group)
                group_stats = stats.(group);
                if isfield(group_stats, 'isSpatiallyTuned') && group_stats.isSpatiallyTuned
                    if isfield(group_stats, 'peakPosition') && isfield(group_stats, 'peakRate')
                        peak_pos = group_stats.peakPosition;
                        % Adjust peak position for short trials (plotted at 60-120 cm)
                        if strcmp(group, 'short')
                            peak_pos = peak_pos + 60;
                        end
                        peak_rate = group_stats.peakRate;
                        color = group_colors.(group);
                        
                        % Plot large asterisk above the peak
                        text(peak_pos, peak_rate + 0.05*max_firing_rate, '*', ...
                             'FontSize', 24, 'FontWeight', 'bold', 'Color', color, ...
                             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                        
                        % Plot horizontal bar showing the detected field extent
                        if isfield(group_stats, 'fieldStart') && isfield(group_stats, 'fieldEnd')
                            field_start = group_stats.fieldStart;
                            field_end = group_stats.fieldEnd;
                            % Adjust for short trials
                            if strcmp(group, 'short')
                                field_start = field_start + 60;
                                field_end = field_end + 60;
                            end
                            % Draw thick horizontal line with offset so blue/red don't overlap
                            % Short trials (red) at 95%, long trials (blue) at 92%
                            if strcmp(group, 'short')
                                bar_y = 0.95 * max_firing_rate;
                            else
                                bar_y = 0.92 * max_firing_rate;
                            end
                            plot([field_start, field_end], [bar_y, bar_y], '-', ...
                                 'LineWidth', 4, 'Color', color, 'HandleVisibility', 'off');
                        end
                    end
                end
            end
        end
    end
    
    % Add "NS" label if no groups are spatially tuned
    if ~isempty(stats)
        any_significant = false;
        for g = 1:length(group_names)
            group = group_names{g};
            if isfield(stats, group) && isfield(stats.(group), 'isSpatiallyTuned')
                if stats.(group).isSpatiallyTuned
                    any_significant = true;
                    break;
                end
            end
        end
        if ~any_significant
            % Add "NS" in top right corner
            text(0.95, 0.95, 'NS', 'Units', 'normalized', ...
                 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5 0.5 0.5], ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end
    end
    
    hold off;
    xlabel('Position (cm)');
    ylabel('Firing rate (Hz)');
    grid on;
    xlim([0, 120]);
    ylim([0, max_firing_rate]);
    set(gca, 'XTick', 0:20:120);
    
    %% Panel 3 (col 3): Shuffle distribution histograms combined in one vertical plot
    if ~isempty(stats)
        % Combined vertical histogram spanning the full column
        subplot(3, 4, 3);
        hold on;
        
        % Determine common x-axis range
        all_info = [];
        if isfield(stats, 'long') && isfield(stats.long, 'infoNull')
            all_info = [all_info, stats.long.infoNull(:)', stats.long.info];
        end
        if isfield(stats, 'short') && isfield(stats.short, 'infoNull')
            all_info = [all_info, stats.short.infoNull(:)', stats.short.info];
        end
        
        % Remove NaN values before computing edges
        all_info = all_info(~isnan(all_info));
        
        if ~isempty(all_info) && (max(all_info) > min(all_info))
            % Create common bin edges
            edges = linspace(min(all_info), max(all_info), 31);
            bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
            bin_width = edges(2) - edges(1);
            
            % Plot long trials histogram if available
            if isfield(stats, 'long') && isfield(stats.long, 'infoNull') && isfield(stats.long, 'info')
                infoNull_long = stats.long.infoNull;
                infoObs_long = stats.long.info;
                pVal_long = stats.long.pVal;
                
                % Plot histogram with vertical bars
                [n_long, ~] = histcounts(infoNull_long, edges);
                bar(bin_centers, n_long, 'FaceColor', [0.6, 0.8, 1.0], 'EdgeColor', 'none', 'BarWidth', 0.8, 'FaceAlpha', 0.7);
                
                % Mark observed value with dot (interpolate histogram height at observed position)
                y_long = interp1(bin_centers, n_long, infoObs_long, 'linear', 0);
                scatter(infoObs_long, y_long, 100, [0, 0, 0.8], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                
                % Add p-value text for long trials
                if pVal_long < 0.001
                    p_str_long = 'Long: p < 0.001';
                else
                    p_str_long = sprintf('Long: p = %.3f', pVal_long);
                end
                text(0.98, 0.95, p_str_long, 'Units', 'normalized', 'FontSize', 9, 'Color', [0, 0, 0.8], ...
                     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
            end
            
            % Plot short trials histogram if available
            if isfield(stats, 'short') && isfield(stats.short, 'infoNull') && isfield(stats.short, 'info')
                infoNull_short = stats.short.infoNull;
                infoObs_short = stats.short.info;
                pVal_short = stats.short.pVal;
                
                % Plot histogram with vertical bars (overlapping with long)
                [n_short, ~] = histcounts(infoNull_short, edges);
                bar(bin_centers, n_short, 'FaceColor', [1.0, 0.6, 0.6], 'EdgeColor', 'none', 'BarWidth', 0.6, 'FaceAlpha', 0.7);
                
                % Mark observed value with dot (interpolate histogram height at observed position)
                y_short = interp1(bin_centers, n_short, infoObs_short, 'linear', 0);
                scatter(infoObs_short, y_short, 100, [0.8, 0, 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                
                % Add p-value text for short trials
                if pVal_short < 0.001
                    p_str_short = 'Short: p < 0.001';
                else
                    p_str_short = sprintf('Short: p = %.3f', pVal_short);
                end
                text(0.98, 0.85, p_str_short, 'Units', 'normalized', 'FontSize', 9, 'Color', [0.8, 0, 0], ...
                     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
            end
            
            xlabel('Skaggs Info (bits/spike)', 'FontSize', 9);
            ylabel('Count', 'FontSize', 9);
            box off;
            grid on;
        else
            text(0.5, 0.5, 'No shuffle data', 'HorizontalAlignment', 'center', 'FontSize', 10);
            axis off;
        end
        
        hold off;
    end
    
    %% Panel 4 (col 4): X-normalized comparison
    subplot(3, 4, 4);
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
        end
    end
    
    hold off;
    xlabel('Track position (%)');
    ylabel('Firing rate (Hz)');
    grid on;
    xlim([0, 100]);
    ylim([0, max_firing_rate]);
    set(gca, 'XTick', 0:10:100);
    
    %% Panel 5 (row 2, col 1): Speed tuning curve
    if ~isempty(tuning_data) && isstruct(tuning_data)
        ax_tuning = subplot(3, 4, 5);
        
        % Extract tuning curve data
        fr = nanmean(tuning_data.tuning, 2);
        sd = nanstd(tuning_data.tuning, [], 2);
        n = sum(~isnan(tuning_data.tuning), 2);
        x = tuning_data.bin_centers;
        
        % Determine if significant
        is_significant = tuning_data.shuffled.p < 0.05;
        
        % Set colors based on significance
        if is_significant
            data_color = 'k';
            fit_color = 'k';
        else
            data_color = [0.6, 0.6, 0.6];
            fit_color = [0.6, 0.6, 0.6];
        end
        
        % Plot shuffled fits in background (light grey, dashed)
        n_shuffs_to_plot = min(4, size(tuning_data.shuffled.beta_shuff, 1));
        hold on;
        for i = 1:n_shuffs_to_plot
            % Linear polynomial fit for shuffled data
            if length(tuning_data.shuffled.beta_shuff(i, :)) >= 2
                x_fit = linspace(min(x), max(x), 100);
                f = polyval(tuning_data.shuffled.beta_shuff(i, :), x_fit);
                plot(ax_tuning, x_fit, f, '--', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 1.5);
            end
        end
        
        % Plot mean firing rate with error bars (SEM)
        errorbar(ax_tuning, x, fr, sd./sqrt(n), 'o', 'Color', data_color, 'MarkerFaceColor', data_color, 'MarkerSize', 6, 'LineWidth', 1.5);
        
        % Plot true fit (linear polynomial) - black if significant, gray if not
        if isfield(tuning_data.shuffled, 'beta') && length(tuning_data.shuffled.beta) >= 2
            x_fit = linspace(min(x), max(x), 100);
            f = polyval(tuning_data.shuffled.beta, x_fit);
            plot(ax_tuning, x_fit, f, '-', 'Color', fit_color, 'LineWidth', 3.5);
        end
        hold off;
        
        xlabel('Speed (cm/s)');
        ylabel('Firing rate (Hz)');
        grid on;
        xlim([-5, max(x)+5]);
        
        % Add significance indicator
        if is_significant
            % Add asterisk in top-right if significant
            text(0.95, 0.95, '*', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold', ...
                 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        else
            % Add "NS" in top-right if not significant
            text(0.95, 0.95, 'NS', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', ...
                 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end
    end
    
    %% Panel 6 (row 2, col 2): Speed tuning shuffle histogram
    if ~isempty(tuning_data) && isstruct(tuning_data)
        ax_hist = subplot(3, 4, 6);
        
        % Create histogram of shuffled r values
        [n, c] = histcounts(tuning_data.shuffled.r_shuff, 50);
        bin_centers = (c(1:end-1) + c(2:end)) / 2;
        bar(ax_hist, bin_centers, n, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'BarWidth', 1);
        
        hold on;
        % Plot observed r value as a dot
        r_obs = tuning_data.shuffled.r;
        y_at_r = interp1(bin_centers, n, r_obs, 'linear', 0);
        scatter(ax_hist, r_obs, y_at_r, 100, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold off;
        
        xlabel('r');
        ylabel('Count');
        xlim([-0.5, 0.5]);
        set(gca, 'PlotBoxAspectRatio', [2, 1, 1]);
        
        % Add p-value text
        if tuning_data.shuffled.p < 0.001
            p_str = 'p < 0.001';
        else
            p_str = sprintf('p = %.3f', tuning_data.shuffled.p);
        end
        text(0.98, 0.95, p_str, 'Units', 'normalized', 'FontSize', 9, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
    end
    
    %% Panel 7 (row 2, col 3): Acceleration tuning curve
    if ~isempty(accel_tuning_data) && isstruct(accel_tuning_data)
        ax_accel_tuning = subplot(3, 4, 7);
        
        % Extract tuning curve data
        fr = nanmean(accel_tuning_data.tuning, 2);
        sd = nanstd(accel_tuning_data.tuning, [], 2);
        n = sum(~isnan(accel_tuning_data.tuning), 2);
        x = accel_tuning_data.bin_centers;
        
        % Determine if significant
        is_significant = accel_tuning_data.shuffled.p < 0.05;
        
        % Set colors based on significance
        if is_significant
            data_color = 'k';
            fit_color = 'k';
        else
            data_color = [0.6, 0.6, 0.6];
            fit_color = [0.6, 0.6, 0.6];
        end
        
        % Plot shuffled fits in background (light grey, dashed)
        n_shuffs_to_plot = min(4, size(accel_tuning_data.shuffled.beta_shuff, 1));
        hold on;
        for i = 1:n_shuffs_to_plot
            % Linear polynomial fit for shuffled data
            if length(accel_tuning_data.shuffled.beta_shuff(i, :)) >= 2
                x_fit = linspace(min(x), max(x), 100);
                f = polyval(accel_tuning_data.shuffled.beta_shuff(i, :), x_fit);
                plot(ax_accel_tuning, x_fit, f, '--', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 1.5);
            end
        end
        
        % Plot mean firing rate with error bars (SEM)
        errorbar(ax_accel_tuning, x, fr, sd./sqrt(n), 'o', 'Color', data_color, 'MarkerFaceColor', data_color, 'MarkerSize', 6, 'LineWidth', 1.5);
        
        % Plot true fit (linear polynomial) - black if significant, gray if not
        if isfield(accel_tuning_data.shuffled, 'beta') && length(accel_tuning_data.shuffled.beta) >= 2
            x_fit = linspace(min(x), max(x), 100);
            f = polyval(accel_tuning_data.shuffled.beta, x_fit);
            plot(ax_accel_tuning, x_fit, f, '-', 'Color', fit_color, 'LineWidth', 3.5);
        end
        hold off;
        
        xlabel('Acceleration (cm/s^2)');
        ylabel('Firing rate (Hz)');
        grid on;
        xlim([min(x)-5, max(x)+5]);
        
        % Add significance indicator
        if is_significant
            % Add asterisk in top-right if significant
            text(0.95, 0.95, '*', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold', ...
                 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        else
            % Add "NS" in top-right if not significant
            text(0.95, 0.95, 'NS', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', ...
                 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end
    end
    
    %% Panel 8 (row 2, col 4): Acceleration tuning shuffle histogram
    if ~isempty(accel_tuning_data) && isstruct(accel_tuning_data)
        ax_accel_hist = subplot(3, 4, 8);
        
        % Create histogram of shuffled r values
        [n, c] = histcounts(accel_tuning_data.shuffled.r_shuff, 50);
        bin_centers = (c(1:end-1) + c(2:end)) / 2;
        bar(ax_accel_hist, bin_centers, n, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'BarWidth', 1);
        
        hold on;
        % Plot observed r value as a dot
        r_obs = accel_tuning_data.shuffled.r;
        y_at_r = interp1(bin_centers, n, r_obs, 'linear', 0);
        scatter(ax_accel_hist, r_obs, y_at_r, 100, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold off;
        
        xlabel('r');
        ylabel('Count');
        xlim([-0.5, 0.5]);
        set(gca, 'PlotBoxAspectRatio', [2, 1, 1]);
        
        % Add p-value text
        if accel_tuning_data.shuffled.p < 0.001
            p_str = 'p < 0.001';
        else
            p_str = sprintf('p = %.3f', accel_tuning_data.shuffled.p);
        end
        text(0.98, 0.95, p_str, 'Units', 'normalized', 'FontSize', 9, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
    end
    
    %% Row 3: 2D Contour Plots (Velocity × Position and Acceleration × Position)
    if ~isempty(rate_maps_2d)
        % Panel 9: Velocity × Position (long trials) with trial-based transparency
        if isfield(rate_maps_2d, 'vel_long') && ~isempty(rate_maps_2d.vel_long)
            ax_vel_long = subplot(3, 4, 9);
            hold(ax_vel_long, 'on');
            
            % Get bin edges and data
            pos_edges = rate_maps_2d.pos_bin_edges_long;
            vel_edges = rate_maps_2d.vel_bin_edges_long;
            rate_data = rate_maps_2d.vel_long;
            
            % Get colormap limits
            cmin = min(rate_data(:));
            cmax = max(rate_data(:));
            cmap = jet(256);
            
            % Create patches for each bin with variable height
            for p = 1:length(pos_edges)-1
                for v = 1:length(vel_edges)-1
                    rate_val = rate_data(p, v);
                    
                    % Skip NaN values
                    if isnan(rate_val)
                        continue;
                    end
                    
                    % Map rate to color
                    if cmax > cmin
                        color_idx = round(((rate_val - cmin) / (cmax - cmin)) * 255) + 1;
                        color_idx = max(1, min(256, color_idx));
                    else
                        color_idx = 1;
                    end
                    face_color = cmap(color_idx, :);
                    
                    % Determine alpha
                    if isfield(rate_maps_2d, 'vel_trial_counts_long') && isfield(rate_maps_2d, 'vel_n_trials_long')
                        alpha_val = rate_maps_2d.vel_trial_counts_long(p, v) / rate_maps_2d.vel_n_trials_long;
                    else
                        alpha_val = 1;
                    end
                    
                    % Create rectangle patch
                    x = [pos_edges(p), pos_edges(p+1), pos_edges(p+1), pos_edges(p)];
                    y = [vel_edges(v), vel_edges(v), vel_edges(v+1), vel_edges(v+1)];
                    patch(ax_vel_long, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
                end
            end
            
            hold(ax_vel_long, 'off');
            colormap(ax_vel_long, 'jet');
            caxis(ax_vel_long, [cmin, cmax]);
            set(gca, 'Color', [1 1 1]);
            
            xlabel('Position (cm)');
            ylabel('Velocity (cm/s)');
        end
        
        % Panel 10: Velocity × Position (short trials) with trial-based transparency
        if isfield(rate_maps_2d, 'vel_short') && ~isempty(rate_maps_2d.vel_short)
            ax_vel_short = subplot(3, 4, 10);
            hold(ax_vel_short, 'on');
            
            % Get bin edges and data
            pos_edges = rate_maps_2d.pos_bin_edges_short;
            vel_edges = rate_maps_2d.vel_bin_edges_short;
            rate_data = rate_maps_2d.vel_short;
            
            % Get colormap limits
            cmin = min(rate_data(:));
            cmax = max(rate_data(:));
            cmap = jet(256);
            
            % Create patches for each bin with variable height
            for p = 1:length(pos_edges)-1
                for v = 1:length(vel_edges)-1
                    rate_val = rate_data(p, v);
                    
                    % Skip NaN values
                    if isnan(rate_val)
                        continue;
                    end
                    
                    % Map rate to color
                    if cmax > cmin
                        color_idx = round(((rate_val - cmin) / (cmax - cmin)) * 255) + 1;
                        color_idx = max(1, min(256, color_idx));
                    else
                        color_idx = 1;
                    end
                    face_color = cmap(color_idx, :);
                    
                    % Determine alpha
                    if isfield(rate_maps_2d, 'vel_trial_counts_short') && isfield(rate_maps_2d, 'vel_n_trials_short')
                        alpha_val = rate_maps_2d.vel_trial_counts_short(p, v) / rate_maps_2d.vel_n_trials_short;
                    else
                        alpha_val = 1;
                    end
                    
                    % Create rectangle patch
                    x = [pos_edges(p), pos_edges(p+1), pos_edges(p+1), pos_edges(p)];
                    y = [vel_edges(v), vel_edges(v), vel_edges(v+1), vel_edges(v+1)];
                    patch(ax_vel_short, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
                end
            end
            
            hold(ax_vel_short, 'off');
            colormap(ax_vel_short, 'jet');
            caxis(ax_vel_short, [cmin, cmax]);
            set(gca, 'Color', [1 1 1]);
            
            xlabel('Position (cm)');
            ylabel('Velocity (cm/s)');
        end
        
        % Panel 11: Acceleration × Position (long trials) with trial-based transparency
        if isfield(rate_maps_2d, 'accel_long') && ~isempty(rate_maps_2d.accel_long)
            ax_accel_long = subplot(3, 4, 11);
            hold(ax_accel_long, 'on');
            
            % Get bin edges and data
            pos_edges = rate_maps_2d.pos_bin_edges_long;
            accel_edges = rate_maps_2d.accel_bin_edges_long;
            rate_data = rate_maps_2d.accel_long;
            
            % Get colormap limits
            cmin = min(rate_data(:));
            cmax = max(rate_data(:));
            cmap = jet(256);
            
            % Create patches for each bin with variable height
            for p = 1:length(pos_edges)-1
                for a = 1:length(accel_edges)-1
                    rate_val = rate_data(p, a);
                    
                    % Skip NaN values
                    if isnan(rate_val)
                        continue;
                    end
                    
                    % Map rate to color
                    if cmax > cmin
                        color_idx = round(((rate_val - cmin) / (cmax - cmin)) * 255) + 1;
                        color_idx = max(1, min(256, color_idx));
                    else
                        color_idx = 1;
                    end
                    face_color = cmap(color_idx, :);
                    
                    % Determine alpha
                    if isfield(rate_maps_2d, 'accel_trial_counts_long') && isfield(rate_maps_2d, 'accel_n_trials_long')
                        alpha_val = rate_maps_2d.accel_trial_counts_long(p, a) / rate_maps_2d.accel_n_trials_long;
                    else
                        alpha_val = 1;
                    end
                    
                    % Create rectangle patch
                    x = [pos_edges(p), pos_edges(p+1), pos_edges(p+1), pos_edges(p)];
                    y = [accel_edges(a), accel_edges(a), accel_edges(a+1), accel_edges(a+1)];
                    patch(ax_accel_long, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
                end
            end
            
            hold(ax_accel_long, 'off');
            colormap(ax_accel_long, 'jet');
            caxis(ax_accel_long, [cmin, cmax]);
            set(gca, 'Color', [1 1 1]);
            
            xlabel('Position (cm)');
            ylabel('Accel (cm/s^2)');
        end
        
        % Panel 12: Acceleration × Position (short trials) with trial-based transparency
        if isfield(rate_maps_2d, 'accel_short') && ~isempty(rate_maps_2d.accel_short)
            ax_accel_short = subplot(3, 4, 12);
            hold(ax_accel_short, 'on');
            
            % Get bin edges and data
            pos_edges = rate_maps_2d.pos_bin_edges_short;
            accel_edges = rate_maps_2d.accel_bin_edges_short;
            rate_data = rate_maps_2d.accel_short;
            
            % Get colormap limits
            cmin = min(rate_data(:));
            cmax = max(rate_data(:));
            cmap = jet(256);
            
            % Create patches for each bin with variable height
            for p = 1:length(pos_edges)-1
                for a = 1:length(accel_edges)-1
                    rate_val = rate_data(p, a);
                    
                    % Skip NaN values
                    if isnan(rate_val)
                        continue;
                    end
                    
                    % Map rate to color
                    if cmax > cmin
                        color_idx = round(((rate_val - cmin) / (cmax - cmin)) * 255) + 1;
                        color_idx = max(1, min(256, color_idx));
                    else
                        color_idx = 1;
                    end
                    face_color = cmap(color_idx, :);
                    
                    % Determine alpha
                    if isfield(rate_maps_2d, 'accel_trial_counts_short') && isfield(rate_maps_2d, 'accel_n_trials_short')
                        alpha_val = rate_maps_2d.accel_trial_counts_short(p, a) / rate_maps_2d.accel_n_trials_short;
                    else
                        alpha_val = 1;
                    end
                    
                    % Create rectangle patch
                    x = [pos_edges(p), pos_edges(p+1), pos_edges(p+1), pos_edges(p)];
                    y = [accel_edges(a), accel_edges(a), accel_edges(a+1), accel_edges(a+1)];
                    patch(ax_accel_short, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
                end
            end
            
            hold(ax_accel_short, 'off');
            colormap(ax_accel_short, 'jet');
            caxis(ax_accel_short, [cmin, cmax]);
            set(gca, 'Color', [1 1 1]);
            
            xlabel('Position (cm)');
            ylabel('Accel (cm/s^2)');
        end
    end
    
    % Overall figure title (replace underscores with spaces to prevent subscript rendering)
    probe_id_display = strrep(probe_id, '_', ' ');
    
    % Adjust subplot positions to add space at top for title and increase spacing between plots
    % With 3x4 grid, we have better natural spacing, so reduce compression
    all_axes = findall(fig, 'Type', 'axes');
    for ax_idx = 1:length(all_axes)
        ax = all_axes(ax_idx);
        pos = get(ax, 'Position');
        if ~isempty(stats)
            % With statistics: gentle compression with good spacing
            pos(4) = pos(4) * 0.88;  % Minimal compression for better readability
            pos(2) = pos(2) * 0.92 + 0.03;  % Shift up slightly for title space
        else
            % Without statistics: maintain readability with minimal adjustment
            pos(4) = pos(4) * 0.90;  % Minimal compression
            pos(2) = pos(2) * 0.90 + 0.02;  % Shift up slightly
        end
        set(ax, 'Position', pos);
    end
    
    % Create title using annotation positioned WAY at top with large margin below
    annotation('textbox', [0, 0.97, 1, 0.03], 'String', sprintf('Cluster %d - %s', cluster_id, probe_id_display), ...
               'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', ...
               'EdgeColor', 'none', 'VerticalAlignment', 'middle');
    
    %% Bottom area: Statistical test results (if available)
    if ~isempty(stats)
        
        % Create invisible axis at the bottom for text annotation
        ax_text = axes('Position', [0.1, -0.02, 0.8, 0.03], 'Visible', 'off');
        
        % Prepare text content
        text_content = {};
        text_content{end+1} = '\fontsize{7}\bfStatistical Test Results:\rm';
        text_content{end+1} = ' ';
        
        group_fields = fieldnames(stats);
        for g_idx = 1:length(group_fields)
            group = group_fields{g_idx};
            s = stats.(group);
            
            % Find group label
            group_display = group;
            for gl_idx = 1:length(group_names)
                if strcmp(group_names{gl_idx}, group)
                    group_display = group_labels{gl_idx};
                    break;
                end
            end
            
            % Calculate displayed peak position (add 60 for short trials to show absolute position)
            peak_pos_display = s.peakPosition;
            if strcmp(group, 'short')
                peak_pos_display = s.peakPosition + 60;
            end
            
            % Replace underscores with spaces to prevent subscript rendering
            reason_display = strrep(s.reason, '_', ' ');
            
            % Group header
            if s.isSpatiallyTuned
                text_content{end+1} = sprintf('\\bf%s: PASS (spatially tuned)\\rm', group_display);
            else
                text_content{end+1} = sprintf('\\bf%s: FAIL (%s)\\rm', group_display, reason_display);
            end
            
            % Test details (include field size in cm for easier interpretation)
            field_size_cm = s.fieldBins * 2;  % 2cm bins
            text_content{end+1} = sprintf('  Spikes: %d  |  Peak: %.2f Hz @ %.1f cm  |  Info: %.4f (p=%.4f)  |  Field: %d bins (%.0f cm)  |  Num fields: %d', ...
                                         s.n_spikes, s.peakRate, peak_pos_display, s.info, s.pVal, s.fieldBins, field_size_cm, s.numFields);
            text_content{end+1} = ' ';
        end
        
        % Display text in the dedicated text axis at the bottom
        text(ax_text, 0.01, 0.5, text_content, 'FontName', 'FixedWidth', 'FontSize', 6, ...
             'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
             'Interpreter', 'tex');
    end
end
