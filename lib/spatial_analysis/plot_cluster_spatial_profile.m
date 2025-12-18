function fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, group_names, group_labels, group_colors, probe_id, stats, tuning_data, accel_tuning_data, rate_maps_2d, dist_comparison, ttg_data)
% PLOT_CLUSTER_SPATIAL_PROFILE Create a combined figure for a single cluster
%
%   fig = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, ...
%                                      group_names, group_labels, group_colors, probe_id, stats, tuning_data, accel_tuning_data, rate_maps_2d, dist_comparison, ttg_data)
%
%   Creates a 4-row, 4-column figure showing:
%       Row 1:
%           - Panel 1: Combined position-based raster for both long and short trials
%           - Panel 2: Smoothed traces with median and IQR shading
%           - Panel 3: Shuffle distribution histograms (long trials top, short trials bottom)
%           - Panel 4: X-normalized comparison with percentage x-axis
%       Row 2:
%           - Panel 5: Time-to-goal raster for both long and short trials
%           - Panel 6: Time-to-goal firing rate with quartiles
%           - Panel 7: (empty/reserved)
%           - Panel 8: Normalized time-to-goal comparison
%       Row 3:
%           - Panel 9: Speed tuning curve (combined long+short trials)
%           - Panel 10: Speed tuning shuffle histogram
%           - Panel 11: Acceleration tuning curve (combined long+short trials)
%           - Panel 12: Acceleration tuning shuffle histogram
%       Row 4:
%           - Panel 13: Velocity × Position contour (long trials)
%           - Panel 14: Velocity × Position contour (short trials)
%           - Panel 15: Acceleration × Position contour (long trials)
%           - Panel 16: Acceleration × Position contour (short trials)
%
%   Inputs:
%       cluster_id         - Numeric ID of the cluster
%       bin_centers_by_group - Struct with fields 'long' and 'short' containing bin centers
%       rate_data          - Struct with fields for each group containing:
%                              .Q1_smooth, .Q2_smooth, .Q3_smooth
%                              .rate_per_trial (raw rates, n_trials x n_bins)
%                              .rate_per_trial_smooth (smoothed rates, n_trials x n_bins)
%                              .spike_positions (cell array of spike positions per trial)
%       group_names        - Cell array {'long', 'short'}
%       group_labels       - Cell array with descriptive labels for legend
%       group_colors       - Struct with RGB colors for each group
%       probe_id           - String identifier for the probe (for title)
%       stats              - (Optional) Struct with statistical test results per group
%       tuning_data        - (Optional) Struct with speed tuning curve data (from data.load_tuning_curves)
%       accel_tuning_data  - (Optional) Struct with acceleration tuning curve data (from data.load_tuning_curves_acceleration)
%       rate_maps_2d       - (Optional) Struct with 2D rate maps for position×velocity and position×acceleration
%       dist_comparison    - (Optional) Struct with distribution comparison results (.absolute and .relative)
%       ttg_data           - (Optional) Struct with time-to-goal data per group
%
%   Outputs:
%       fig - Figure handle

    % Handle missing arguments
    if nargin < 13
        ttg_data = [];
    end
    if nargin < 8
        stats = [];
    end
    if nargin < 9
        tuning_data = [];
    end
    if nargin < 12
        dist_comparison = [];
    end
    if nargin < 10
        accel_tuning_data = [];
    end
    if nargin < 11
        rate_maps_2d = [];
    end

    fig = figure('Position', [50, 50, 2400, 1800], 'Visible', 'off');
    
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
    subplot(4, 4, 1);
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
        
        % Define spatial bins for raster (use 2 cm bins from 0-120 cm)
        spatial_bin_size = 2;  % cm
        spatial_bin_edges = 0:spatial_bin_size:120;
        spatial_bin_centers = spatial_bin_edges(1:end-1) + spatial_bin_size/2;
        n_spatial_bins = length(spatial_bin_centers);
        
        % First pass: find global maximum spike count across all trials and bins
        global_max_spike_count = 0;
        for t = 1:length(all_trials_data)
            trial_data = all_trials_data(t);
            if ~isempty(trial_data.spike_positions)
                % Adjust spike positions for short trials (shift to 60-120 cm range)
                spike_pos = trial_data.spike_positions;
                if strcmp(trial_data.group, 'short')
                    spike_pos = spike_pos + 60;
                end
                
                % Bin spike positions for this trial
                spike_counts = histcounts(spike_pos, spatial_bin_edges);
                trial_max = max(spike_counts);
                if trial_max > global_max_spike_count
                    global_max_spike_count = trial_max;
                end
            end
        end
        
        % Second pass: plot trials with normalization across all trials
        for t = 1:length(all_trials_data)
            trial_data = all_trials_data(t);
            if ~isempty(trial_data.spike_positions)
                % Adjust spike positions for short trials (shift to 60-120 cm range)
                spike_pos = trial_data.spike_positions;
                if strcmp(trial_data.group, 'short')
                    spike_pos = spike_pos + 60;
                end
                
                % Bin spike positions for this trial
                spike_counts = histcounts(spike_pos, spatial_bin_edges);
                
                % Plot one vertical bar per bin, with thickness proportional to spike count
                if global_max_spike_count > 0
                    % Normalize spike counts to bar thickness (0 to 0.8 units) using global maximum
                    normalized_counts = spike_counts / global_max_spike_count * 0.8;
                    
                    for b = 1:n_spatial_bins
                        if spike_counts(b) > 0
                            % Bar thickness proportional to spike count (normalized across all trials)
                            bar_half_height = normalized_counts(b) / 2;
                            plot([spatial_bin_centers(b), spatial_bin_centers(b)], ...
                                 [t - bar_half_height, t + bar_half_height], ...
                                 'Color', trial_data.color, 'LineWidth', 2);
                        end
                    end
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
    subplot(4, 4, 2);
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
    
    % Add absolute position distribution comparison result
    if ~isempty(dist_comparison) && isfield(dist_comparison, 'absolute') && isstruct(dist_comparison.absolute)
        if isfield(dist_comparison.absolute, 'same_distribution') && ~isnan(dist_comparison.absolute.same_distribution)
            text_x_pos = 0.02;
            text_y_pos = 0.95;
            
            if dist_comparison.absolute.same_distribution
                abs_str = 'Abs: Same';
            else
                abs_str = 'Abs: Different';
            end
            if isfield(dist_comparison.absolute, 'p_value')
                if dist_comparison.absolute.p_value < 0.001
                    abs_str = sprintf('%s (p<0.001)', abs_str);
                else
                    abs_str = sprintf('%s (p=%.3f)', abs_str, dist_comparison.absolute.p_value);
                end
            end
            text(text_x_pos, text_y_pos, abs_str, 'Units', 'normalized', ...
                 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0 0 0], ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
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
        subplot(4, 4, 3);
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
    subplot(4, 4, 4);
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
    
    % Add relative position distribution comparison result
    if ~isempty(dist_comparison) && isfield(dist_comparison, 'relative') && isstruct(dist_comparison.relative)
        if isfield(dist_comparison.relative, 'same_distribution') && ~isnan(dist_comparison.relative.same_distribution)
            text_x_pos = 0.02;
            text_y_pos = 0.95;
            
            if dist_comparison.relative.same_distribution
                rel_str = 'Rel: Same';
            else
                rel_str = 'Rel: Different';
            end
            if isfield(dist_comparison.relative, 'p_value')
                if dist_comparison.relative.p_value < 0.001
                    rel_str = sprintf('%s (p<0.001)', rel_str);
                else
                    rel_str = sprintf('%s (p=%.3f)', rel_str, dist_comparison.relative.p_value);
                end
            end
            text(text_x_pos, text_y_pos, rel_str, 'Units', 'normalized', ...
                 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0 0 0], ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end
    end
    
    hold off;
    xlabel('Track position (%)');
    ylabel('Firing rate (Hz)');
    grid on;
    xlim([0, 100]);
    ylim([0, max_firing_rate]);
    set(gca, 'XTick', 0:10:100);
    
    %% NEW ROW 2: Time-to-Goal Analysis
    
    %% Panel 5 (row 2, col 1): Normalized time-to-goal raster (0-100%)
    if ~isempty(ttg_data)
        subplot(4, 4, 5);
        hold on;
        
        % Collect all trials with their global indices and group info
        all_trials_ttg = [];
        
        % Collect long trials
        if isfield(ttg_data, 'long') && isfield(ttg_data.long, 'spike_ttg_norm')
            spike_ttg_norm = ttg_data.long.spike_ttg_norm;
            global_indices = ttg_data.long.global_trial_indices;
            n_long_trials = length(spike_ttg_norm);
            color = group_colors.long;
            
            for t = 1:n_long_trials
                all_trials_ttg = [all_trials_ttg; struct(...
                    'global_idx', global_indices(t), ...
                    'spike_ttg_norm', {spike_ttg_norm{t}}, ...
                    'color', color, ...
                    'group', 'long')]; %#ok<AGROW>
            end
        end
        
        % Collect short trials
        if isfield(ttg_data, 'short') && isfield(ttg_data.short, 'spike_ttg_norm')
            spike_ttg_norm = ttg_data.short.spike_ttg_norm;
            global_indices = ttg_data.short.global_trial_indices;
            n_short_trials = length(spike_ttg_norm);
            color = group_colors.short;
            
            for t = 1:n_short_trials
                all_trials_ttg = [all_trials_ttg; struct(...
                    'global_idx', global_indices(t), ...
                    'spike_ttg_norm', {spike_ttg_norm{t}}, ...
                    'color', color, ...
                    'group', 'short')]; %#ok<AGROW>
            end
        end
        
        % Sort by global trial index
        if ~isempty(all_trials_ttg)
            [~, sort_idx] = sort([all_trials_ttg.global_idx]);
            all_trials_ttg = all_trials_ttg(sort_idx);
            
            % Get normalized TTG bin edges from data (0-100%)
            ttg_bin_edges = ttg_data.long.bin_edges_norm;
            ttg_bin_centers = ttg_data.long.bin_centers_norm;
            n_ttg_bins = length(ttg_bin_centers);
            
            % First pass: find global maximum spike count across all trials and bins
            global_max_spike_count = 0;
            for t = 1:length(all_trials_ttg)
                trial_data = all_trials_ttg(t);
                if ~isempty(trial_data.spike_ttg_norm)
                    spike_counts = histcounts(trial_data.spike_ttg_norm, ttg_bin_edges);
                    trial_max = max(spike_counts);
                    if trial_max > global_max_spike_count
                        global_max_spike_count = trial_max;
                    end
                end
            end
            
            % Second pass: plot trials with normalization
            for t = 1:length(all_trials_ttg)
                trial_data = all_trials_ttg(t);
                if ~isempty(trial_data.spike_ttg_norm)
                    % Bin spike TTG (normalized) for this trial
                    spike_counts = histcounts(trial_data.spike_ttg_norm, ttg_bin_edges);
                    
                    % Plot one vertical bar per bin, with thickness proportional to spike count
                    if global_max_spike_count > 0
                        % Normalize spike counts to bar thickness (0 to 0.8 units)
                        normalized_counts = spike_counts / global_max_spike_count * 0.8;
                        
                        for b = 1:n_ttg_bins
                            if spike_counts(b) > 0
                                bar_half_height = normalized_counts(b) / 2;
                                plot([ttg_bin_centers(b), ttg_bin_centers(b)], ...
                                     [t - bar_half_height, t + bar_half_height], ...
                                     'Color', trial_data.color, 'LineWidth', 2);
                            end
                        end
                    end
                end
            end
            
            set(gca, 'YDir', 'reverse', 'YLim', [0, length(all_trials_ttg)+1]);
            if length(all_trials_ttg) > 1
                set(gca, 'YTick', [1, length(all_trials_ttg)]);
            end
        end
        
        hold off;
        xlabel('Time to goal (%)');
        ylabel('Trial # (chronological)');
        grid on;
        xlim([0, 100]);
        set(gca, 'XDir', 'reverse');  % Invert x-axis so 0 is on the right
        set(gca, 'XTick', 0:20:100);
    end
    
    %% Panel 6 (row 2, col 2): Time-to-goal firing rate with NORMALIZED binning (0-100%)
    if ~isempty(ttg_data)
        subplot(4, 4, 6);
        hold on;
        
        % Compute firing rates binned by normalized time-to-goal for each group
        for g = 1:length(group_names)
            group = group_names{g};
            
            if ~isfield(ttg_data, group)
                continue;
            end
            
            spike_ttg_norm = ttg_data.(group).spike_ttg_norm;
            bin_edges = ttg_data.(group).bin_edges_norm;  % 0-100%
            bin_centers = ttg_data.(group).bin_centers_norm;
            n_bins = length(bin_centers);
            n_trials = length(spike_ttg_norm);
            trial_durations = ttg_data.(group).trial_durations;  % Duration of each trial in seconds
            
            % Compute spike counts and occupancy per trial per bin
            spike_counts_per_trial = nan(n_trials, n_bins);
            occupancy_per_trial = nan(n_trials, n_bins);
            
            bin_width_percent = bin_edges(2) - bin_edges(1);  % e.g., 2.5%
            
            for t = 1:n_trials
                if trial_durations(t) > 0
                    % Bin duration in seconds for this trial
                    bin_duration_sec = trial_durations(t) * (bin_width_percent / 100);
                    
                    % Compute spike counts
                    if ~isempty(spike_ttg_norm{t})
                        spike_counts_per_trial(t, :) = histcounts(spike_ttg_norm{t}, bin_edges);
                    else
                        spike_counts_per_trial(t, :) = 0;
                    end
                    
                    % Occupancy: each bin has the same occupancy (bin_duration_sec) since we're
                    % using normalized time - every trial covers 0-100%
                    occupancy_per_trial(t, :) = bin_duration_sec;
                end
            end
            
            % Smooth spike counts and occupancy with Gaussian kernel
            sigma_percent = 13;  % percent
            sigma_bins = sigma_percent / bin_width_percent;
            
            spike_counts_smooth = zeros(size(spike_counts_per_trial));
            occupancy_smooth = zeros(size(occupancy_per_trial));
            
            for t = 1:n_trials
                spike_counts_smooth(t, :) = imgaussfilt(spike_counts_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
                occupancy_smooth(t, :) = imgaussfilt(occupancy_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
            end
            
            % Compute firing rate: smoothed spikes / smoothed occupancy
            rate_per_trial_smooth = spike_counts_smooth ./ occupancy_smooth;
            rate_per_trial_smooth(occupancy_smooth == 0) = nan;  % Avoid division by zero
            
            % Compute quartiles across trials
            Q1_smooth = quantile(rate_per_trial_smooth, 0.25, 1);
            Q2_smooth = quantile(rate_per_trial_smooth, 0.50, 1);  % Median
            Q3_smooth = quantile(rate_per_trial_smooth, 0.75, 1);
            
            color = group_colors.(group);
            
            % Plot Q1-Q3 shaded area
            fill_x = [bin_centers, fliplr(bin_centers)];
            fill_y = [Q1_smooth, fliplr(Q3_smooth)];
            fill(fill_x, fill_y, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
            
            % Plot individual smoothed trials as thin lines
            for t = 1:n_trials
                plot(bin_centers, rate_per_trial_smooth(t, :), '-', 'LineWidth', 0.5, 'Color', [color, 0.3], ...
                     'HandleVisibility', 'off');
            end
            
            % Plot median (Q2) as thick solid line
            plot(bin_centers, Q2_smooth, '-', 'LineWidth', 2.5, 'Color', color, ...
                 'HandleVisibility', 'off');
        end
        
        hold off;
        xlabel('Time to goal (%)');
        ylabel('Firing rate (Hz)');
        grid on;
        xlim([0, 100]);
        set(gca, 'XDir', 'reverse');  % Invert x-axis so 0 is on the right
        set(gca, 'XTick', 0:20:100);
    end
    
    %% Panel 7 (row 2, col 3): Skaggs information for normalized TTG
    if ~isempty(ttg_data)
        subplot(4, 4, 7);
        hold on;
        
        % Determine common x-axis range for both groups
        all_info = [];
        
        for g = 1:length(group_names)
            group = group_names{g};
            
            if ~isfield(ttg_data, group) || ~isfield(ttg_data.(group), 'spike_ttg_norm')
                continue;
            end
            
            spike_ttg_norm = ttg_data.(group).spike_ttg_norm;
            bin_edges = ttg_data.(group).bin_edges_norm;
            bin_centers = ttg_data.(group).bin_centers_norm;
            n_bins = length(bin_centers);
            n_trials = length(spike_ttg_norm);
            
            % Compute rate matrix (trials x bins)
            rate_matrix = nan(n_trials, n_bins);
            for t = 1:n_trials
                if ~isempty(spike_ttg_norm{t})
                    spike_counts = histcounts(spike_ttg_norm{t}, bin_edges);
                    rate_matrix(t, :) = spike_counts;
                end
            end
            
            % Compute mean rate across bins
            mean_rate_per_bin = nanmean(rate_matrix, 1);
            overall_mean_rate = nanmean(mean_rate_per_bin);
            
            % Compute occupancy (fraction of trials in each bin)
            occupancy = sum(~isnan(rate_matrix), 1) / n_trials;
            
            % Compute Skaggs information
            info_observed = 0;
            for b = 1:n_bins
                if occupancy(b) > 0 && mean_rate_per_bin(b) > 0 && overall_mean_rate > 0
                    ratio = mean_rate_per_bin(b) / overall_mean_rate;
                    info_observed = info_observed + occupancy(b) * ratio * log2(ratio);
                end
            end
            
            % Shuffle test (100 shuffles)
            n_shuffles = 100;
            info_shuffled = nan(n_shuffles, 1);
            
            for shuf = 1:n_shuffles
                % Shuffle the rate matrix (break trial-bin associations)
                rate_matrix_shuf = rate_matrix(:);
                rate_matrix_shuf = rate_matrix_shuf(randperm(length(rate_matrix_shuf)));
                rate_matrix_shuf = reshape(rate_matrix_shuf, size(rate_matrix));
                
                % Compute shuffled Skaggs info
                mean_rate_shuf = nanmean(rate_matrix_shuf, 1);
                overall_mean_shuf = nanmean(mean_rate_shuf);
                
                info_shuf = 0;
                for b = 1:n_bins
                    if occupancy(b) > 0 && mean_rate_shuf(b) > 0 && overall_mean_shuf > 0
                        ratio = mean_rate_shuf(b) / overall_mean_shuf;
                        info_shuf = info_shuf + occupancy(b) * ratio * log2(ratio);
                    end
                end
                info_shuffled(shuf) = info_shuf;
            end
            
            % Compute p-value
            p_value = sum(info_shuffled >= info_observed) / n_shuffles;
            
            % Store for plotting
            all_info = [all_info, info_shuffled(:)', info_observed];
            
            % Plot histogram
            color = group_colors.(group);
            [n_hist, ~] = histcounts(info_shuffled, 30);
            edges_hist = linspace(min(info_shuffled), max(info_shuffled), 31);
            bin_centers_hist = (edges_hist(1:end-1) + edges_hist(2:end)) / 2;
            
            bar(bin_centers_hist, n_hist, 'FaceColor', color, 'EdgeColor', 'none', ...
                'FaceAlpha', 0.6, 'BarWidth', 0.8);
            
            % Mark observed value
            y_obs = interp1(bin_centers_hist, n_hist, info_observed, 'linear', 0);
            scatter(info_observed, y_obs, 100, color, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            
            % Add p-value text
            if strcmp(group, 'long')
                text_y = 0.95;
            else
                text_y = 0.85;
            end
            
            if p_value < 0.001
                p_str = sprintf('%s: p < 0.001', group_labels{g});
            else
                p_str = sprintf('%s: p = %.3f', group_labels{g}, p_value);
            end
            text(0.98, text_y, p_str, 'Units', 'normalized', 'FontSize', 9, ...
                 'Color', color, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
        end
        
        hold off;
        xlabel('Skaggs Info (bits/spike)');
        ylabel('Count');
        grid on;
        box off;
    end
    
    %% Panel 8 (row 2, col 4): Empty (reserved for future use)
    % Panels 5-7 in Row 2 show: normalized TTG raster, firing rate, and Skaggs info
    % Panel 8 intentionally left empty - Row 2 uses only 3 columns
    
    %% Panel 9 (row 3, col 1): Speed tuning curve
    if ~isempty(tuning_data) && isstruct(tuning_data)
        ax_tuning = subplot(4, 4, 9);
        
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
    
    %% Panel 10 (row 3, col 2): Speed tuning shuffle histogram
    if ~isempty(tuning_data) && isstruct(tuning_data)
        ax_hist = subplot(4, 4, 10);
        
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
    
    %% Panel 11 (row 3, col 3): Acceleration tuning curve
    if ~isempty(accel_tuning_data) && isstruct(accel_tuning_data)
        ax_accel_tuning = subplot(4, 4, 11);
        
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
    
    %% Panel 12 (row 3, col 4): Acceleration tuning shuffle histogram
    if ~isempty(accel_tuning_data) && isstruct(accel_tuning_data)
        ax_accel_hist = subplot(4, 4, 12);
        
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
    
    %% Row 4: 2D Contour Plots (Velocity × Position and Acceleration × Position)
    if ~isempty(rate_maps_2d)
        % Panel 13: Velocity × Position (long trials) with trial-based transparency
        if isfield(rate_maps_2d, 'vel_long') && ~isempty(rate_maps_2d.vel_long)
            ax_vel_long = subplot(4, 4, 13);
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
            if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
                caxis(ax_vel_long, [cmin, cmax]);
            end
            set(gca, 'Color', [1 1 1]);
            
            xlabel('Position (cm)');
            ylabel('Velocity (cm/s)');
        end
        
        % Panel 14: Velocity × Position (short trials) with trial-based transparency
        if isfield(rate_maps_2d, 'vel_short') && ~isempty(rate_maps_2d.vel_short)
            ax_vel_short = subplot(4, 4, 14);
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
            if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
                caxis(ax_vel_short, [cmin, cmax]);
            end
            set(gca, 'Color', [1 1 1]);
            
            xlabel('Position (cm)');
            ylabel('Velocity (cm/s)');
        end
        
        % Panel 15: Acceleration × Position (long trials) with trial-based transparency
        if isfield(rate_maps_2d, 'accel_long') && ~isempty(rate_maps_2d.accel_long)
            ax_accel_long = subplot(4, 4, 15);
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
            if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
                caxis(ax_accel_long, [cmin, cmax]);
            end
            set(gca, 'Color', [1 1 1]);
            
            xlabel('Position (cm)');
            ylabel('Accel (cm/s^2)');
        end
        
        % Panel 16: Acceleration × Position (short trials) with trial-based transparency
        if isfield(rate_maps_2d, 'accel_short') && ~isempty(rate_maps_2d.accel_short)
            ax_accel_short = subplot(4, 4, 16);
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
            if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
                caxis(ax_accel_short, [cmin, cmax]);
            end
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
