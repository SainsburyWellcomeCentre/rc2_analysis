function [fig, ttg_fields] = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, group_names, group_labels, group_colors, probe_id, stats, tuning_data, accel_tuning_data, rate_maps_2d, dist_comparison, ttg_data, ttg_stats)
% PLOT_CLUSTER_SPATIAL_PROFILE Create a combined figure for a single cluster
%
%   [fig, ttg_fields] = plot_cluster_spatial_profile(cluster_id, bin_centers_by_group, rate_data, ...
%                                      group_names, group_labels, group_colors, probe_id, stats, tuning_data, accel_tuning_data, rate_maps_2d, dist_comparison, ttg_data, ttg_stats)
%
%   Creates a 5-row, 4-column portrait figure showing:
%       Row 1: Position Tuning
%           - Panel 1: Combined position-based raster for both long and short trials
%           - Panel 2: Smoothed traces with median and IQR shading
%           - Panel 3: Shuffle distribution histograms (long trials top, short trials bottom)
%           - Panel 4: X-normalized comparison with percentage x-axis
%       Row 2: Time-to-Goal Tuning
%           - Panel 5: Time-to-goal raster for both long and short trials
%           - Panel 6: Time-to-goal firing rate with quartiles
%           - Panel 7: TTG information shuffle distribution histograms
%           - Panel 8: Kruskal-Wallis test summary
%       Row 3: Velocity & Acceleration Tuning Curves
%           - Panel 9: Velocity tuning curve (combined long+short)
%           - Panel 10: Velocity shuffle histogram
%           - Panel 11: Acceleration tuning curve (combined long+short)
%           - Panel 12: Acceleration shuffle histogram
%       Row 4: Position × Velocity/Acceleration Maps
%           - Panel 13: Velocity × Position (long trials)
%           - Panel 14: Velocity × Position (short trials)
%           - Panel 15: Acceleration × Position (long trials)
%           - Panel 16: Acceleration × Position (short trials)
%       Row 5: TTG × Velocity/Acceleration Maps
%           - Panel 17: Velocity × TTG (long trials)
%           - Panel 18: Velocity × TTG (short trials)
%           - Panel 19: Acceleration × TTG (long trials)
%           - Panel 20: Acceleration × TTG (short trials)
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
%       stats              - (Optional) Struct with spatial tuning statistical test results per group
%       tuning_data        - (Optional) Struct with speed tuning curve data (from data.load_tuning_curves)
%       accel_tuning_data  - (Optional) Struct with acceleration tuning curve data (from data.load_tuning_curves_acceleration)
%       rate_maps_2d       - (Optional) Struct with 2D rate maps for position×velocity and position×acceleration
%       dist_comparison    - (Optional) Struct with distribution comparison results (.absolute and .relative)
%       ttg_data           - (Optional) Struct with time-to-goal data per group
%       ttg_stats          - (Optional) Struct with TTG tuning statistical test results per group
%
%   Outputs:
%       fig - Figure handle
%       ttg_fields - Struct with TTG field information per group:
%                    .long.n_fields - number of detected fields
%                    .long.fields - cell array of field bin indices
%                    .short.n_fields - number of detected fields  
%                    .short.fields - cell array of field bin indices

    % Initialize ttg_fields output
    ttg_fields = struct();
    
    % Handle missing arguments
    if nargin < 14
        ttg_stats = [];
    end
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

    fig = figure('Position', [50, 50, 2000, 2600], 'Visible', 'off');
    
    % Set tighter subplot spacing
    set(fig, 'DefaultAxesFontSize', 8);  % Smaller default font size
    
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
    subplot(5, 4, 1, 'Parent', fig);
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
        
        % Define spatial bins for heatmap (use 2 cm bins from 0-120 cm)
        spatial_bin_size = 2;  % cm
        spatial_bin_edges = 0:spatial_bin_size:120;
        spatial_bin_centers = spatial_bin_edges(1:end-1) + spatial_bin_size/2;
        n_spatial_bins = length(spatial_bin_centers);
        n_trials = length(all_trials_data);
        
        % Build spike count matrix (trials x spatial bins)
        spike_count_matrix = zeros(n_trials, n_spatial_bins);
        trial_colors = zeros(n_trials, 3);  % Store colors for each trial
        
        for t = 1:n_trials
            trial_data = all_trials_data(t);
            trial_colors(t, :) = trial_data.color;  % Store color
            
            if ~isempty(trial_data.spike_positions)
                % Adjust spike positions for short trials (shift to 60-120 cm range)
                spike_pos = trial_data.spike_positions;
                if strcmp(trial_data.group, 'short')
                    spike_pos = spike_pos + 60;
                end
                
                % Bin spike positions for this trial
                spike_count_matrix(t, :) = histcounts(spike_pos, spatial_bin_edges);
            end
        end
        
        % Create RGB image where each row uses its trial color (blue or red)
        % Normalize spike counts using GLOBAL maximum across all trials and bins
        % (same normalization as the original histogram)
        global_max_spike_count = max(spike_count_matrix(:));
        if global_max_spike_count == 0
            global_max_spike_count = 1;  % Avoid division by zero
        end
        normalized_matrix = spike_count_matrix / global_max_spike_count;
        
        % Create RGB image: each row blends from white (0 spikes) to trial color (max spikes)
        rgb_image = ones(n_trials, n_spatial_bins, 3);
        for t = 1:n_trials
            for c = 1:3  % RGB channels
                % Blend from white (1) to trial color based on normalized spike count
                rgb_image(t, :, c) = 1 - normalized_matrix(t, :) * (1 - trial_colors(t, c));
            end
        end
        
        % Display RGB image
        image(spatial_bin_centers, 1:n_trials, rgb_image);
        
        set(gca, 'YDir', 'reverse', 'YLim', [0.5, n_trials+0.5]);
        if n_trials > 1
            set(gca, 'YTick', [1, n_trials]);
        end
    end
    
    hold off;
    xlabel('Position (cm)');
    ylabel('Trial # (chronological)');
    grid on;
    xlim([0, 120]);
    set(gca, 'XTick', 0:20:120);
    
    %% Panel 2 (col 2): Smoothed traces with median and IQR shading
    subplot(5, 4, 2, 'Parent', fig);
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
    
    hold off;
    xlabel('Position (cm)');
    ylabel('Firing rate (Hz)');
    grid on;
    xlim([0, 120]);
    ylim([0, max_firing_rate]);
    set(gca, 'XTick', 0:20:120);
    
    %% Panel 3 (row 1, col 3): Shuffle distribution histograms for SPATIAL tuning (Skaggs info)
    if ~isempty(stats)
        % Combined vertical histogram spanning the full column
        subplot(5, 4, 3, 'Parent', fig);
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
    
    %% Panel 4 (row 1, col 4): X-normalized comparison
    subplot(5, 4, 4, 'Parent', fig);
    hold on;
    
    % Common normalized grid (0 to 1, displayed as 0% to 100%)
    x_norm = linspace(0, 1, 100);
    x_percent = x_norm * 100;  % For display as percentage
    
    % Plot downsampled long data with individual trials
    if isfield(rate_data, 'long_downsampled')
        group = 'long_downsampled';
        pos = rate_data.(group).bin_centers(:);
        Q1_smooth = rate_data.(group).Q1_smooth(:);
        Q2_smooth = rate_data.(group).Q2_smooth(:);
        Q3_smooth = rate_data.(group).Q3_smooth(:);
        rate_per_trial_smooth = rate_data.(group).rate_per_trial_smooth;
        
        % Skip if all NaN
        if ~all(isnan(Q2_smooth))
            % Ensure pos and Q2_smooth have same length
            n_bins_data = length(Q2_smooth);
            if length(pos) ~= n_bins_data
                % Trim or pad pos to match Q2_smooth
                if length(pos) > n_bins_data
                    pos = pos(1:n_bins_data);
                else
                    % Pad with extrapolated values
                    pos_step = pos(end) - pos(end-1);
                    n_pad = n_bins_data - length(pos);
                    pos = [pos; (pos(end) + (1:n_pad)' * pos_step)];
                end
            end
            
            % Normalize position (start=0, end=1)
            pos_norm = (pos - min(pos)) / (max(pos) - min(pos));
            
            % Keep only finite points
            idx = isfinite(pos_norm) & isfinite(Q2_smooth);
            
            if sum(idx) >= 2
                pos_n = pos_norm(idx);
                Q1_n = Q1_smooth(idx);
                Q2_n = Q2_smooth(idx);
                Q3_n = Q3_smooth(idx);
                
                % Interpolate to common grid
                Q1_i = interp1(pos_n, Q1_n, x_norm, 'linear', 'extrap');
                Q2_i = interp1(pos_n, Q2_n, x_norm, 'linear', 'extrap');
                Q3_i = interp1(pos_n, Q3_n, x_norm, 'linear', 'extrap');
                
                color = group_colors.long;  % Use long color (blue)
                
                % Plot individual trials as thin lines
                n_trials = size(rate_per_trial_smooth, 1);
                n_bins_trial = size(rate_per_trial_smooth, 2);
                for t = 1:n_trials
                    trial_rate = rate_per_trial_smooth(t, :);
                    % Create position array matching trial_rate dimensions
                    if length(pos_norm) == n_bins_trial
                        pos_norm_trial = pos_norm(:)';  % Ensure row vector
                    else
                        % If mismatch, use first n_bins_trial elements or pad
                        pos_norm_trial = pos_norm(1:min(n_bins_trial, length(pos_norm)));
                        pos_norm_trial = pos_norm_trial(:)';  % Ensure row vector
                        if length(pos_norm_trial) < n_bins_trial
                            % Pad if needed
                            pos_norm_trial = [pos_norm_trial, nan(1, n_bins_trial - length(pos_norm_trial))];
                        end
                    end
                    idx_trial = isfinite(pos_norm_trial) & isfinite(trial_rate);
                    if sum(idx_trial) >= 2
                        trial_rate_i = interp1(pos_norm_trial(idx_trial), trial_rate(idx_trial), x_norm, 'linear', 'extrap');
                        plot(x_percent, trial_rate_i, '-', 'LineWidth', 0.5, 'Color', [color, 0.3], ...
                             'HandleVisibility', 'off');
                    end
                end
                
                % Plot Q1-Q3 shaded area
                fill_x = [x_percent, fliplr(x_percent)];
                fill_y = [Q1_i, fliplr(Q3_i)];
                fill(fill_x, fill_y, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                     'HandleVisibility', 'off');
                
                % Plot median as thick solid line
                plot(x_percent, Q2_i, '-', 'LineWidth', 2, 'Color', color, ...
                     'HandleVisibility', 'off');
            end
        end
    end
    
    % Plot short data with individual trials
    if isfield(rate_data, 'short')
        group = 'short';
        pos = bin_centers_by_group.(group)(:);
        Q1_smooth = rate_data.(group).Q1_smooth(:);
        Q2_smooth = rate_data.(group).Q2_smooth(:);
        Q3_smooth = rate_data.(group).Q3_smooth(:);
        rate_per_trial_smooth = rate_data.(group).rate_per_trial_smooth;
        
        % Skip if all NaN
        if ~all(isnan(Q2_smooth))
            % Ensure pos and Q2_smooth have same length
            n_bins_data = length(Q2_smooth);
            if length(pos) ~= n_bins_data
                % Trim or pad pos to match Q2_smooth
                if length(pos) > n_bins_data
                    pos = pos(1:n_bins_data);
                else
                    % Pad with extrapolated values
                    pos_step = pos(end) - pos(end-1);
                    n_pad = n_bins_data - length(pos);
                    pos = [pos; (pos(end) + (1:n_pad)' * pos_step)];
                end
            end
            
            % Normalize position (start=0, end=1)
            pos_norm = (pos - min(pos)) / (max(pos) - min(pos));
            
            % Keep only finite points
            idx = isfinite(pos_norm) & isfinite(Q2_smooth);
            
            if sum(idx) >= 2
                pos_n = pos_norm(idx);
                Q1_n = Q1_smooth(idx);
                Q2_n = Q2_smooth(idx);
                Q3_n = Q3_smooth(idx);
                
                % Interpolate to common grid
                Q1_i = interp1(pos_n, Q1_n, x_norm, 'linear', 'extrap');
                Q2_i = interp1(pos_n, Q2_n, x_norm, 'linear', 'extrap');
                Q3_i = interp1(pos_n, Q3_n, x_norm, 'linear', 'extrap');
                
                color = group_colors.(group);
                
                % Plot individual trials as thin lines
                n_trials = size(rate_per_trial_smooth, 1);
                n_bins_trial = size(rate_per_trial_smooth, 2);
                for t = 1:n_trials
                    trial_rate = rate_per_trial_smooth(t, :);
                    % Create position array matching trial_rate dimensions
                    if length(pos_norm) == n_bins_trial
                        pos_norm_trial = pos_norm(:)';  % Ensure row vector
                    else
                        % If mismatch, use first n_bins_trial elements or pad
                        pos_norm_trial = pos_norm(1:min(n_bins_trial, length(pos_norm)));
                        pos_norm_trial = pos_norm_trial(:)';  % Ensure row vector
                        if length(pos_norm_trial) < n_bins_trial
                            % Pad if needed
                            pos_norm_trial = [pos_norm_trial, nan(1, n_bins_trial - length(pos_norm_trial))];
                        end
                    end
                    idx_trial = isfinite(pos_norm_trial) & isfinite(trial_rate);
                    if sum(idx_trial) >= 2
                        trial_rate_i = interp1(pos_norm_trial(idx_trial), trial_rate(idx_trial), x_norm, 'linear', 'extrap');
                        plot(x_percent, trial_rate_i, '-', 'LineWidth', 0.5, 'Color', [color, 0.3], ...
                             'HandleVisibility', 'off');
                    end
                end
                
                % Plot Q1-Q3 shaded area
                fill_x = [x_percent, fliplr(x_percent)];
                fill_y = [Q1_i, fliplr(Q3_i)];
                fill(fill_x, fill_y, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                     'HandleVisibility', 'off');
                
                % Plot median as thick solid line
                plot(x_percent, Q2_i, '-', 'LineWidth', 2, 'Color', color, ...
                     'HandleVisibility', 'off');
            end
        end
    end
    
    % Compute correlation between downsampled long and short medians
    if isfield(rate_data, 'long_downsampled') && isfield(rate_data, 'short')
        pos_long = rate_data.long_downsampled.bin_centers(:);
        Q2_long = rate_data.long_downsampled.Q2_smooth(:);
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
    
    %% NEW ROW 2: Time-to-Goal Analysis
    
    %% Panel 5 (row 2, col 1): Normalized time-to-goal raster (0-100%)
    if ~isempty(ttg_data)
        subplot(5, 4, 5, 'Parent', fig);
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
            n_trials = length(all_trials_ttg);
            
            % Build spike count matrix (trials x TTG bins)
            spike_count_matrix = zeros(n_trials, n_ttg_bins);
            trial_colors = zeros(n_trials, 3);  % Store colors for each trial
            
            for t = 1:n_trials
                trial_data = all_trials_ttg(t);
                trial_colors(t, :) = trial_data.color;  % Store color
                
                if ~isempty(trial_data.spike_ttg_norm)
                    % Bin spike TTG (normalized) for this trial
                    spike_count_matrix(t, :) = histcounts(trial_data.spike_ttg_norm, ttg_bin_edges);
                end
            end
            
            % Create RGB image where each row uses its trial color (blue or red)
            % Normalize spike counts using GLOBAL maximum across all trials and bins
            global_max_spike_count = max(spike_count_matrix(:));
            if global_max_spike_count == 0
                global_max_spike_count = 1;  % Avoid division by zero
            end
            normalized_matrix = spike_count_matrix / global_max_spike_count;
            
            % Create RGB image: each row blends from white (0 spikes) to trial color (max spikes)
            rgb_image = ones(n_trials, n_ttg_bins, 3);
            for t = 1:n_trials
                for c = 1:3  % RGB channels
                    % Blend from white (1) to trial color based on normalized spike count
                    rgb_image(t, :, c) = 1 - normalized_matrix(t, :) * (1 - trial_colors(t, c));
                end
            end
            
            % Display RGB image
            image(ttg_bin_centers, 1:n_trials, rgb_image);
            
            set(gca, 'YDir', 'reverse', 'YLim', [0.5, n_trials+0.5]);
            if n_trials > 1
                set(gca, 'YTick', [1, n_trials]);
            end
        end
        
        hold off;
        xlabel('Time to goal (%)');
        ylabel('Trial # (chronological)');
        grid on;
        xlim([0, 100]);
        set(gca, 'XDir', 'reverse');  % 100% (far from goal) on left, 0% (at goal) on right
        set(gca, 'XTick', 0:20:100);
    end
    
    %% Panel 6 (row 2, col 2): Time-to-goal firing rate with NORMALIZED binning (0-100%)
    if ~isempty(ttg_data)
        subplot(5, 4, 6, 'Parent', fig);
        hold on;
        
        % Store Q3 values to compute y-axis max
        all_Q3_ttg = [];
        
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
            sigma_percent = 8;  % percent (4 bins * 2% per bin)
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
            
            % Store Q3 for y-axis calculation
            all_Q3_ttg = [all_Q3_ttg, Q3_smooth];
            
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
        
        % Determine Y-axis limit first to ensure fields are plotted correctly
        y_max_limit = 10; % Default
        if ~isempty(all_Q3_ttg)
            max_Q3_ttg = max(all_Q3_ttg);
            if ~isnan(max_Q3_ttg) && max_Q3_ttg > 0
                y_max_limit = max_Q3_ttg * 1.2;
            end
        end
        ylim([0, y_max_limit]);

        % Compute and mark firing fields for each group (simple peak detection)
        % Use the same field detection logic as spatial tuning
        % Threshold: 70% of the range (min to max) of the smoothed median firing rate
        % Minimum width: 4 contiguous bins (approx 10% of trial duration)
        fieldFrac = 0.7;  
        minFieldBins = 4;  
        
        for g = 1:length(group_names)
            group = group_names{g};
            
            % Initialize field info for this group
            ttg_fields.(group).n_fields = 0;
            ttg_fields.(group).fields = {};
            
            if ~isfield(ttg_data, group)
                continue;
            end
            
            % Get data for this group
            spike_ttg_norm = ttg_data.(group).spike_ttg_norm;
            bin_edges = ttg_data.(group).bin_edges_norm;
            bin_centers = ttg_data.(group).bin_centers_norm;
            n_bins = length(bin_centers);
            n_trials = length(spike_ttg_norm);
            trial_durations = ttg_data.(group).trial_durations;
            
            % Recompute smoothed median rate (Q2) for field detection
            spike_counts_per_trial = nan(n_trials, n_bins);
            occupancy_per_trial = nan(n_trials, n_bins);
            bin_width_percent = bin_edges(2) - bin_edges(1);
            
            for t = 1:n_trials
                if trial_durations(t) > 0
                    bin_duration_sec = trial_durations(t) * (bin_width_percent / 100);
                    if ~isempty(spike_ttg_norm{t})
                        spike_counts_per_trial(t, :) = histcounts(spike_ttg_norm{t}, bin_edges);
                    else
                        spike_counts_per_trial(t, :) = 0;
                    end
                    occupancy_per_trial(t, :) = bin_duration_sec;
                end
            end
            
            % Smooth
            sigma_percent = 8;
            sigma_bins = sigma_percent / bin_width_percent;
            spike_counts_smooth = zeros(size(spike_counts_per_trial));
            occupancy_smooth = zeros(size(occupancy_per_trial));
            for t = 1:n_trials
                spike_counts_smooth(t, :) = imgaussfilt(spike_counts_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
                occupancy_smooth(t, :) = imgaussfilt(occupancy_per_trial(t, :), sigma_bins, 'Padding', 'replicate');
            end
            rate_per_trial_smooth = spike_counts_smooth ./ occupancy_smooth;
            rate_per_trial_smooth(occupancy_smooth == 0) = nan;
            Q2_smooth = quantile(rate_per_trial_smooth, 0.50, 1);
            
            % Detect all fields using threshold-based detection
            minRate = min(Q2_smooth);
            maxRate = max(Q2_smooth);
            thresh = minRate + fieldFrac * (maxRate - minRate);
            above = find(Q2_smooth >= thresh);
            
            if ~isempty(above)
                % Find contiguous segments
                d = diff(above);
                breaks = [0, find(d > 1), length(above)];
                numFields = length(breaks) - 1;
                
                % Filter fields: only keep those with at least minFieldBins contiguous bins
                valid_fields = [];
                for k = 1:numFields
                    idx = (breaks(k)+1) : breaks(k+1);
                    segBins = above(idx);
                    
                    if length(segBins) >= minFieldBins
                        valid_fields = [valid_fields; k]; %#ok<AGROW>
                    end
                end
                
                % Store and plot valid fields
                ttg_fields.(group).n_fields = length(valid_fields);
                
                for i = 1:length(valid_fields)
                    k = valid_fields(i);
                    idx = (breaks(k)+1) : breaks(k+1);
                    segBins = above(idx);
                    
                    % Store field bin indices
                    ttg_fields.(group).fields{i} = segBins;
                    
                    % Get start and end positions
                    field_start = bin_centers(segBins(1));
                    field_end = bin_centers(segBins(end));
                    
                    % Draw horizontal line for this field with condition-specific color
                    % Offset: short trials at 95%, long trials at 92%
                    if strcmp(group, 'short')
                        bar_y_frac = 0.95;
                        field_color = [0.8, 0, 0];  % Red for short trials
                    else
                        bar_y_frac = 0.92;
                        field_color = [0, 0, 0.8];  % Blue for long trials
                    end
                    
                    % Use the fixed y_max_limit for consistent plotting
                    bar_y = bar_y_frac * y_max_limit;
                    
                    % Plot with condition-specific color
                    plot([field_start, field_end], [bar_y, bar_y], '-', ...
                         'LineWidth', 4, 'Color', field_color, 'HandleVisibility', 'off');
                end
            end
        end
        
        hold off;
        xlabel('Time to goal (%)');
        ylabel('Firing rate (Hz)');
        grid on;
        xlim([0, 100]);
        set(gca, 'XDir', 'reverse');  % 100% (far from goal) on left, 0% (at goal) on right
        set(gca, 'XTick', 0:20:100);
    end
    
    %% Panel 7 (row 2, col 3): TTG information shuffle distribution histograms
    if ~isempty(ttg_stats)
        subplot(5, 4, 7, 'Parent', fig);
        hold on;
        
        % Collect all info values to determine common x-axis
        all_info = [];
        if isfield(ttg_stats, 'long') && isfield(ttg_stats.long, 'infoNull')
            all_info = [all_info, ttg_stats.long.infoNull(:)', ttg_stats.long.info];
        end
        if isfield(ttg_stats, 'short') && isfield(ttg_stats.short, 'infoNull')
            all_info = [all_info, ttg_stats.short.infoNull(:)', ttg_stats.short.info];
        end
        
        % Remove NaN values before computing edges
        all_info = all_info(~isnan(all_info));
        
        if ~isempty(all_info) && (max(all_info) > min(all_info))
            % Create common bin edges
            edges = linspace(min(all_info), max(all_info), 31);
            bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
            
            % Plot long trials histogram if available
            if isfield(ttg_stats, 'long') && isfield(ttg_stats.long, 'infoNull') && isfield(ttg_stats.long, 'info')
                infoNull_long = ttg_stats.long.infoNull;
                infoObs_long = ttg_stats.long.info;
                pVal_long = ttg_stats.long.pVal;
                
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
            if isfield(ttg_stats, 'short') && isfield(ttg_stats.short, 'infoNull') && isfield(ttg_stats.short, 'info')
                infoNull_short = ttg_stats.short.infoNull;
                infoObs_short = ttg_stats.short.info;
                pVal_short = ttg_stats.short.pVal;
                
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
            
            xlabel('TTG Info (bits/spike)', 'FontSize', 9);
            ylabel('Count', 'FontSize', 9);
            box off;
            grid on;
        else
            text(0.5, 0.5, 'No shuffle data', 'HorizontalAlignment', 'center', 'FontSize', 10);
            axis off;
        end
        
        hold off;
    end
    
    %% Panel 8 (row 2, col 4): Kruskal-Wallis test results summary
    if ~isempty(dist_comparison) && isstruct(dist_comparison)
        ax_text = subplot(5, 4, 8, 'Parent', fig);
        axis(ax_text, 'off');
        
        % Prepare text for each test
        test_names = {'Absolute (60-120cm)', 'Relative (normalized)', 'Absolute Shifted (0-60cm)', 'TTG'};
        test_fields = {'absolute', 'relative', 'absolute_shifted', 'ttg'};
        
        % Starting y position
        y_pos = 0.95;
        line_spacing = 0.20;
        
        % Title
        text(0.5, y_pos, 'Kruskal-Wallis Tests', 'Units', 'normalized', ...
            'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        y_pos = y_pos - line_spacing;
        
        % Display each test result
        for i = 1:length(test_fields)
            field = test_fields{i};
            if isfield(dist_comparison, field) && isstruct(dist_comparison.(field))
                result = dist_comparison.(field);
                
                % Check for valid result
                if isfield(result, 'p_value') && ~isnan(result.p_value)
                    % Determine if distributions are same or different
                    if result.same_distribution
                        result_str = 'same';
                        color = [0.8, 0, 0];  % Red for same (user cares about this)
                    else
                        result_str = 'diff';
                        color = [0.5, 0.5, 0.5];  % Gray for different
                    end
                    
                    % Display test name and result
                    text(0.05, y_pos, test_names{i}, 'Units', 'normalized', ...
                        'FontSize', 7.5, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
                    text(0.95, y_pos, result_str, 'Units', 'normalized', ...
                        'FontSize', 7.5, 'HorizontalAlignment', 'right', 'Color', color, 'FontWeight', 'bold');
                else
                    % No valid result
                    text(0.05, y_pos, test_names{i}, 'Units', 'normalized', ...
                        'FontSize', 7.5, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
                    text(0.95, y_pos, 'N/A', 'Units', 'normalized', ...
                        'FontSize', 7.5, 'HorizontalAlignment', 'right', 'Color', [0.5, 0.5, 0.5]);
                end
            else
                % Field doesn't exist
                text(0.05, y_pos, test_names{i}, 'Units', 'normalized', ...
                    'FontSize', 7.5, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
                text(0.95, y_pos, 'N/A', 'Units', 'normalized', ...
                    'FontSize', 7.5, 'HorizontalAlignment', 'right', 'Color', [0.5, 0.5, 0.5]);
            end
            
            y_pos = y_pos - line_spacing;
        end
    end
    
    %% Panel 9 (row 3, col 1): Velocity tuning curve
    if ~isempty(tuning_data) && isstruct(tuning_data)
        ax_tuning = subplot(5, 4, 9, 'Parent', fig);
        
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
            % Quadratic polynomial fit for shuffled data
            if length(tuning_data.shuffled.beta_shuff(i, :)) >= 3
                x_fit = linspace(min(x), max(x), 100);
                f = polyval(tuning_data.shuffled.beta_shuff(i, :), x_fit);
                plot(ax_tuning, x_fit, f, '--', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 1.5);
            end
        end
        
        % Plot mean firing rate with error bars (SEM)
        errorbar(ax_tuning, x, fr, sd./sqrt(n), 'o', 'Color', data_color, 'MarkerFaceColor', data_color, 'MarkerSize', 6, 'LineWidth', 1.5);
        
        % Plot true fit (quadratic polynomial) - black if significant, gray if not
        if isfield(tuning_data.shuffled, 'beta') && length(tuning_data.shuffled.beta) >= 3
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
        
        % Add beta parameter values as formula
        if isfield(tuning_data.shuffled, 'beta') && length(tuning_data.shuffled.beta) >= 3
            beta = tuning_data.shuffled.beta;
            beta_str = sprintf('%.2fx^2 %+.2fx %+.2f', beta(1), beta(2), beta(3));
            text(0.05, 0.05, beta_str, 'Units', 'normalized', 'FontSize', 7, ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontWeight', 'normal', ...
                 'BackgroundColor', [1 1 1 0.7], 'Interpreter', 'tex');
        end
    end
    
    %% Panel 10 (row 3, col 2): Velocity tuning shuffle histogram
    if ~isempty(tuning_data) && isstruct(tuning_data)
        ax_hist = subplot(5, 4, 10, 'Parent', fig);
        
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
        ax_accel_tuning = subplot(5, 4, 11, 'Parent', fig);
        
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
            % Quadratic polynomial fit for shuffled data
            if length(accel_tuning_data.shuffled.beta_shuff(i, :)) >= 3
                x_fit = linspace(min(x), max(x), 100);
                f = polyval(accel_tuning_data.shuffled.beta_shuff(i, :), x_fit);
                plot(ax_accel_tuning, x_fit, f, '--', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 1.5);
            end
        end
        
        % Plot mean firing rate with error bars (SEM)
        errorbar(ax_accel_tuning, x, fr, sd./sqrt(n), 'o', 'Color', data_color, 'MarkerFaceColor', data_color, 'MarkerSize', 6, 'LineWidth', 1.5);
        
        % Plot true fit (quadratic polynomial) - black if significant, gray if not
        if isfield(accel_tuning_data.shuffled, 'beta') && length(accel_tuning_data.shuffled.beta) >= 3
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
        
        % Add beta parameter values as formula
        if isfield(accel_tuning_data.shuffled, 'beta') && length(accel_tuning_data.shuffled.beta) >= 3
            beta = accel_tuning_data.shuffled.beta;
            beta_str = sprintf('%.2fx^2 %+.2fx %+.2f', beta(1), beta(2), beta(3));
            text(0.05, 0.05, beta_str, 'Units', 'normalized', 'FontSize', 7, ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontWeight', 'normal', ...
                 'BackgroundColor', [1 1 1 0.7], 'Interpreter', 'tex');
        end
    end
    
    %% Panel 12 (row 3, col 4): Acceleration tuning shuffle histogram
    if ~isempty(accel_tuning_data) && isstruct(accel_tuning_data)
        ax_accel_hist = subplot(5, 4, 12, 'Parent', fig);
        
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
    
    %% Panel 13 (row 4, col 1): Velocity × Position contour (long trials)
    if isfield(rate_maps_2d, 'vel_long') && ~isempty(rate_maps_2d.vel_long)
        ax_vel_long = subplot(5, 4, 13, 'Parent', fig);
        hold(ax_vel_long, 'on');
        
        % Get bin edges and data
        pos_edges = rate_maps_2d.pos_bin_edges_long;
        vel_edges = rate_maps_2d.vel_bin_edges_long;
        rate_data = rate_maps_2d.vel_long;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
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
        colormap(ax_vel_long, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_vel_long, [cmin, cmax]);
            colorbar(ax_vel_long);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Position (cm)');
        ylabel('Velocity (cm/s)');
        title('Vel x Pos (Long)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end
    
    %% Panel 14 (row 4, col 2): Velocity × Position contour (short trials)
    if isfield(rate_maps_2d, 'vel_short') && ~isempty(rate_maps_2d.vel_short)
        ax_vel_short = subplot(5, 4, 14, 'Parent', fig);
        hold(ax_vel_short, 'on');
        
        % Get bin edges and data
        pos_edges = rate_maps_2d.pos_bin_edges_short;
        vel_edges = rate_maps_2d.vel_bin_edges_short;
        rate_data = rate_maps_2d.vel_short;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
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
        colormap(ax_vel_short, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_vel_short, [cmin, cmax]);
            colorbar(ax_vel_short);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Position (cm)');
        ylabel('Velocity (cm/s)');
        title('Vel x Pos (Short)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end
    
    %% Panel 15 (row 4, col 3): Acceleration × Position contour (long trials)
    if isfield(rate_maps_2d, 'accel_long') && ~isempty(rate_maps_2d.accel_long)
        ax_accel_long = subplot(5, 4, 15, 'Parent', fig);
        hold(ax_accel_long, 'on');
        
        % Get bin edges and data
        pos_edges = rate_maps_2d.pos_bin_edges_long;
        accel_edges = rate_maps_2d.accel_bin_edges_long;
        rate_data = rate_maps_2d.accel_long;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
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
        colormap(ax_accel_long, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_accel_long, [cmin, cmax]);
            colorbar(ax_accel_long);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Position (cm)');
        ylabel('Accel (cm/s^2)');
        title('Accel x Pos (Long)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end
    
    %% Panel 16 (row 4, col 4): Acceleration × Position contour (short trials)
    if isfield(rate_maps_2d, 'accel_short') && ~isempty(rate_maps_2d.accel_short)
        ax_accel_short = subplot(5, 4, 16, 'Parent', fig);
        hold(ax_accel_short, 'on');
        
        % Get bin edges and data
        pos_edges = rate_maps_2d.pos_bin_edges_short;
        accel_edges = rate_maps_2d.accel_bin_edges_short;
        rate_data = rate_maps_2d.accel_short;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
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
        colormap(ax_accel_short, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_accel_short, [cmin, cmax]);
            colorbar(ax_accel_short);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Position (cm)');
        ylabel('Accel (cm/s^2)');
        title('Accel x Pos (Short)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end
    
    %% Panel 17 (row 5, col 1): Velocity × Relative TTG contour (long trials)
    if isfield(rate_maps_2d, 'ttg_vel_long') && ~isempty(rate_maps_2d.ttg_vel_long)
        ax_ttg_vel_long = subplot(5, 4, 17, 'Parent', fig);
        hold(ax_ttg_vel_long, 'on');
        
        % Get bin edges and data
        ttg_edges = rate_maps_2d.ttg_bin_edges_long;
        vel_edges = rate_maps_2d.vel_bin_edges_long;
        rate_data = rate_maps_2d.ttg_vel_long;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
        cmap = jet(256);
        
        % Create patches for each bin with variable height
        for t = 1:length(ttg_edges)-1
            for v = 1:length(vel_edges)-1
                rate_val = rate_data(t, v);
                
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
                if isfield(rate_maps_2d, 'ttg_vel_trial_counts_long') && isfield(rate_maps_2d, 'ttg_vel_n_trials_long')
                    alpha_val = rate_maps_2d.ttg_vel_trial_counts_long(t, v) / rate_maps_2d.ttg_vel_n_trials_long;
                else
                    alpha_val = 1;
                end
                
                % Create rectangle patch
                x = [ttg_edges(t), ttg_edges(t+1), ttg_edges(t+1), ttg_edges(t)];
                y = [vel_edges(v), vel_edges(v), vel_edges(v+1), vel_edges(v+1)];
                patch(ax_ttg_vel_long, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
            end
        end
        
        hold(ax_ttg_vel_long, 'off');
        colormap(ax_ttg_vel_long, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_ttg_vel_long, [cmin, cmax]);
            colorbar(ax_ttg_vel_long);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Time to Goal (%)');
        ylabel('Velocity (cm/s)');
        set(gca, 'XDir', 'reverse');
        title('Vel x TTG (Long)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end
    
    %% Panel 18 (row 5, col 2): Velocity × Relative TTG contour (short trials)
    if isfield(rate_maps_2d, 'ttg_vel_short') && ~isempty(rate_maps_2d.ttg_vel_short)
        ax_ttg_vel_short = subplot(5, 4, 18, 'Parent', fig);
        hold(ax_ttg_vel_short, 'on');
        
        % Get bin edges and data
        ttg_edges = rate_maps_2d.ttg_bin_edges_short;
        vel_edges = rate_maps_2d.vel_bin_edges_short;
        rate_data = rate_maps_2d.ttg_vel_short;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
        cmap = jet(256);
        
        % Create patches for each bin with variable height
        for t = 1:length(ttg_edges)-1
            for v = 1:length(vel_edges)-1
                rate_val = rate_data(t, v);
                
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
                if isfield(rate_maps_2d, 'ttg_vel_trial_counts_short') && isfield(rate_maps_2d, 'ttg_vel_n_trials_short')
                    alpha_val = rate_maps_2d.ttg_vel_trial_counts_short(t, v) / rate_maps_2d.ttg_vel_n_trials_short;
                else
                    alpha_val = 1;
                end
                
                % Create rectangle patch
                x = [ttg_edges(t), ttg_edges(t+1), ttg_edges(t+1), ttg_edges(t)];
                y = [vel_edges(v), vel_edges(v), vel_edges(v+1), vel_edges(v+1)];
                patch(ax_ttg_vel_short, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
            end
        end
        
        hold(ax_ttg_vel_short, 'off');
        colormap(ax_ttg_vel_short, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_ttg_vel_short, [cmin, cmax]);
            colorbar(ax_ttg_vel_short);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Time to Goal (%)');
        ylabel('Velocity (cm/s)');
        set(gca, 'XDir', 'reverse');
        title('Vel x TTG (Short)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end
    
    %% Panel 19 (row 5, col 3): Acceleration × Relative TTG contour (long trials)
    if isfield(rate_maps_2d, 'ttg_accel_long') && ~isempty(rate_maps_2d.ttg_accel_long)
        ax_ttg_accel_long = subplot(5, 4, 19, 'Parent', fig);
        hold(ax_ttg_accel_long, 'on');
        
        % Get bin edges and data
        ttg_edges = rate_maps_2d.ttg_bin_edges_long;
        accel_edges = rate_maps_2d.accel_bin_edges_long;
        rate_data = rate_maps_2d.ttg_accel_long;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
        cmap = jet(256);
        
        % Create patches for each bin
        for t = 1:length(ttg_edges)-1
            for a = 1:length(accel_edges)-1
                rate_val = rate_data(t, a);
                
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
                if isfield(rate_maps_2d, 'ttg_accel_trial_counts_long') && isfield(rate_maps_2d, 'ttg_accel_n_trials_long')
                    alpha_val = rate_maps_2d.ttg_accel_trial_counts_long(t, a) / rate_maps_2d.ttg_accel_n_trials_long;
                else
                    alpha_val = 1;
                end
                
                % Create rectangle patch
                x = [ttg_edges(t), ttg_edges(t+1), ttg_edges(t+1), ttg_edges(t)];
                y = [accel_edges(a), accel_edges(a), accel_edges(a+1), accel_edges(a+1)];
                patch(ax_ttg_accel_long, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
            end
        end
        
        hold(ax_ttg_accel_long, 'off');
        colormap(ax_ttg_accel_long, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_ttg_accel_long, [cmin, cmax]);
            colorbar(ax_ttg_accel_long);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Time to Goal (%)');
        ylabel('Accel (cm/s^2)');
        set(gca, 'XDir', 'reverse');
        title('Accel x TTG (Long)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end
    
    %% Panel 20 (row 5, col 4): Acceleration × Relative TTG contour (short trials)
    if isfield(rate_maps_2d, 'ttg_accel_short') && ~isempty(rate_maps_2d.ttg_accel_short)
        ax_ttg_accel_short = subplot(5, 4, 20, 'Parent', fig);
        hold(ax_ttg_accel_short, 'on');
        
        % Get bin edges and data
        ttg_edges = rate_maps_2d.ttg_bin_edges_short;
        accel_edges = rate_maps_2d.accel_bin_edges_short;
        rate_data = rate_maps_2d.ttg_accel_short;
        
        % Get colormap limits
        cmin = min(rate_data(:), [], 'omitnan');
        cmax = max(rate_data(:), [], 'omitnan');
        cmap = jet(256);
        
        % Create patches for each bin
        for t = 1:length(ttg_edges)-1
            for a = 1:length(accel_edges)-1
                rate_val = rate_data(t, a);
                
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
                if isfield(rate_maps_2d, 'ttg_accel_trial_counts_short') && isfield(rate_maps_2d, 'ttg_accel_n_trials_short')
                    alpha_val = rate_maps_2d.ttg_accel_trial_counts_short(t, a) / rate_maps_2d.ttg_accel_n_trials_short;
                else
                    alpha_val = 1;
                end
                
                % Create rectangle patch
                x = [ttg_edges(t), ttg_edges(t+1), ttg_edges(t+1), ttg_edges(t)];
                y = [accel_edges(a), accel_edges(a), accel_edges(a+1), accel_edges(a+1)];
                patch(ax_ttg_accel_short, x, y, face_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_val);
            end
        end
        
        hold(ax_ttg_accel_short, 'off');
        colormap(ax_ttg_accel_short, 'parula');
        if ~isnan(cmin) && ~isnan(cmax) && cmax > cmin
            caxis(ax_ttg_accel_short, [cmin, cmax]);
            colorbar(ax_ttg_accel_short);
        end
        set(gca, 'Color', [1 1 1]);
        
        xlabel('Time to Goal (%)');
        ylabel('Accel (cm/s^2)');
        set(gca, 'XDir', 'reverse');
        title('Accel x TTG (Short)', 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    end

    
    % Overall figure title (replace underscores with spaces to prevent subscript rendering)
    probe_id_display = strrep(probe_id, '_', ' ');
    
    % Adjust subplot positions for better use of space with tighter margins
    % and proper spacing between rows/columns
    all_axes = findall(fig, 'Type', 'axes');
    
    % Define new layout parameters
    left_margin = 0.05;      % Reduced from default
    right_margin = 0.02;     % Reduced from default
    bottom_margin = 0.06;    % Increased space at bottom for statistics text
    top_margin = 0.04;       % Reduced space at top
    
    % Spacing between subplots
    h_spacing = 0.06;        % Horizontal spacing between columns
    v_spacing = 0.035;       % Vertical spacing between rows (white gap to avoid label overlap)
    
    % Calculate available space
    total_width = 1 - left_margin - right_margin;
    total_height = 1 - top_margin - bottom_margin;
    
    % Calculate subplot dimensions
    n_cols = 4;
    n_rows = 5;
    subplot_width = (total_width - (n_cols - 1) * h_spacing) / n_cols;
    subplot_height = (total_height - (n_rows - 1) * v_spacing) / n_rows;
    
    % Create a mapping of original positions to new positions
    % We need to identify which subplot each axis is based on its original position
    for ax_idx = 1:length(all_axes)
        ax = all_axes(ax_idx);
        pos = get(ax, 'Position');
        
        % Skip very small axes (likely the stats text axis)
        if pos(4) < 0.05 || pos(3) < 0.05
            continue;
        end
        
        % Identify subplot row and column from ORIGINAL position
        % MATLAB's default subplot positions (approximately):
        % Columns: start at ~0.13, spaced by ~0.2125
        % Rows: start at ~0.775 for row 1, decrease by ~0.157 per row
        
        % More robust: find closest match to expected subplot positions
        % Expected left positions for columns 1-4
        expected_lefts = [0.13, 0.3425, 0.555, 0.7675];
        % Expected bottom positions for rows 1-5
        expected_bottoms = [0.7814, 0.6244, 0.4674, 0.3104, 0.1534];
        
        % Find which column (1-4)
        [~, col] = min(abs(pos(1) - expected_lefts));
        
        % Find which row (1-5)
        [~, row] = min(abs(pos(2) - expected_bottoms));
        
        % Calculate new position
        new_left = left_margin + (col - 1) * (subplot_width + h_spacing);
        new_bottom = bottom_margin + (n_rows - row) * (subplot_height + v_spacing);
        
        % Set new position
        set(ax, 'Position', [new_left, new_bottom, subplot_width, subplot_height]);
        
        % Reduce font sizes
        set(ax, 'FontSize', 7);
    end
    
    % Create title using annotation positioned at top with minimal margin
    annotation('textbox', [0, 0.98, 1, 0.02], 'String', sprintf('Cluster %d - %s', cluster_id, probe_id_display), ...
               'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center', ...
               'EdgeColor', 'none', 'VerticalAlignment', 'middle');
    
    %% Bottom area: Statistical test results (if available)
    if ~isempty(stats)
        
        % Create invisible axis at the bottom for text annotation
        ax_text = axes('Parent', fig, 'Position', [0.05, 0.005, 0.9, 0.05], 'Visible', 'off');
        
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
        text(ax_text, 0.01, 0.5, text_content, 'FontName', 'FixedWidth', 'FontSize', 4, ...
             'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
             'Interpreter', 'tex');
    end
end
