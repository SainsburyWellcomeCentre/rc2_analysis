function fig = plot_time_to_goal_analysis(trial_groups, group_names, group_labels, all_bin_centers_by_group, all_edges_by_group, probe_id, ctl)
% PLOT_TIME_TO_GOAL_ANALYSIS Create time-to-goal analysis figure (4x2 subplots)
%
%   fig = plot_time_to_goal_analysis(trial_groups, group_names, group_labels, ...
%         all_bin_centers_by_group, all_edges_by_group, probe_id, ctl)
%
% Inputs:
%   trial_groups - Struct with 'long' and 'short' fields containing trial arrays
%   group_names - Cell array of group names {'long', 'short'}
%   group_labels - Cell array of display labels for groups
%   all_bin_centers_by_group - Struct with bin centers for each group
%   all_edges_by_group - Struct with bin edges for each group
%   probe_id - String identifier for the probe
%   ctl - Control structure with figure management (ctl.figs.a4figure)
%
% Outputs:
%   fig - Figure handle
%
% The figure contains 4 rows × 2 columns:
%   Row 1: Absolute TTG vs Position (line plot & heatmap)
%   Row 2: Velocity vs Relative TTG (%) (line plot & heatmap)
%   Row 3: Position vs Relative TTG (%) (line plot & heatmap)
%   Row 4: Occupancy vs Relative TTG (%) (line plot & heatmap)

% Create figure
fig = ctl.figs.a4figure('landscape');

% Collect trial indices for chronological ordering
all_idx_vec = [];
for g = 1:2
    trials_struct = trial_groups.(group_names{g});
    for k = 1:length(trials_struct)
        if isfield(trials_struct(k), 'global_trial_idx')
            all_idx_vec(end+1, 1) = double(trials_struct(k).global_trial_idx); %#ok<AGROW>
        end
    end
end
all_idx_vec = all_idx_vec(isfinite(all_idx_vec));
if isempty(all_idx_vec)
    error('No valid global trial indices found for TTG plotting.');
end
max_global_idx = max(all_idx_vec);
trial_cells = cell(max_global_idx, 1);

% Process each group
for g = 1:2
    group = group_names{g};
    trials_struct = trial_groups.(group);
    n_trials = length(trials_struct);
    
    % Storage for plot 1: absolute time-to-goal vs position
    pos_bins{g} = all_bin_centers_by_group.(group);
    edges = all_edges_by_group.(group);
    n_bins = length(pos_bins{g});
    ttg_vs_pos_trials{g} = nan(n_trials, n_bins);
    
    % Storage for plots 2 & 3: velocity and position vs relative TTG
    ttg_bin_edges = 0:2:100;  % 2% bins
    ttg_bin_centers = ttg_bin_edges(1:end-1) + 1;
    n_ttg_bins = length(ttg_bin_centers);
    vel_vs_ttg_trials{g} = nan(n_trials, n_ttg_bins);
    pos_vs_ttg_trials{g} = nan(n_trials, n_ttg_bins);
    occ_vs_ttg_trials{g} = nan(n_trials, n_ttg_bins);
    
    for k = 1:n_trials
        trial = trials_struct(k).trial;
        motion_mask = trial.motion_mask();
        pos = trial.position(motion_mask);
        vel = trial.velocity(motion_mask);
        time = trial.probe_t(motion_mask);
        
        % Store trial info for heatmaps
        trial_info.group = group;
        if isfield(trials_struct(k), 'global_trial_idx')
            trial_info.global_idx = double(trials_struct(k).global_trial_idx);
        else
            trial_info.global_idx = k;
        end
        if strcmp(group, 'long')
            trial_info.color = [0 0 0.8];  % blue
        else
            trial_info.color = [0.8 0 0];  % red
        end
        trial_info.pos_bin_centers = pos_bins{g};
        
        if ~isempty(pos) && ~isempty(time) && length(time) > 1
            % Plot 1: Compute absolute time-to-goal (seconds) vs position
            time_to_goal = time(end) - time;
            ttg_binned = nan(1, n_bins);
            for b = 1:n_bins
                in_bin = pos >= edges(b) & pos < edges(b+1);
                if sum(in_bin) > 0
                    ttg_binned(b) = nanmean(time_to_goal(in_bin));
                end
            end
            ttg_vs_pos_trials{g}(k, :) = ttg_binned;
            trial_info.ttg_vs_pos = ttg_binned;
            
            % Plots 2 & 3: Normalize TTG to percentage based on MOTION TIME
            % Bin by time (not samples) to ensure constant occupancy per bin
            
            % Create continuous time vector by removing gaps (stationary periods)
            % Count actual motion duration from number of samples and sampling rate
            n_samples = length(time);
            sampling_rate = 10000;  % Hz
            trial_duration = n_samples / sampling_rate;  % Total motion duration in seconds
            
            % Create continuous time from 0 to trial_duration
            continuous_time = linspace(0, trial_duration, n_samples)';
            
            if trial_duration > 0
                % Compute time-to-goal for each motion sample (time from sample to end)
                time_to_goal = trial_duration - continuous_time;
                % Normalize to percentage of motion duration
                ttg_percent = (time_to_goal / trial_duration) * 100;
                
                vel_binned = nan(1, n_ttg_bins);
                pos_binned = nan(1, n_ttg_bins);
                occ_binned = nan(1, n_ttg_bins);
                
                % Constant occupancy per bin based on motion duration
                bin_width_percent = ttg_bin_edges(2) - ttg_bin_edges(1);  % 2%
                occ_per_bin = trial_duration * (bin_width_percent / 100);  % Constant occupancy in seconds
                
                % Debug: Check motion mask and data availability
                full_duration = trial.probe_t(end) - trial.probe_t(1);
                sampling_rate = length(trial.probe_t) / full_duration;
                expected_samples = full_duration * sampling_rate;
                coverage_pct = (length(time) / expected_samples) * 100;
                
                fprintf('      Trial %d (Group: %s): motion_mask has %d samples, motion duration=%.3fs, full duration=%.3fs (%.1f%% coverage)\n', ...
                    k, group, length(time), trial_duration, full_duration, coverage_pct);
                fprintf('        TTG time-based range: [%.1f%%, %.1f%%], constant occupancy per bin: %.4fs\n', ...
                    min(ttg_percent), max(ttg_percent), occ_per_bin);
                
                for b = 1:n_ttg_bins
                    % Find samples in this TTG time bin
                    in_bin = ttg_percent >= ttg_bin_edges(b) & ttg_percent < ttg_bin_edges(b+1);
                    n_samples_in_bin = sum(in_bin);
                    
                    % Occupancy is ALWAYS constant (based on bin width in time)
                    occ_binned(b) = occ_per_bin;
                    
                    if n_samples_in_bin > 0
                        % Map velocity and position from motion samples in this time bin
                        vel_binned(b) = nanmean(vel(in_bin));
                        pos_binned(b) = nanmean(pos(in_bin));
                    else
                        % Debug: Report bins with no motion samples
                        fprintf('        WARNING: TTG bin %d (%.1f-%.1f%%) has NO motion samples!\n', ...
                            b, ttg_bin_edges(b), ttg_bin_edges(b+1));
                    end
                end
                
                % Debug: Report how many bins are missing
                n_missing_vel = sum(isnan(vel_binned));
                n_missing_pos = sum(isnan(pos_binned));
                n_missing_occ = sum(isnan(occ_binned));
                if n_missing_vel > 0 || n_missing_pos > 0 || n_missing_occ > 0
                    fprintf('        Missing bins - Vel: %d, Pos: %d, Occ: %d (out of %d)\n', ...
                        n_missing_vel, n_missing_pos, n_missing_occ, n_ttg_bins);
                end
                
                vel_vs_ttg_trials{g}(k, :) = vel_binned;
                occ_vs_ttg_trials{g}(k, :) = occ_binned;
                % Adjust position for short trials (add 60cm offset for display)
                if strcmp(group, 'short')
                    pos_vs_ttg_trials{g}(k, :) = pos_binned + 60;
                else
                    pos_vs_ttg_trials{g}(k, :) = pos_binned;
                end
                
                trial_info.vel_vs_ttg = vel_binned;
                trial_info.pos_vs_ttg = pos_binned + (strcmp(group, 'short') * 60);
                trial_info.occ_vs_ttg = occ_binned;
                trial_info.has_ttg_data = true;
            else
                trial_info.has_ttg_data = false;
            end
        else
            trial_info.has_ttg_data = false;
        end
        
        trial_cells{trial_info.global_idx} = trial_info;
    end
end

% Collapse non-empty cells preserving original experiment order
trial_cells = trial_cells(~cellfun(@isempty, trial_cells));
all_trials_data = vertcat(trial_cells{:});

% Common bins for heatmaps
pos_bin_centers_common = 0:2:120;  % 2cm bins
ttg_bin_centers = 1:2:99;  % 2% bins for TTG
n_total_trials = length(all_trials_data);
n_ttg_bins = length(ttg_bin_centers);

% === Row 1, Col 1: Absolute TTG vs Position (line plot) ===
subplot(4, 2, 1, 'Parent', fig);
hold on;
for g = 1:2
    ttg_pos = ttg_vs_pos_trials{g};
    pos_axis = pos_bins{g};
    % Plot individual trials
    for k = 1:size(ttg_pos, 1)
        if g == 1
            plot(pos_axis, ttg_pos(k, :), 'Color', [0.5 0.5 1], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        else
            plot(pos_axis, ttg_pos(k, :), 'Color', [1 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
    % Plot mean with thick line
    mean_ttg_pos = nanmean(ttg_pos, 1);
    if g == 1
        plot(pos_axis, mean_ttg_pos, 'Color', [0 0 0.8], 'LineWidth', 2, 'DisplayName', group_labels{g});
    else
        plot(pos_axis, mean_ttg_pos, 'Color', [0.8 0 0], 'LineWidth', 2, 'DisplayName', group_labels{g});
    end
end
hold off;
xlabel('Position (cm)');
ylabel('Time to Goal (s)');
grid on;

% === Row 1, Col 2: Absolute TTG vs Position (heatmap) ===
subplot(4, 2, 2, 'Parent', fig);
ttg_heatmap = nan(n_total_trials, length(pos_bin_centers_common));
for t = 1:n_total_trials
    if all_trials_data(t).has_ttg_data && isfield(all_trials_data(t), 'ttg_vs_pos')
        centers = all_trials_data(t).pos_bin_centers;
        values = all_trials_data(t).ttg_vs_pos;
        ttg_heatmap(t, :) = interp1(centers, values, pos_bin_centers_common, 'linear', nan);
    end
end
ttg_heatmap_rgb = create_colored_heatmap(ttg_heatmap, all_trials_data);
image(pos_bin_centers_common, 1:n_total_trials, ttg_heatmap_rgb);
set(gca, 'YDir', 'reverse');
ylabel('Trial # (chronological)');
xlabel('Position (cm)');
xlim([min(pos_bin_centers_common), max(pos_bin_centers_common)]);
title('Time to Goal (s) Heatmap');

% === Row 2: Velocity vs Relative TTG ===
plot_ttg_variable(fig, 3, 4, vel_vs_ttg_trials, ttg_bin_centers, all_trials_data, ...
    group_labels, 'Velocity (cm/s)', 'Velocity vs Relative Time-to-Goal', 'Velocity Heatmap');

% === Row 3: Position vs Relative TTG ===
plot_ttg_variable(fig, 5, 6, pos_vs_ttg_trials, ttg_bin_centers, all_trials_data, ...
    group_labels, 'Position (cm)', 'Position vs Relative Time-to-Goal', 'Position Heatmap');

% === Row 4: Occupancy vs Relative TTG ===
plot_ttg_variable(fig, 7, 8, occ_vs_ttg_trials, ttg_bin_centers, all_trials_data, ...
    group_labels, 'Occupancy (s)', 'Occupancy vs Relative Time-to-Goal', 'Occupancy Heatmap');

% Add figure title
FigureTitle(fig, sprintf('Time-to-Goal Analysis - %s', probe_id));

end

% Helper function to plot a TTG variable (line plot + heatmap)
function plot_ttg_variable(fig, line_idx, heatmap_idx, data_trials, ttg_bin_centers, ...
    all_trials_data, group_labels, ylabel_str, line_title, heatmap_title)

n_total_trials = length(all_trials_data);
n_ttg_bins = length(ttg_bin_centers);

% Line plot
subplot(4, 2, line_idx, 'Parent', fig);
hold on;
for g = 1:2
    data = data_trials{g};
    % Plot individual trials
    for k = 1:size(data, 1)
        if g == 1
            plot(ttg_bin_centers, data(k, :), 'Color', [0.5 0.5 1], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        else
            plot(ttg_bin_centers, data(k, :), 'Color', [1 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
    % Plot mean with thick line
    mean_data = nanmean(data, 1);
    if g == 1
        plot(ttg_bin_centers, mean_data, 'Color', [0 0 0.8], 'LineWidth', 2, 'DisplayName', group_labels{g});
    else
        plot(ttg_bin_centers, mean_data, 'Color', [0.8 0 0], 'LineWidth', 2, 'DisplayName', group_labels{g});
    end
end
hold off;
xlabel('Time to Goal (%)');
ylabel(ylabel_str);
grid on;
xlim([0 100]);
set(gca, 'XDir', 'reverse');
title(line_title);

% Heatmap
subplot(4, 2, heatmap_idx, 'Parent', fig);
heatmap_data = nan(n_total_trials, n_ttg_bins);
field_name = get_field_name_from_trials(all_trials_data, line_idx);
for t = 1:n_total_trials
    if all_trials_data(t).has_ttg_data && isfield(all_trials_data(t), field_name)
        heatmap_data(t, :) = all_trials_data(t).(field_name);
    end
end
heatmap_rgb = create_colored_heatmap(heatmap_data, all_trials_data);
image(ttg_bin_centers, 1:n_total_trials, heatmap_rgb);
set(gca, 'YDir', 'reverse');
set(gca, 'XDir', 'reverse');
ylabel('Trial # (chronological)');
xlabel('Time to Goal (%)');
xlim([0 100]);
title(heatmap_title);

end

% Helper to get field name from subplot index
function field_name = get_field_name_from_trials(~, line_idx)
switch line_idx
    case 3
        field_name = 'vel_vs_ttg';
    case 5
        field_name = 'pos_vs_ttg';
    case 7
        field_name = 'occ_vs_ttg';
    otherwise
        field_name = '';
end
end

% Helper function to create colored heatmap
function rgb_image = create_colored_heatmap(heatmap_data, all_trials_data)
% Normalize heatmap data
data_min = nanmin(heatmap_data(:));
data_max = nanmax(heatmap_data(:));
if isnan(data_min) || isnan(data_max) || data_max == data_min
    data_min = 0;
    data_max = 1;
end
data_norm = (heatmap_data - data_min) / (data_max - data_min);
data_norm(isnan(data_norm)) = 0;

% Create RGB image with group-specific colors
rgb_image = ones(size(data_norm, 1), size(data_norm, 2), 3);
for t = 1:size(data_norm, 1)
    base_color = all_trials_data(t).color;
    for c = 1:3
        rgb_image(t, :, c) = 1 - data_norm(t, :) * (1 - base_color(c));
    end
end
end
