function fig = plot_average_velocity_occupancy_profiles(trial_groups, group_names, group_labels, ...
                                                        analyzer, clusters, probe_id, ctl)
% PLOT_AVERAGE_VELOCITY_OCCUPANCY_PROFILES Create 2x2 figure with velocity and occupancy plots
%
% SYNTAX:
%   fig = plot_average_velocity_occupancy_profiles(trial_groups, group_names, group_labels, ...
%                                                   analyzer, clusters, probe_id, ctl)
%
% INPUTS:
%   trial_groups              - Struct with 'long' and 'short' trial arrays
%   group_names               - Cell array of group names {'long', 'short'}
%   group_labels              - Cell array of group labels for display
%   analyzer                  - SpatialTuningAnalyzer instance
%   clusters                  - Array of cluster objects
%   probe_id                  - String identifier for probe
%   ctl                       - RC2Analysis controller instance
%
% OUTPUTS:
%   fig                       - Figure handle
%
% DESCRIPTION:
%   Creates a 2x2 figure with:
%   - Average velocity profile with individual trials
%   - Velocity heatmap across all trials
%   - Average occupancy profile with individual trials
%   - Occupancy heatmap across all trials

    fig = ctl.figs.a4figure('landscape');
    
    % Reconstruct per-trial velocity and occupancy data
    trial_vel_by_group = struct();
    trial_occ_by_group = struct();
    
    % Collect data for heatmaps (preserving chronological order)
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
        % Fallback if global indices are missing
        max_global_idx = 0;
        for g = 1:2
            max_global_idx = max_global_idx + length(trial_groups.(group_names{g}));
        end
    else
        max_global_idx = max(all_idx_vec);
    end
    
    trial_cells = cell(max_global_idx, 1);
    
    % Process trials for each group
    for g = 1:2
        group = group_names{g};
        trials_struct = trial_groups.(group);
        edges = analyzer.bin_config.(group).edges;
        bin_centers = analyzer.bin_config.(group).centers;
        n_bins = analyzer.bin_config.(group).n_bins;
        
        n_trials = length(trials_struct);
        
        % Get per-trial occupancy and velocity directly from analyzer
        % (stored for the first cluster)
        vel_per_trial = analyzer.firing_rates.(group).vel_per_trial{1};
        occ_per_trial = analyzer.firing_rates.(group).occ_per_trial{1};
        
        for k = 1:n_trials
            
            % Store for heatmap
            trial_info = struct();
            trial_info.group = group;
            if isfield(trials_struct(k), 'global_trial_idx')
                trial_info.global_idx = double(trials_struct(k).global_trial_idx);
            else
                trial_info.global_idx = k;
            end
            
            if strcmp(group, 'long')
                trial_info.color = [0 0 0.8];  % blue
                trial_info.bin_centers = bin_centers;
            else
                trial_info.color = [0.8 0 0];  % red
                trial_info.bin_centers = bin_centers + 60;
            end
            trial_info.vel = vel_per_trial(k, :);
            trial_info.occ = occ_per_trial(k, :);
            
            if trial_info.global_idx > 0 && trial_info.global_idx <= max_global_idx
                trial_cells{trial_info.global_idx} = trial_info;
            end
        end
        
        trial_vel_by_group.(group) = vel_per_trial;
        trial_occ_by_group.(group) = occ_per_trial;
    end
    
    % Collapse non-empty cells
    trial_cells = trial_cells(~cellfun(@isempty, trial_cells));
    all_trials_data = vertcat(trial_cells{:});
    n_total_trials = length(all_trials_data);
    
    % Common bins for heatmaps (0 to 120 cm)
    pos_bin_centers_common = 0:2:120;
    
    % --- Subplot 1: Average velocity with individual trials ---
    subplot(2, 2, 1, 'Parent', fig);
    hold on;
    for g = 1:2
        group = group_names{g};
        bin_centers = analyzer.bin_config.(group).plot_centers;
        vel_trials = trial_vel_by_group.(group);
        
        % Plot individual trials with thin, semi-transparent lines
        for k = 1:size(vel_trials, 1)
            if g == 1  % long trials - blue
                plot(bin_centers, vel_trials(k, :), 'Color', [0.5 0.5 1], ...
                     'LineWidth', 0.5, 'HandleVisibility', 'off');
            else  % short trials - red
                plot(bin_centers, vel_trials(k, :), 'Color', [1 0.5 0.5], ...
                     'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
        end
        
        % Plot average with thick line (use first cluster's avg_velocity)
        avg_vel = analyzer.firing_rates.(group).avg_velocity(1, :);
        if g == 1  % long trials - blue
            plot(bin_centers, avg_vel, 'Color', [0 0 0.8], 'LineWidth', 2, ...
                 'DisplayName', group_labels{g});
        else  % short trials - red
            plot(bin_centers, avg_vel, 'Color', [0.8 0 0], 'LineWidth', 2, ...
                 'DisplayName', group_labels{g});
        end
    end
    hold off;
    xlabel('Position (cm)');
    ylabel('Avg Velocity (cm/s)');
    grid on;
    
    % --- Subplot 2: Velocity Heatmap ---
    subplot(2, 2, 2, 'Parent', fig);
    vel_heatmap = create_trial_heatmap(all_trials_data, pos_bin_centers_common, 'vel');
    image(pos_bin_centers_common, 1:n_total_trials, vel_heatmap);
    set(gca, 'YDir', 'reverse');
    ylabel('Trial #');
    xlabel('Position (cm)');
    xlim([min(pos_bin_centers_common), max(pos_bin_centers_common)]);
    title('Velocity Heatmap');
    
    % --- Subplot 3: Average occupancy with individual trials ---
    subplot(2, 2, 3, 'Parent', fig);
    hold on;
    for g = 1:2
        group = group_names{g};
        bin_centers = analyzer.bin_config.(group).plot_centers;
        occ_trials = trial_occ_by_group.(group);
        
        % Plot individual trials with thin, semi-transparent lines
        for k = 1:size(occ_trials, 1)
            if g == 1  % long trials - blue
                plot(bin_centers, occ_trials(k, :), 'Color', [0.5 0.5 1], ...
                     'LineWidth', 0.5, 'HandleVisibility', 'off');
            else  % short trials - red
                plot(bin_centers, occ_trials(k, :), 'Color', [1 0.5 0.5], ...
                     'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
        end
        
        % Plot average with thick line (use first cluster's occupancy)
        avg_occ = analyzer.firing_rates.(group).occupancy(1, :);
        if g == 1  % long trials - blue
            plot(bin_centers, avg_occ, 'Color', [0 0 0.8], 'LineWidth', 2, ...
                 'DisplayName', group_labels{g});
        else  % short trials - red
            plot(bin_centers, avg_occ, 'Color', [0.8 0 0], 'LineWidth', 2, ...
                 'DisplayName', group_labels{g});
        end
    end
    hold off;
    xlabel('Position (cm)');
    ylabel('Occupancy (s)');
    grid on;
    
    % --- Subplot 4: Occupancy Heatmap ---
    subplot(2, 2, 4, 'Parent', fig);
    occ_heatmap = create_trial_heatmap(all_trials_data, pos_bin_centers_common, 'occ');
    image(pos_bin_centers_common, 1:n_total_trials, occ_heatmap);
    set(gca, 'YDir', 'reverse');
    ylabel('Trial #');
    xlabel('Position (cm)');
    xlim([min(pos_bin_centers_common), max(pos_bin_centers_common)]);
    title('Occupancy Heatmap');
    
    % Add figure title
    FigureTitle(fig, sprintf('Average Profiles - %s', probe_id));
end


function rgb_image = create_trial_heatmap(all_trials_data, pos_bin_centers_common, field_name)
% CREATE_TRIAL_HEATMAP Generate RGB heatmap image from trial data
%
% INPUTS:
%   all_trials_data         - Struct array with trial info
%   pos_bin_centers_common  - Common position bins for interpolation
%   field_name              - Field to plot ('vel' or 'occ')
%
% OUTPUTS:
%   rgb_image               - RGB image array for display

    n_trials = length(all_trials_data);
    heatmap_data = nan(n_trials, length(pos_bin_centers_common));
    
    % Interpolate data to common bins
    for t = 1:n_trials
        centers = all_trials_data(t).bin_centers;
        values = all_trials_data(t).(field_name);
        heatmap_data(t, :) = interp1(centers, values, pos_bin_centers_common, 'linear', nan);
    end
    
    % Normalize data
    data_min = nanmin(heatmap_data(:));
    data_max = nanmax(heatmap_data(:));
    if isnan(data_min) || isnan(data_max) || data_max == data_min
        data_min = 0;
        data_max = 1;
    end
    data_norm = (heatmap_data - data_min) / (data_max - data_min);
    data_norm(isnan(data_norm)) = 0;
    
    % Create RGB image
    rgb_image = ones(size(data_norm, 1), size(data_norm, 2), 3);
    for t = 1:n_trials
        base_color = all_trials_data(t).color;
        for c = 1:3
            rgb_image(t, :, c) = 1 - data_norm(t, :) * (1 - base_color(c));
        end
    end
end
