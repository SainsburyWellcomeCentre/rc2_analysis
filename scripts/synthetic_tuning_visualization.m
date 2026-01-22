% SYNTHETIC_TUNING_VISUALIZATION
% Generate synthetic firing rates based on different tuning properties and
% visualize how they manifest across different representations (absolute space,
% relative space, velocity, acceleration).
%
% This script allows you to specify one primary tuning characteristic
% (e.g., velocity tuning) and see how it appears in all other representations.
% It uses real velocity profiles from experimental data to create realistic
% synthetic neural responses.

%% ==================== PARAMETERS ====================

% --- Data Selection ---
probe_id = 'CAA-1123476a_rec1';  % Which probe's data to use
trial_index_long = 'all';   % Which long trial's velocity profile to use (1-based index or 'all')
trial_index_short = 'all';  % Which short trial's velocity profile to use (1-based index or 'all')

% --- Primary Tuning Type ---
% Choose which tuning characteristic to start from
% Options: 'absolute_space', 'relative_space', 'velocity', 'acceleration'
primary_tuning = 'acceleration';

% --- Absolute Space Tuning Parameters (for absolute positions on track) ---
% Long trials: 0-120 cm, Short trials: 60-120 cm
% If mu=60, cell fires at center of long trials and start of short trials
abs_space = struct();
abs_space.enabled = true;
abs_space.type = 'gaussian';  % 'gaussian', 'linear', 'flat'
abs_space.mu = 60;            % Peak position in absolute coordinates (cm)
abs_space.sigma = 20;         % Width of Gaussian (cm)
abs_space.amplitude = 50;     % Peak firing rate (Hz)
abs_space.baseline = 2;       % Baseline firing rate (Hz)

% --- Relative Space Tuning Parameters (0-100% of track) ---
rel_space = struct();
rel_space.enabled = true;
rel_space.type = 'gaussian';  % 'gaussian', 'linear', 'flat'
rel_space.mu = 10;           % Peak position (fraction 0-1)
rel_space.sigma = 15;        % Width of Gaussian (fraction)
rel_space.amplitude = 50;     % Peak firing rate (Hz)
rel_space.baseline = 2;       % Baseline firing rate (Hz)

% --- Velocity Tuning Parameters ---
vel_tuning = struct();
vel_tuning.enabled = true;
vel_tuning.type = 'gaussian';  % 'linear_positive', 'linear_negative', 'gaussian', 'sigmoid'
vel_tuning.mu = 1;            % Peak velocity (cm/s) for Gaussian
vel_tuning.sigma = 1;         % Width for Gaussian (cm/s)
vel_tuning.amplitude = 100;     % Peak firing rate (Hz)
vel_tuning.baseline = 10;       % Baseline firing rate (Hz)
vel_tuning.slope = 2;        % Slope for linear (Hz per cm/s)
vel_tuning.midpoint = 10;      % Midpoint for sigmoid (cm/s)
vel_tuning.steepness = 0.5;    % Steepness for sigmoid

% --- Acceleration Tuning Parameters ---
accel_tuning = struct();
accel_tuning.enabled = true;
accel_tuning.type = 'u_shaped_positive';  % 'linear_positive', 'linear_negative', 'sigmoid', 
                                           % 'u_shaped_positive', 'u_shaped_negative'
accel_tuning.amplitude = 30;     % Amplitude (Hz)
accel_tuning.baseline = 1;      % Baseline firing rate (Hz)
accel_tuning.slope = 0.5;        % Slope for linear (Hz per cm/s²)
accel_tuning.midpoint = 0;       % Midpoint for sigmoid or U-shape (cm/s²)
accel_tuning.steepness = 0.05;   % Steepness for sigmoid
accel_tuning.curvature = 0.1;   % Curvature for U-shape (Hz per (cm/s²)²)

%% ==================== SETUP ====================

% Initialize controller
ctl = RC2Analysis();

% Get path configuration
cfg = path_config();

% Define cache file path - search in common spatial_firing_rate output directories
cache_filename = sprintf('%s_spatial_analysis_cache.mat', probe_id);
cache_search_dirs = {
    fullfile(cfg.figure_dir, 'spatial_firing_rate', 'ambient_light'),
    fullfile(cfg.figure_dir, 'spatial_firing_rate', 'training_running')
};

% Find cache file
cache_filepath = '';
for d = 1:length(cache_search_dirs)
    test_path = fullfile(cache_search_dirs{d}, cache_filename);
    if exist(test_path, 'file')
        cache_filepath = test_path;
        fprintf('Found cache file in: %s\n', cache_search_dirs{d});
        break;
    end
end

if isempty(cache_filepath)
    error('Cache file not found: %s\nSearched in:\n  %s\n\nPlease run spatial_firing_rate_profile.m first to generate the cache file.', ...
        cache_filename, strjoin(cache_search_dirs, '\n  '));
end

fprintf('Loading cache from: %s\n', cache_filepath);

%% ==================== LOAD DATA ====================

% Load cache - the cache contains direct velocity and position data per trial
cache = load(cache_filepath);

% Get trial velocity and position data from cache
% The cache stores velocity per trial in all_vel_per_trial_by_group
if ~isfield(cache, 'all_vel_per_trial_by_group')
    error('Cache file does not contain velocity data. Please regenerate cache with spatial_firing_rate_profile.m');
end

vel_per_trial = cache.all_vel_per_trial_by_group;
if ~isfield(vel_per_trial, trial_type)
    error('Trial type ''%s'' not found in cache', trial_type);
end

% Get position data (stored as spike positions per trial, but we can reconstruct from bin centers)
% Actually, we need to load the formatted data to get proper trial information
fprintf('Loading formatted data for %s...\n', probe_id);
data = ctl.load_formatted_data(probe_id);
sessions = data.motion_sessions();

% Collect long and short trials separately
long_trials = [];
short_trials = [];
global_trial_idx = 0;
for s = 1:length(sessions)
    sess = sessions{s};
    trials_in_session = sess.trials;  % Use .trials method to get trial objects
    for t = 1:length(trials_in_session)
        % Convert to aligned trial (same as collect_and_group_trials.m)
        trial = trials_in_session{t}.to_aligned;
        global_trial_idx = global_trial_idx + 1;
        
        % Check if trial has motion data
        motion_mask = trial.motion_mask();
        if ~any(motion_mask)
            continue;
        end
        
        pos = trial.position(motion_mask);
        if isempty(pos)
            continue;
        end
        
        % Classify trial type based on maximum position
        max_pos = max(pos);
        if max_pos > cache.position_threshold
            long_trials = [long_trials; trial]; %#ok<AGROW>
        else
            short_trials = [short_trials; trial]; %#ok<AGROW>
        end
    end
end

% Determine if we're analyzing single or multiple trials
use_all_long = ischar(trial_index_long) && strcmpi(trial_index_long, 'all');
use_all_short = ischar(trial_index_short) && strcmpi(trial_index_short, 'all');

% Validate trial indices
if ~use_all_long && trial_index_long > length(long_trials)
    error('Long trial index %d exceeds number of available long trials (%d)', ...
        trial_index_long, length(long_trials));
end
if ~use_all_short && trial_index_short > length(short_trials)
    error('Short trial index %d exceeds number of available short trials (%d)', ...
        trial_index_short, length(short_trials));
end

% Extract data for LONG trial(s)
if use_all_long
    % Extract data from all long trials
    long_trial_data = cell(length(long_trials), 1);
    for i = 1:length(long_trials)
        trial = long_trials(i);
        motion_mask = trial.motion_mask();
        long_trial_data{i}.position = trial.position(motion_mask);
        long_trial_data{i}.velocity = trial.velocity(motion_mask);
        long_trial_data{i}.acceleration = trial.acceleration(motion_mask);
        long_trial_data{i}.time = trial.probe_t(motion_mask);
    end
    fprintf('Loaded %d LONG trials\n', length(long_trials));
    % For compatibility, use first trial as reference
    position_long = long_trial_data{1}.position;
    velocity_long = long_trial_data{1}.velocity;
    acceleration_long = long_trial_data{1}.acceleration;
    time_long = long_trial_data{1}.time;
else
    % Extract data for single long trial
    trial_long = long_trials(trial_index_long);
    motion_mask_long = trial_long.motion_mask();
    position_long = trial_long.position(motion_mask_long);
    velocity_long = trial_long.velocity(motion_mask_long);
    acceleration_long = trial_long.acceleration(motion_mask_long);
    time_long = trial_long.probe_t(motion_mask_long);
    
    fprintf('Loaded LONG trial %d: %d samples, duration %.2f s\n', ...
        trial_index_long, length(time_long), time_long(end) - time_long(1));
    fprintf('  Position: %.1f - %.1f cm, Velocity: %.1f - %.1f cm/s, Accel: %.1f - %.1f cm/s²\n', ...
        min(position_long), max(position_long), min(velocity_long), max(velocity_long), ...
        min(acceleration_long), max(acceleration_long));
end

% Extract data for SHORT trial(s)
short_trial_start_position = 60;  % cm - offset for absolute coordinates
if use_all_short
    % Extract data from all short trials
    short_trial_data = cell(length(short_trials), 1);
    for i = 1:length(short_trials)
        trial = short_trials(i);
        motion_mask = trial.motion_mask();
        short_trial_data{i}.position = trial.position(motion_mask);
        short_trial_data{i}.position_absolute = trial.position(motion_mask) + short_trial_start_position;
        short_trial_data{i}.velocity = trial.velocity(motion_mask);
        short_trial_data{i}.acceleration = trial.acceleration(motion_mask);
        short_trial_data{i}.time = trial.probe_t(motion_mask);
    end
    fprintf('Loaded %d SHORT trials\n', length(short_trials));
    % For compatibility, use first trial as reference
    position_short = short_trial_data{1}.position;
    position_short_absolute = short_trial_data{1}.position_absolute;
    velocity_short = short_trial_data{1}.velocity;
    acceleration_short = short_trial_data{1}.acceleration;
    time_short = short_trial_data{1}.time;
else
    % Extract data for single short trial
    trial_short = short_trials(trial_index_short);
    motion_mask_short = trial_short.motion_mask();
    position_short = trial_short.position(motion_mask_short);
    velocity_short = trial_short.velocity(motion_mask_short);
    acceleration_short = trial_short.acceleration(motion_mask_short);
    time_short = trial_short.probe_t(motion_mask_short);
    position_short_absolute = position_short + short_trial_start_position;
    
    fprintf('Loaded SHORT trial %d: %d samples, duration %.2f s\n', ...
        trial_index_short, length(time_short), time_short(end) - time_short(1));
    fprintf('  Position: %.1f - %.1f cm, Velocity: %.1f - %.1f cm/s, Accel: %.1f - %.1f cm/s²\n', ...
        min(position_short), max(position_short), min(velocity_short), max(velocity_short), ...
        min(acceleration_short), max(acceleration_short));
    fprintf('  Short trial absolute positions: %.1f - %.1f cm (offset by %.1f cm)\n', ...
        min(position_short_absolute), max(position_short_absolute), short_trial_start_position);
end

%% ==================== GENERATE SYNTHETIC FIRING RATES ====================

noise_level = 0.1;  % 10% noise

% Generate firing rates for LONG trial(s)
if use_all_long
    % Generate firing rates for all long trials
    for i = 1:length(long_trial_data)
        pos = long_trial_data{i}.position;
        vel = long_trial_data{i}.velocity;
        accel = long_trial_data{i}.acceleration;
        
        switch primary_tuning
            case 'absolute_space'
                fr = generate_spatial_tuning(pos, abs_space);
            case 'relative_space'
                rel_pos = (pos - min(pos)) / (max(pos) - min(pos));
                fr = generate_spatial_tuning(rel_pos * 100, rel_space);
            case 'velocity'
                fr = generate_velocity_tuning(vel, vel_tuning);
            case 'acceleration'
                fr = generate_acceleration_tuning(accel, accel_tuning);
            otherwise
                error('Unknown primary tuning type: %s', primary_tuning);
        end
        
        % Add noise
        fr = fr .* (1 + noise_level * randn(size(fr)));
        fr(fr < 0) = 0;
        long_trial_data{i}.firing_rate = fr;
    end
    % Reference for compatibility
    firing_rate_long = long_trial_data{1}.firing_rate;
else
    % Generate firing rate for single long trial
    switch primary_tuning
        case 'absolute_space'
            firing_rate_long = generate_spatial_tuning(position_long, abs_space);
        case 'relative_space'
            rel_position_long = (position_long - min(position_long)) / (max(position_long) - min(position_long));
            firing_rate_long = generate_spatial_tuning(rel_position_long * 100, rel_space);
        case 'velocity'
            firing_rate_long = generate_velocity_tuning(velocity_long, vel_tuning);
        case 'acceleration'
            firing_rate_long = generate_acceleration_tuning(acceleration_long, accel_tuning);
        otherwise
            error('Unknown primary tuning type: %s', primary_tuning);
    end
    
    % Add noise
    firing_rate_long = firing_rate_long .* (1 + noise_level * randn(size(firing_rate_long)));
    firing_rate_long(firing_rate_long < 0) = 0;
end

% Generate firing rates for SHORT trial(s)
if use_all_short
    % Generate firing rates for all short trials
    for i = 1:length(short_trial_data)
        pos = short_trial_data{i}.position;
        pos_abs = short_trial_data{i}.position_absolute;
        vel = short_trial_data{i}.velocity;
        accel = short_trial_data{i}.acceleration;
        
        switch primary_tuning
            case 'absolute_space'
                fr = generate_spatial_tuning(pos_abs, abs_space);
            case 'relative_space'
                rel_pos = (pos - min(pos)) / (max(pos) - min(pos));
                fr = generate_spatial_tuning(rel_pos * 100, rel_space);
            case 'velocity'
                fr = generate_velocity_tuning(vel, vel_tuning);
            case 'acceleration'
                fr = generate_acceleration_tuning(accel, accel_tuning);
            otherwise
                error('Unknown primary tuning type: %s', primary_tuning);
        end
        
        % Add noise
        fr = fr .* (1 + noise_level * randn(size(fr)));
        fr(fr < 0) = 0;
        short_trial_data{i}.firing_rate = fr;
    end
    % Reference for compatibility
    firing_rate_short = short_trial_data{1}.firing_rate;
else
    % Generate firing rate for single short trial
    switch primary_tuning
        case 'absolute_space'
            firing_rate_short = generate_spatial_tuning(position_short_absolute, abs_space);
        case 'relative_space'
            rel_position_short = (position_short - min(position_short)) / (max(position_short) - min(position_short));
            firing_rate_short = generate_spatial_tuning(rel_position_short * 100, rel_space);
        case 'velocity'
            firing_rate_short = generate_velocity_tuning(velocity_short, vel_tuning);
        case 'acceleration'
            firing_rate_short = generate_acceleration_tuning(acceleration_short, accel_tuning);
        otherwise
            error('Unknown primary tuning type: %s', primary_tuning);
    end
    
    % Add noise
    firing_rate_short = firing_rate_short .* (1 + noise_level * randn(size(firing_rate_short)));
    firing_rate_short(firing_rate_short < 0) = 0;
end

%% ==================== CREATE FIGURE ====================

fig = figure('Position', [100, 100, 1200, 900], 'Color', 'w');

% Color scheme
color_long = [0, 0.4470, 0.7410];    % Blue
color_short = [0.8500, 0.3250, 0.0980];  % Red

% Panel 1 (Row 1, Col 1): Parameters
ax1 = subplot(3, 6, 1:2);
axis off;
% Adjust position to add spacing
pos1 = get(ax1, 'Position');
pos1(3) = pos1(3) * 0.80;  % Reduce width to 80% to create spacing
set(ax1, 'Position', pos1);

% Build description based on primary tuning
switch primary_tuning
    case 'absolute_space'
        tuning_desc = sprintf(['Neuron tuned to absolute position %.0f cm ', ...
            'with %.0f cm spread (Gaussian). Fires at same track ', ...
            'location regardless of trial type.'], abs_space.mu, abs_space.sigma);
    case 'relative_space'
        tuning_desc = sprintf(['Neuron tuned to %.0f%% of track length ', ...
            'with %.0f%% spread (Gaussian). Fires at same relative ', ...
            'position across trial types.'], rel_space.mu*100, rel_space.sigma*100);
    case 'velocity'
        if strcmp(vel_tuning.type, 'gaussian')
            tuning_desc = sprintf(['Neuron tuned to velocity %.1f cm/s ', ...
                'with %.1f cm/s spread (Gaussian). Fires when mouse ', ...
                'runs at specific speed.'], vel_tuning.mu, vel_tuning.sigma);
        else
            tuning_desc = sprintf('Neuron with %s velocity tuning. ', ...
                strrep(vel_tuning.type, '_', ' '));
        end
    case 'acceleration'
        tuning_desc = sprintf(['Neuron with %s acceleration tuning.\n', ...
            'Fires based on changes in running speed.'], ...
            strrep(accel_tuning.type, '_', ' '));
    otherwise
        tuning_desc = '';
end

param_text = {
    '\bf{Synthetic Tuning Parameters}',
    '',
    sprintf('\\bf{Probe:} %s', strrep(probe_id, '_', '\_')),
    sprintf('\\bf{Long Trial:} %s, \\bf{Short Trial:} %s', ...
        char(string(trial_index_long)), char(string(trial_index_short))),
    '',
    sprintf('\\bf{Primary Tuning:} %s', strrep(primary_tuning, '_', ' ')),
    '',
    tuning_desc
};

text(0.05, 0.95, param_text, 'VerticalAlignment', 'top', 'FontSize', 11, ...
    'Interpreter', 'tex');

% Panel 2 (Row 1, Col 2): Long trial running profile
ax2 = subplot(3, 6, [3 4]);
pos2 = get(ax2, 'Position');
pos2(1) = pos2(1) + 0.03;  % Shift right to create gap
set(ax2, 'Position', pos2);
hold on;
if use_all_long
    % Plot individual trials only
    for i = 1:length(long_trial_data)
        t = long_trial_data{i}.time - long_trial_data{i}.time(1);  % Start from 0
        yyaxis left
        plot(t, long_trial_data{i}.position, '-', 'Color', [0, 0, 0, 0.5], 'LineWidth', 0.8);
        yyaxis right
        plot(t, long_trial_data{i}.velocity, '-', 'Color', [color_long, 0.5], 'LineWidth', 0.8);
    end
else
    % Single trial
    yyaxis left
    plot(time_long - time_long(1), position_long, '-', 'Color', 'k', 'LineWidth', 1.5);
    yyaxis right
    plot(time_long - time_long(1), velocity_long, '-', 'Color', color_long, 'LineWidth', 1);
end
yyaxis left
ylabel('Position (cm)', 'Color', 'k');
ylim([0, 120]);
ax = gca;
ax.YAxis(1).Color = 'k';
yyaxis right
ylabel('Velocity (cm/s)', 'Color', color_long);
ax.YAxis(2).Color = color_long;
xlabel('Time (s)');
title('Long Trial Running Profile');
grid on;
hold off;

% Panel 3 (Row 1, Col 3): Short trial running profile
ax3 = subplot(3, 6, [5 6]);
pos3 = get(ax3, 'Position');
pos3(1) = pos3(1) + 0.08;  % Shift right to create gap
set(ax3, 'Position', pos3);
hold on;
if use_all_short
    % Plot individual trials only
    for i = 1:length(short_trial_data)
        t = short_trial_data{i}.time - short_trial_data{i}.time(1);  % Start from 0
        yyaxis left
        plot(t, short_trial_data{i}.position_absolute, '-', 'Color', [0, 0, 0, 0.5], 'LineWidth', 0.8);
        yyaxis right
        plot(t, short_trial_data{i}.velocity, '-', 'Color', [color_short, 0.5], 'LineWidth', 0.8);
    end
else
    % Single trial
    yyaxis left
    plot(time_short - time_short(1), position_short_absolute, '-', 'Color', 'k', 'LineWidth', 1.5);
    yyaxis right
    plot(time_short - time_short(1), velocity_short, '-', 'Color', color_short, 'LineWidth', 1);
end
yyaxis left
ylabel('Position (cm)', 'Color', 'k');
ylim([0, 120]);
ax = gca;
ax.YAxis(1).Color = 'k';
yyaxis right
ylabel('Velocity (cm/s)', 'Color', color_short);
ax.YAxis(2).Color = color_short;
xlabel('Time (s)');
title('Short Trial Running Profile');
grid on;
hold off;

% Panel 4 (Row 2, Col 1): Absolute space tuning
subplot(3, 6, 7:9);
hold on;

% Use bin configuration from cache to match spatial_firing_rate_profile.m
pos_bin_size = cache.bin_size_cm;  % Use cached bin size (typically 2 cm)
pos_bins = 0:pos_bin_size:120;
pos_centers = pos_bins(1:end-1) + pos_bin_size/2;

% Process LONG trials
if use_all_long
    % Plot individual trials (thin) and average (thick)
    all_binned_long = zeros(length(long_trial_data), length(pos_centers));
    for trial_idx = 1:length(long_trial_data)
        binned_rate = zeros(size(pos_centers));
        binned_count = zeros(size(pos_centers));
        pos = long_trial_data{trial_idx}.position;
        fr = long_trial_data{trial_idx}.firing_rate;
        for i = 1:length(pos)
            bin_idx = find(pos(i) >= pos_bins(1:end-1) & pos(i) < pos_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_binned_long(trial_idx, :) = binned_rate;
        plot(pos_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_long, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_binned_long = mean(all_binned_long, 1);
    valid_avg_long = avg_binned_long > 0;
    plot(pos_centers(valid_avg_long), avg_binned_long(valid_avg_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2.5, 'DisplayName', 'Long (avg)');
else
    % Single trial
    binned_rate_long = zeros(size(pos_centers));
    binned_count_long = zeros(size(pos_centers));
    for i = 1:length(position_long)
        bin_idx = find(position_long(i) >= pos_bins(1:end-1) & position_long(i) < pos_bins(2:end), 1);
        if ~isempty(bin_idx)
            binned_rate_long(bin_idx) = binned_rate_long(bin_idx) + firing_rate_long(i);
            binned_count_long(bin_idx) = binned_count_long(bin_idx) + 1;
        end
    end
    valid_bins_long = binned_count_long > 0;
    binned_rate_long(valid_bins_long) = binned_rate_long(valid_bins_long) ./ binned_count_long(valid_bins_long);
    plot(pos_centers(valid_bins_long), binned_rate_long(valid_bins_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
end

% Process SHORT trials
if use_all_short
    % Plot individual trials (thin) and average (thick)
    all_binned_short = zeros(length(short_trial_data), length(pos_centers));
    for trial_idx = 1:length(short_trial_data)
        binned_rate = zeros(size(pos_centers));
        binned_count = zeros(size(pos_centers));
        pos = short_trial_data{trial_idx}.position_absolute;
        fr = short_trial_data{trial_idx}.firing_rate;
        for i = 1:length(pos)
            bin_idx = find(pos(i) >= pos_bins(1:end-1) & pos(i) < pos_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_binned_short(trial_idx, :) = binned_rate;
        plot(pos_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_short, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_binned_short = mean(all_binned_short, 1);
    valid_avg_short = avg_binned_short > 0;
    plot(pos_centers(valid_avg_short), avg_binned_short(valid_avg_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2.5, 'DisplayName', 'Short (avg)');
else
    % Single trial
    binned_rate_short = zeros(size(pos_centers));
    binned_count_short = zeros(size(pos_centers));
    for i = 1:length(position_short)
        pos_absolute = position_short_absolute(i);
        bin_idx = find(pos_absolute >= pos_bins(1:end-1) & pos_absolute < pos_bins(2:end), 1);
        if ~isempty(bin_idx)
            binned_rate_short(bin_idx) = binned_rate_short(bin_idx) + firing_rate_short(i);
            binned_count_short(bin_idx) = binned_count_short(bin_idx) + 1;
        end
    end
    valid_bins_short = binned_count_short > 0;
    binned_rate_short(valid_bins_short) = binned_rate_short(valid_bins_short) ./ binned_count_short(valid_bins_short);
    plot(pos_centers(valid_bins_short), binned_rate_short(valid_bins_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');
end

xlabel('Absolute Position (cm)');
ylabel('Firing Rate (Hz)');
title('Absolute Space Tuning');
xlim([0, 120]);
grid on;
hold off;

% Panel 5 (Row 2, Col 2): Relative space tuning
subplot(3, 6, 10:12);
hold on;

% Bin firing rate by relative position (0-100%)
rel_bin_size = 2;  % 2% bins for relative position
rel_bins = 0:rel_bin_size:100;
rel_centers = rel_bins(1:end-1) + rel_bin_size/2;

% Process LONG trials
if use_all_long
    % Plot individual trials (thin) and average (thick)
    all_rel_binned_long = zeros(length(long_trial_data), length(rel_centers));
    for trial_idx = 1:length(long_trial_data)
        binned_rate = zeros(size(rel_centers));
        binned_count = zeros(size(rel_centers));
        pos = long_trial_data{trial_idx}.position;
        rel_pos = (pos - min(pos)) / (max(pos) - min(pos)) * 100;
        fr = long_trial_data{trial_idx}.firing_rate;
        for i = 1:length(rel_pos)
            bin_idx = find(rel_pos(i) >= rel_bins(1:end-1) & rel_pos(i) < rel_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_rel_binned_long(trial_idx, :) = binned_rate;
        plot(rel_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_long, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_rel_binned_long = mean(all_rel_binned_long, 1);
    valid_avg_long = avg_rel_binned_long > 0;
    plot(rel_centers(valid_avg_long), avg_rel_binned_long(valid_avg_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2.5, 'DisplayName', 'Long (avg)');
else
    % Single trial
    rel_pos_long = (position_long - min(position_long)) / (max(position_long) - min(position_long)) * 100;
    rel_binned_rate_long = zeros(size(rel_centers));
    rel_binned_count_long = zeros(size(rel_centers));
    for i = 1:length(rel_pos_long)
        bin_idx = find(rel_pos_long(i) >= rel_bins(1:end-1) & rel_pos_long(i) < rel_bins(2:end), 1);
        if ~isempty(bin_idx)
            rel_binned_rate_long(bin_idx) = rel_binned_rate_long(bin_idx) + firing_rate_long(i);
            rel_binned_count_long(bin_idx) = rel_binned_count_long(bin_idx) + 1;
        end
    end
    rel_valid_bins_long = rel_binned_count_long > 0;
    rel_binned_rate_long(rel_valid_bins_long) = rel_binned_rate_long(rel_valid_bins_long) ./ rel_binned_count_long(rel_valid_bins_long);
    plot(rel_centers(rel_valid_bins_long), rel_binned_rate_long(rel_valid_bins_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
end

% Process SHORT trials
if use_all_short
    % Plot individual trials (thin) and average (thick)
    all_rel_binned_short = zeros(length(short_trial_data), length(rel_centers));
    for trial_idx = 1:length(short_trial_data)
        binned_rate = zeros(size(rel_centers));
        binned_count = zeros(size(rel_centers));
        pos = short_trial_data{trial_idx}.position;
        rel_pos = (pos - min(pos)) / (max(pos) - min(pos)) * 100;
        fr = short_trial_data{trial_idx}.firing_rate;
        for i = 1:length(rel_pos)
            bin_idx = find(rel_pos(i) >= rel_bins(1:end-1) & rel_pos(i) < rel_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_rel_binned_short(trial_idx, :) = binned_rate;
        plot(rel_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_short, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_rel_binned_short = mean(all_rel_binned_short, 1);
    valid_avg_short = avg_rel_binned_short > 0;
    plot(rel_centers(valid_avg_short), avg_rel_binned_short(valid_avg_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2.5, 'DisplayName', 'Short (avg)');
else
    % Single trial
    rel_pos_short = (position_short - min(position_short)) / (max(position_short) - min(position_short)) * 100;
    rel_binned_rate_short = zeros(size(rel_centers));
    rel_binned_count_short = zeros(size(rel_centers));
    for i = 1:length(rel_pos_short)
        bin_idx = find(rel_pos_short(i) >= rel_bins(1:end-1) & rel_pos_short(i) < rel_bins(2:end), 1);
        if ~isempty(bin_idx)
            rel_binned_rate_short(bin_idx) = rel_binned_rate_short(bin_idx) + firing_rate_short(i);
            rel_binned_count_short(bin_idx) = rel_binned_count_short(bin_idx) + 1;
        end
    end
    rel_valid_bins_short = rel_binned_count_short > 0;
    rel_binned_rate_short(rel_valid_bins_short) = rel_binned_rate_short(rel_valid_bins_short) ./ rel_binned_count_short(rel_valid_bins_short);
    plot(rel_centers(rel_valid_bins_short), rel_binned_rate_short(rel_valid_bins_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');
end

xlabel('Relative Position (%)');
ylabel('Firing Rate (Hz)');
title('Relative Space Tuning');
xlim([0, 100]);
grid on;
hold off;

% Panel 6 (Row 3, Col 1): Velocity tuning
subplot(3, 6, 13:15);
hold on;

% Determine common velocity bins across both trials from the actual data
if use_all_long && use_all_short
    all_vels = [];
    for i = 1:length(long_trial_data)
        all_vels = [all_vels; long_trial_data{i}.velocity(:)]; %#ok<AGROW>
    end
    for i = 1:length(short_trial_data)
        all_vels = [all_vels; short_trial_data{i}.velocity(:)]; %#ok<AGROW>
    end
else
    all_vels = [velocity_long(:); velocity_short(:)];
end
vel_bin_size = 5;  % cm/s bins
vel_bins = floor(min(all_vels)):vel_bin_size:ceil(max(all_vels));
vel_centers = vel_bins(1:end-1) + vel_bin_size/2;

% Process LONG trials
if use_all_long
    % Plot individual trials (thin) and average (thick)
    all_vel_binned_long = zeros(length(long_trial_data), length(vel_centers));
    for trial_idx = 1:length(long_trial_data)
        binned_rate = zeros(size(vel_centers));
        binned_count = zeros(size(vel_centers));
        vel = long_trial_data{trial_idx}.velocity;
        fr = long_trial_data{trial_idx}.firing_rate;
        for i = 1:length(vel)
            bin_idx = find(vel(i) >= vel_bins(1:end-1) & vel(i) < vel_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_vel_binned_long(trial_idx, :) = binned_rate;
        plot(vel_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_long, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_vel_binned_long = mean(all_vel_binned_long, 1);
    valid_avg_long = avg_vel_binned_long > 0;
    plot(vel_centers(valid_avg_long), avg_vel_binned_long(valid_avg_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2.5, 'DisplayName', 'Long (avg)');
else
    % Single trial
    vel_binned_rate_long = zeros(size(vel_centers));
    vel_binned_count_long = zeros(size(vel_centers));
    for i = 1:length(velocity_long)
        bin_idx = find(velocity_long(i) >= vel_bins(1:end-1) & velocity_long(i) < vel_bins(2:end), 1);
        if ~isempty(bin_idx)
            vel_binned_rate_long(bin_idx) = vel_binned_rate_long(bin_idx) + firing_rate_long(i);
            vel_binned_count_long(bin_idx) = vel_binned_count_long(bin_idx) + 1;
        end
    end
    vel_valid_bins_long = vel_binned_count_long > 0;
    vel_binned_rate_long(vel_valid_bins_long) = vel_binned_rate_long(vel_valid_bins_long) ./ vel_binned_count_long(vel_valid_bins_long);
    plot(vel_centers(vel_valid_bins_long), vel_binned_rate_long(vel_valid_bins_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
end

% Process SHORT trials
if use_all_short
    % Plot individual trials (thin) and average (thick)
    all_vel_binned_short = zeros(length(short_trial_data), length(vel_centers));
    for trial_idx = 1:length(short_trial_data)
        binned_rate = zeros(size(vel_centers));
        binned_count = zeros(size(vel_centers));
        vel = short_trial_data{trial_idx}.velocity;
        fr = short_trial_data{trial_idx}.firing_rate;
        for i = 1:length(vel)
            bin_idx = find(vel(i) >= vel_bins(1:end-1) & vel(i) < vel_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_vel_binned_short(trial_idx, :) = binned_rate;
        plot(vel_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_short, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_vel_binned_short = mean(all_vel_binned_short, 1);
    valid_avg_short = avg_vel_binned_short > 0;
    plot(vel_centers(valid_avg_short), avg_vel_binned_short(valid_avg_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2.5, 'DisplayName', 'Short (avg)');
else
    % Single trial
    vel_binned_rate_short = zeros(size(vel_centers));
    vel_binned_count_short = zeros(size(vel_centers));
    for i = 1:length(velocity_short)
        bin_idx = find(velocity_short(i) >= vel_bins(1:end-1) & velocity_short(i) < vel_bins(2:end), 1);
        if ~isempty(bin_idx)
            vel_binned_rate_short(bin_idx) = vel_binned_rate_short(bin_idx) + firing_rate_short(i);
            vel_binned_count_short(bin_idx) = vel_binned_count_short(bin_idx) + 1;
        end
    end
    vel_valid_bins_short = vel_binned_count_short > 0;
    vel_binned_rate_short(vel_valid_bins_short) = vel_binned_rate_short(vel_valid_bins_short) ./ vel_binned_count_short(vel_valid_bins_short);
    plot(vel_centers(vel_valid_bins_short), vel_binned_rate_short(vel_valid_bins_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');
end

xlabel('Velocity (cm/s)');
ylabel('Firing Rate (Hz)');
title('Velocity Tuning');
xlim([0, 30]);
grid on;
hold off;

% Panel 7 (Row 3, Col 2): Acceleration tuning
subplot(3, 6, 16:18);
hold on;

% Determine common acceleration bins across both trials from the actual data
if use_all_long && use_all_short
    all_accels = [];
    for i = 1:length(long_trial_data)
        all_accels = [all_accels; long_trial_data{i}.acceleration(:)]; %#ok<AGROW>
    end
    for i = 1:length(short_trial_data)
        all_accels = [all_accels; short_trial_data{i}.acceleration(:)]; %#ok<AGROW>
    end
else
    all_accels = [acceleration_long(:); acceleration_short(:)];
end
accel_bin_size = 0.5;  % cm/s² bins
accel_bins = floor(min(all_accels)):accel_bin_size:ceil(max(all_accels));
accel_centers = accel_bins(1:end-1) + accel_bin_size/2;

% Process LONG trials
if use_all_long
    % Plot individual trials (thin) and average (thick)
    all_accel_binned_long = zeros(length(long_trial_data), length(accel_centers));
    for trial_idx = 1:length(long_trial_data)
        binned_rate = zeros(size(accel_centers));
        binned_count = zeros(size(accel_centers));
        accel = long_trial_data{trial_idx}.acceleration;
        fr = long_trial_data{trial_idx}.firing_rate;
        for i = 1:length(accel)
            bin_idx = find(accel(i) >= accel_bins(1:end-1) & accel(i) < accel_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_accel_binned_long(trial_idx, :) = binned_rate;
        plot(accel_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_long, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_accel_binned_long = mean(all_accel_binned_long, 1);
    valid_avg_long = avg_accel_binned_long > 0;
    plot(accel_centers(valid_avg_long), avg_accel_binned_long(valid_avg_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2.5, 'DisplayName', 'Long (avg)');
else
    % Single trial
    accel_binned_rate_long = zeros(size(accel_centers));
    accel_binned_count_long = zeros(size(accel_centers));
    for i = 1:length(acceleration_long)
        bin_idx = find(acceleration_long(i) >= accel_bins(1:end-1) & acceleration_long(i) < accel_bins(2:end), 1);
        if ~isempty(bin_idx)
            accel_binned_rate_long(bin_idx) = accel_binned_rate_long(bin_idx) + firing_rate_long(i);
            accel_binned_count_long(bin_idx) = accel_binned_count_long(bin_idx) + 1;
        end
    end
    accel_valid_bins_long = accel_binned_count_long > 0;
    accel_binned_rate_long(accel_valid_bins_long) = accel_binned_rate_long(accel_valid_bins_long) ./ accel_binned_count_long(accel_valid_bins_long);
    plot(accel_centers(accel_valid_bins_long), accel_binned_rate_long(accel_valid_bins_long), 'o-', ...
        'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
end

% Process SHORT trials
if use_all_short
    % Plot individual trials (thin) and average (thick)
    all_accel_binned_short = zeros(length(short_trial_data), length(accel_centers));
    for trial_idx = 1:length(short_trial_data)
        binned_rate = zeros(size(accel_centers));
        binned_count = zeros(size(accel_centers));
        accel = short_trial_data{trial_idx}.acceleration;
        fr = short_trial_data{trial_idx}.firing_rate;
        for i = 1:length(accel)
            bin_idx = find(accel(i) >= accel_bins(1:end-1) & accel(i) < accel_bins(2:end), 1);
            if ~isempty(bin_idx)
                binned_rate(bin_idx) = binned_rate(bin_idx) + fr(i);
                binned_count(bin_idx) = binned_count(bin_idx) + 1;
            end
        end
        valid_bins = binned_count > 0;
        binned_rate(valid_bins) = binned_rate(valid_bins) ./ binned_count(valid_bins);
        all_accel_binned_short(trial_idx, :) = binned_rate;
        plot(accel_centers(valid_bins), binned_rate(valid_bins), '-', ...
            'Color', [color_short, 0.3], 'LineWidth', 0.5);
    end
    % Average
    avg_accel_binned_short = mean(all_accel_binned_short, 1);
    valid_avg_short = avg_accel_binned_short > 0;
    plot(accel_centers(valid_avg_short), avg_accel_binned_short(valid_avg_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2.5, 'DisplayName', 'Short (avg)');
else
    % Single trial
    accel_binned_rate_short = zeros(size(accel_centers));
    accel_binned_count_short = zeros(size(accel_centers));
    for i = 1:length(acceleration_short)
        bin_idx = find(acceleration_short(i) >= accel_bins(1:end-1) & acceleration_short(i) < accel_bins(2:end), 1);
        if ~isempty(bin_idx)
            accel_binned_rate_short(bin_idx) = accel_binned_rate_short(bin_idx) + firing_rate_short(i);
            accel_binned_count_short(bin_idx) = accel_binned_count_short(bin_idx) + 1;
        end
    end
    accel_valid_bins_short = accel_binned_count_short > 0;
    accel_binned_rate_short(accel_valid_bins_short) = accel_binned_rate_short(accel_valid_bins_short) ./ accel_binned_count_short(accel_valid_bins_short);
    plot(accel_centers(accel_valid_bins_short), accel_binned_rate_short(accel_valid_bins_short), 'o-', ...
        'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');
end

xlabel('Acceleration (cm/s²)');
ylabel('Firing Rate (Hz)');
title('Acceleration Tuning');
xlim([-5, 5]);
grid on;
hold off;

% Overall title
long_str = char(string(trial_index_long));
short_str = char(string(trial_index_short));
sgtitle(sprintf('Synthetic Tuning Visualization - %s (Long: %s, Short: %s)', ...
    strrep(probe_id, '_', '\\_'), long_str, short_str), ...
    'FontSize', 12, 'FontWeight', 'bold');

fprintf('\nFigure created successfully!\n');

%% ==================== HELPER FUNCTIONS ====================

function firing_rate = generate_spatial_tuning(position, params)
    % Generate firing rate from spatial position
    switch params.type
        case 'gaussian'
            firing_rate = params.baseline + params.amplitude * ...
                exp(-((position - params.mu).^2) / (2 * params.sigma^2));
            
        case 'linear'
            % Linear increase from start to end
            firing_rate = params.baseline + params.amplitude * ...
                (position - min(position)) / (max(position) - min(position));
            
        case 'flat'
            firing_rate = params.baseline * ones(size(position));
            
        otherwise
            error('Unknown spatial tuning type: %s', params.type);
    end
end

function firing_rate = generate_velocity_tuning(velocity, params)
    % Generate firing rate from velocity
    switch params.type
        case 'gaussian'
            firing_rate = params.baseline + params.amplitude * ...
                exp(-((velocity - params.mu).^2) / (2 * params.sigma^2));
            
        case 'linear_positive'
            firing_rate = params.baseline + params.slope * velocity;
            
        case 'linear_negative'
            firing_rate = params.baseline - params.slope * velocity;
            
        case 'sigmoid'
            firing_rate = params.baseline + params.amplitude ./ ...
                (1 + exp(-params.steepness * (velocity - params.midpoint)));
            
        otherwise
            error('Unknown velocity tuning type: %s', params.type);
    end
    
    % Ensure non-negative
    firing_rate(firing_rate < 0) = 0;
end

function firing_rate = generate_acceleration_tuning(acceleration, params)
    % Generate firing rate from acceleration
    switch params.type
        case 'linear_positive'
            firing_rate = params.baseline + params.slope * acceleration;
            
        case 'linear_negative'
            firing_rate = params.baseline - params.slope * acceleration;
            
        case 'sigmoid'
            firing_rate = params.baseline + params.amplitude ./ ...
                (1 + exp(-params.steepness * (acceleration - params.midpoint)));
            
        case 'u_shaped_positive'
            % U-shape: high at extremes, low at center
            firing_rate = params.baseline + params.amplitude * ...
                params.curvature * (acceleration - params.midpoint).^2;
            
        case 'u_shaped_negative'
            % Inverted U-shape: low at extremes, high at center
            max_rate = params.baseline + params.amplitude;
            firing_rate = max_rate - params.amplitude * ...
                params.curvature * (acceleration - params.midpoint).^2;
            
        otherwise
            error('Unknown acceleration tuning type: %s', params.type);
    end
    
    % Ensure non-negative
    firing_rate(firing_rate < 0) = 0;
end
