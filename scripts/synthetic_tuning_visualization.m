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
probe_id = 'CAA-1123432a_rec1';  % Which probe's data to use
trial_index_long = 1;   % Which long trial's velocity profile to use (1-based index)
trial_index_short = 2;  % Which short trial's velocity profile to use (1-based index)

% --- Primary Tuning Type ---
% Choose which tuning characteristic to start from
% Options: 'absolute_space', 'relative_space', 'velocity', 'acceleration'
primary_tuning = 'relative_space';

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
vel_tuning.mu = 10;            % Peak velocity (cm/s) for Gaussian
vel_tuning.sigma = 10;         % Width for Gaussian (cm/s)
vel_tuning.amplitude = 50;     % Peak firing rate (Hz)
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

% Validate trial indices
if trial_index_long > length(long_trials)
    error('Long trial index %d exceeds number of available long trials (%d)', ...
        trial_index_long, length(long_trials));
end
if trial_index_short > length(short_trials)
    error('Short trial index %d exceeds number of available short trials (%d)', ...
        trial_index_short, length(short_trials));
end

% Extract data for LONG trial
trial_long = long_trials(trial_index_long);
motion_mask_long = trial_long.motion_mask();
position_long = trial_long.position(motion_mask_long);
velocity_long = trial_long.velocity(motion_mask_long);
acceleration_long = trial_long.acceleration(motion_mask_long);  % Use real acceleration from trial
time_long = trial_long.probe_t(motion_mask_long);

fprintf('Loaded LONG trial %d: %d samples, duration %.2f s\n', ...
    trial_index_long, length(time_long), time_long(end) - time_long(1));
fprintf('  Position: %.1f - %.1f cm, Velocity: %.1f - %.1f cm/s, Accel: %.1f - %.1f cm/s²\n', ...
    min(position_long), max(position_long), min(velocity_long), max(velocity_long), ...
    min(acceleration_long), max(acceleration_long));

% Extract data for SHORT trial
trial_short = short_trials(trial_index_short);
motion_mask_short = trial_short.motion_mask();
position_short = trial_short.position(motion_mask_short);
velocity_short = trial_short.velocity(motion_mask_short);
acceleration_short = trial_short.acceleration(motion_mask_short);  % Use real acceleration from trial
time_short = trial_short.probe_t(motion_mask_short);

fprintf('Loaded SHORT trial %d: %d samples, duration %.2f s\n', ...
    trial_index_short, length(time_short), time_short(end) - time_short(1));
fprintf('  Position: %.1f - %.1f cm, Velocity: %.1f - %.1f cm/s, Accel: %.1f - %.1f cm/s²\n', ...
    min(position_short), max(position_short), min(velocity_short), max(velocity_short), ...
    min(acceleration_short), max(acceleration_short));

% For absolute space analysis, short trials start at 60 cm on the track
% Add offset to convert to absolute track coordinates (long: 0-120, short: 60-120)
short_trial_start_position = 60;  % cm
position_short_absolute = position_short + short_trial_start_position;
fprintf('  Short trial absolute positions: %.1f - %.1f cm (offset by %.1f cm)\n', ...
    min(position_short_absolute), max(position_short_absolute), short_trial_start_position);

%% ==================== GENERATE SYNTHETIC FIRING RATES ====================

% Generate firing rates for LONG trial
switch primary_tuning
    case 'absolute_space'
        % Use actual position values - same tuning curve for both trials
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

% Generate firing rates for SHORT trial
switch primary_tuning
    case 'absolute_space'
        % Use absolute position values (with offset) - same tuning curve for both trials
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

% Add noise to make it more realistic
noise_level = 0.1;  % 10% noise
firing_rate_long = firing_rate_long .* (1 + noise_level * randn(size(firing_rate_long)));
firing_rate_long(firing_rate_long < 0) = 0;  % Ensure non-negative
firing_rate_short = firing_rate_short .* (1 + noise_level * randn(size(firing_rate_short)));
firing_rate_short(firing_rate_short < 0) = 0;  % Ensure non-negative

%% ==================== CREATE FIGURE ====================

fig = figure('Position', [100, 100, 1200, 900], 'Color', 'w');

% Color scheme
color_long = [0, 0.4470, 0.7410];    % Blue
color_short = [0.8500, 0.3250, 0.0980];  % Red

% Panel 1 (Row 1, Col 1): Parameters
subplot(3, 6, 1:2);
axis off;
param_text = {
    '\bf{Synthetic Tuning Parameters}',
    '',
    sprintf('\\bf{Probe:} %s', strrep(probe_id, '_', '\_')),
    sprintf('\\bf{Long Trial:} %d, \\bf{Short Trial:} %d', trial_index_long, trial_index_short),
    sprintf('\\bf{Primary Tuning:} %s', strrep(primary_tuning, '_', ' ')),
    '',
    '\bf{Active Tuning Types:}'
};

if abs_space.enabled
    param_text{end+1} = sprintf('  • Absolute Space (%s, \\mu=%.0f cm, \\sigma=%.0f cm)', ...
        abs_space.type, abs_space.mu, abs_space.sigma);
end
if rel_space.enabled
    param_text{end+1} = sprintf('  • Relative Space (%s, \\mu=%.0f%%, \\sigma=%.0f%%)', ...
        rel_space.type, rel_space.mu*100, rel_space.sigma*100);
end
if vel_tuning.enabled
    param_text{end+1} = sprintf('  • Velocity (%s)', vel_tuning.type);
end
if accel_tuning.enabled
    param_text{end+1} = sprintf('  • Acceleration (%s)', accel_tuning.type);
end

text(0.05, 0.95, param_text, 'VerticalAlignment', 'top', 'FontSize', 9, ...
    'Interpreter', 'tex');

% Panel 2 (Row 1, Col 2): Long trial running profile
subplot(3, 6, 3:4);
hold on;
yyaxis left
plot(time_long, position_long, '-', 'Color', color_long, 'LineWidth', 1.5);
ylabel('Position (cm)');
ylim([0, 120]);

yyaxis right
plot(time_long, velocity_long, '-', 'Color', color_long, 'LineWidth', 1);
ylabel('Velocity (cm/s)');

xlabel('Time (s)');
title('Long Trial Running Profile');
grid on;
hold off;

% Panel 3 (Row 1, Col 3): Short trial running profile
subplot(3, 6, 5:6);
hold on;
yyaxis left
plot(time_short, position_short + 60, '-', 'Color', color_short, 'LineWidth', 1.5);
ylabel('Position (cm)');
ylim([0, 120]);

yyaxis right
plot(time_short, velocity_short, '-', 'Color', color_short, 'LineWidth', 1);
ylabel('Velocity (cm/s)');

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

% Bin firing rate by position for SHORT trial (use absolute positions with offset)
binned_rate_short = zeros(size(pos_centers));
binned_count_short = zeros(size(pos_centers));

for i = 1:length(position_short)
    pos_absolute = position_short_absolute(i);  % Use pre-computed absolute position
    bin_idx = find(pos_absolute >= pos_bins(1:end-1) & pos_absolute < pos_bins(2:end), 1);
    if ~isempty(bin_idx)
        binned_rate_short(bin_idx) = binned_rate_short(bin_idx) + firing_rate_short(i);
        binned_count_short(bin_idx) = binned_count_short(bin_idx) + 1;
    end
end
valid_bins_short = binned_count_short > 0;
binned_rate_short(valid_bins_short) = binned_rate_short(valid_bins_short) ./ binned_count_short(valid_bins_short);

% Plot both
plot(pos_centers(valid_bins_long), binned_rate_long(valid_bins_long), 'o-', ...
    'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
plot(pos_centers(valid_bins_short), binned_rate_short(valid_bins_short), 'o-', ...
    'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');

xlabel('Absolute Position (cm)');
ylabel('Firing Rate (Hz)');
title('Absolute Space Tuning');
legend('Location', 'best');
xlim([0, 120]);
max_rate = max([binned_rate_long(valid_bins_long), binned_rate_short(valid_bins_short)]);
ylim([0, max_rate * 1.2]);
grid on;
hold off;

% Panel 5 (Row 2, Col 2): Relative space tuning
subplot(3, 6, 10:12);
hold on;

% Bin firing rate by relative position (0-100%)
rel_bin_size = 2;  % 2% bins for relative position
rel_bins = 0:rel_bin_size:100;
rel_centers = rel_bins(1:end-1) + rel_bin_size/2;

% Bin firing rate by relative position for LONG trial
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

% Bin firing rate by relative position for SHORT trial
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

% Plot both
plot(rel_centers(rel_valid_bins_long), rel_binned_rate_long(rel_valid_bins_long), 'o-', ...
    'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
plot(rel_centers(rel_valid_bins_short), rel_binned_rate_short(rel_valid_bins_short), 'o-', ...
    'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');

xlabel('Relative Position (%)');
ylabel('Firing Rate (Hz)');
title('Relative Space Tuning');
legend('Location', 'best');
xlim([0, 100]);
max_rate = max([rel_binned_rate_long(rel_valid_bins_long), rel_binned_rate_short(rel_valid_bins_short)]);
ylim([0, max_rate * 1.2]);
grid on;
hold off;

% Panel 6 (Row 3, Col 1): Velocity tuning
subplot(3, 6, 13:15);
hold on;

% Determine common velocity bins across both trials from the actual data
all_velocities = [velocity_long(:); velocity_short(:)];
vel_bin_size = 5;  % cm/s bins
vel_bins = floor(min(all_velocities)):vel_bin_size:ceil(max(all_velocities));
vel_centers = vel_bins(1:end-1) + vel_bin_size/2;

% Bin firing rate by velocity for LONG trial
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

% Bin firing rate by velocity for SHORT trial
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

% Plot both
plot(vel_centers(vel_valid_bins_long), vel_binned_rate_long(vel_valid_bins_long), 'o-', ...
    'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
plot(vel_centers(vel_valid_bins_short), vel_binned_rate_short(vel_valid_bins_short), 'o-', ...
    'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');

xlabel('Velocity (cm/s)');
ylabel('Firing Rate (Hz)');
title('Velocity Tuning');
legend('Location', 'best');
if length(vel_centers) > 1
    xlim([min(vel_centers), max(vel_centers)]);
end
max_rate = max([vel_binned_rate_long(vel_valid_bins_long), vel_binned_rate_short(vel_valid_bins_short)]);
if max_rate > 0
    ylim([0, max_rate * 1.2]);
end
grid on;
hold off;

% Panel 7 (Row 3, Col 2): Acceleration tuning
subplot(3, 6, 16:18);
hold on;

% Determine common acceleration bins across both trials from the actual data
all_accelerations = [acceleration_long(:); acceleration_short(:)];
accel_bin_size = 0.5;  % cm/s² bins (should be in realistic ±5 cm/s² range)
accel_bins = floor(min(all_accelerations)):accel_bin_size:ceil(max(all_accelerations));
accel_centers = accel_bins(1:end-1) + accel_bin_size/2;

% Bin firing rate by acceleration for LONG trial
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

% Bin firing rate by acceleration for SHORT trial
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

% Plot both
plot(accel_centers(accel_valid_bins_long), accel_binned_rate_long(accel_valid_bins_long), 'o-', ...
    'Color', color_long, 'MarkerFaceColor', color_long, 'LineWidth', 2, 'DisplayName', 'Long');
plot(accel_centers(accel_valid_bins_short), accel_binned_rate_short(accel_valid_bins_short), 'o-', ...
    'Color', color_short, 'MarkerFaceColor', color_short, 'LineWidth', 2, 'DisplayName', 'Short');

xlabel('Acceleration (cm/s²)');
ylabel('Firing Rate (Hz)');
title('Acceleration Tuning');
legend('Location', 'best');
if length(accel_centers) > 1
    xlim([min(accel_centers), max(accel_centers)]);
end
max_rate = max([accel_binned_rate_long(accel_valid_bins_long), accel_binned_rate_short(accel_valid_bins_short)]);
if max_rate > 0
    ylim([0, max_rate * 1.2]);
end
grid on;
hold off;

% Overall title
sgtitle(sprintf('Synthetic Tuning Visualization - %s (Long: %d, Short: %d)', ...
    strrep(probe_id, '_', '\_'), trial_index_long, trial_index_short), ...
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
