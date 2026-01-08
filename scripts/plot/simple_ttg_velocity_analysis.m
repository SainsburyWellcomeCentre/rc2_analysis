% Simple script to visualize running profile for one animal
% Creates 4 plots:
%   1. Velocity vs Time
%   2. Velocity vs Position (space bins)
%   3. Velocity vs Absolute Time-to-Goal (TTG)
%   4. Velocity vs Relative Time-to-Goal (TTG normalized by trial duration)
%
% Usage: Adjust experiment_groups and trial selection below

%%
% Configuration
experiment_groups = {'ambient_light'};  % Change this to your experiment group
trial_index = 1;  % Which trial to analyze (1 = first trial)

% Initialize controller
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});

fprintf('Found %d probe(s) for experiment group: %s\n', length(probe_ids), experiment_groups{1});

% Use first probe
if isempty(probe_ids)
    error('No probes found for experiment group: %s', experiment_groups{1});
end

probe_id = probe_ids{1};
fprintf('Using probe: %s\n', probe_id);

% Load data
data = ctl.load_formatted_data(probe_id);
sessions = data.motion_sessions();

if isempty(sessions)
    error('No motion sessions found for probe: %s', probe_id);
end

% Get trials from first session
session = sessions(1);
trials = session.trials();

if isempty(trials)
    error('No trials found in session');
end

if trial_index > length(trials)
    fprintf('Warning: trial_index %d exceeds number of trials (%d), using last trial\n', ...
        trial_index, length(trials));
    trial_index = length(trials);
end

trial = trials(trial_index);
fprintf('Analyzing trial %d/%d\n', trial_index, length(trials));

% Extract data during motion
motion_mask = trial.motion_mask();
time = trial.probe_t(motion_mask);
position = trial.position(motion_mask);
velocity = trial.velocity(motion_mask);

% Compute time-to-goal (absolute and relative)
trial_start = trial.probe_t(1);
trial_end = trial.probe_t(end);
trial_duration = trial_end - trial_start;

ttg_absolute = trial_end - time;  % seconds until trial end
ttg_relative = (ttg_absolute / trial_duration) * 100;  % percentage of trial remaining

% Bin position data (2 cm bins)
bin_size_cm = 2;
position_min = min(position);
position_max = max(position);
position_edges = position_min:bin_size_cm:position_max;
position_bins = (position_edges(1:end-1) + position_edges(2:end)) / 2;

% Compute average velocity per position bin
velocity_per_bin = nan(size(position_bins));
for i = 1:length(position_bins)
    if i < length(position_edges)
        in_bin = position >= position_edges(i) & position < position_edges(i+1);
        if any(in_bin)
            velocity_per_bin(i) = mean(velocity(in_bin));
        end
    end
end

% Create figure with 4 subplots
figure('Position', [100, 100, 1200, 900]);

% Plot 1: Velocity vs Time
subplot(2, 2, 1);
plot(time - trial_start, velocity, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (cm/s)');
title('Velocity vs Time');
grid on;

% Plot 2: Velocity vs Position
subplot(2, 2, 2);
plot(position_bins, velocity_per_bin, 'b-', 'LineWidth', 1.5);
xlabel('Position (cm)');
ylabel('Velocity (cm/s)');
title('Velocity vs Position');
grid on;

% Plot 3: Velocity vs Absolute TTG
subplot(2, 2, 3);
plot(ttg_absolute, velocity, 'r-', 'LineWidth', 1.5);
xlabel('Absolute Time-to-Goal (s)');
ylabel('Velocity (cm/s)');
title('Velocity vs Absolute TTG');
grid on;
set(gca, 'XDir', 'reverse');  % Reverse x-axis so goal is on the right

% Plot 4: Velocity vs Relative TTG
subplot(2, 2, 4);
plot(ttg_relative, velocity, 'm-', 'LineWidth', 1.5);
xlabel('Relative Time-to-Goal (%)');
ylabel('Velocity (cm/s)');
title('Velocity vs Relative TTG');
grid on;
set(gca, 'XDir', 'reverse');  % Reverse x-axis so goal is on the right

% Add overall title
sgtitle(sprintf('Running Profile - Probe: %s, Trial %d', probe_id, trial_index), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nPlot created successfully!\n');
fprintf('Trial info:\n');
fprintf('  Duration: %.2f seconds\n', trial_duration);
fprintf('  Position range: %.1f - %.1f cm\n', position_min, position_max);
fprintf('  Velocity range: %.1f - %.1f cm/s\n', min(velocity), max(velocity));
fprintf('  Mean velocity: %.1f cm/s\n', mean(velocity));
