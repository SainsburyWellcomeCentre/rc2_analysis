% Plot lick events per trial for each probe
%
% Creates event plots showing mouse licks with:
%   - Rows: one per probe
%   - Y-axis: trial number (with lick events marked)
%   - Column 1: Time (seconds from trial start)
%   - Column 2: Position (continuous space, cm)
%   - Column 3: Time-to-goal (continuous TTG, normalized 0-100%)
%
% The lick channel is a voltage signal that spikes when the mouse licks.
% Lick events are detected using a threshold crossing approach.
%
% Requires:
%   - Trial.lick property (available from Trial class)
%   - Trial.position() method for computing position
%   - compute_time_to_goal_tuning.m for TTG computation

%%
close all;

% Configuration
experiment_groups   = {'ambient_light'};
trial_group_labels  = {'RT'};
save_figs           = true;
overwrite           = true;
figure_dir          = {'lick_events', 'ambient_light'};

% Lick detection parameters
lick_threshold      = 1.0;  % Voltage threshold for lick detection (V)
min_lick_interval   = 0.05; % Minimum time between licks (s) to avoid double-counting

%%
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));

for pid = 1:length(probe_ids)
    fprintf('\nProcessing probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    
    % Load data
    data = ctl.load_formatted_data(probe_ids{pid});
    
    % Get all trials for the specified trial group
    all_trials = {};
    for kk = 1:length(trial_group_labels)
        trials = data.get_trials_with_trial_group_label(trial_group_labels{kk});
        all_trials = [all_trials; trials];
    end
    
    n_trials = length(all_trials);
    fprintf('  Found %d trials\n', n_trials);
    
    % --- Plot example raw lick signals for first probe only ---
    if pid == 1
        fprintf('  Creating example raw lick signal plots...\n');
        fig_examples = ctl.figs.a4figure('portrait');
        
        % Select up to 12 example trials
        n_examples = min(12, n_trials);
        example_indices = round(linspace(1, n_trials, n_examples));
        
        for ex_idx = 1:n_examples
            trial_idx = example_indices(ex_idx);
            trial = all_trials{trial_idx};
            
            % Get lick signal and time vector
            lick_signal = trial.lick;
            trial_time = trial.probe_t;
            trial_time_rel = trial_time - trial_time(1);  % Relative to trial start
            
            % Plot raw lick signal
            subplot(n_examples, 1, ex_idx);
            plot(trial_time_rel, lick_signal, 'k-', 'LineWidth', 0.5);
            hold on;
            
            % Draw threshold line
            yline(lick_threshold, 'r--', 'LineWidth', 1.5, 'Label', sprintf('Threshold = %.1fV', lick_threshold));
            
            % Detect and mark lick events
            lick_events_abs = detect_lick_events(lick_signal, lick_threshold, trial_time, min_lick_interval);
            if ~isempty(lick_events_abs)
                lick_events_rel = lick_events_abs - trial_time(1);
                for lick_t = lick_events_rel'
                    xline(lick_t, 'b-', 'LineWidth', 0.5, 'Alpha', 0.5);
                end
            end
            
            hold off;
            ylabel('Lick signal (V)');
            if ex_idx == n_examples
                xlabel('Time from trial start (s)');
            end
            title(sprintf('Trial %d (ID: %d)', trial_idx, trial.trial_id));
            grid on;
        end
        
        % Add figure title
        FigureTitle(fig_examples, sprintf('Example Raw Lick Signals - %s', probe_ids{pid}));
        
        % Save figure
        if save_figs
            ctl.figs.save_fig_to_join();
        end
    end
end

% --- Create combined figure with all probes (rows) and 3 columns ---
fprintf('\nCreating combined lick event plots for all probes...\n');
n_probes = length(probe_ids);
fig_combined = ctl.figs.a4figure('landscape');

% Track global max values for consistent axis limits
global_max_time = 0;
global_max_position = 0;

for pid = 1:length(probe_ids)
    fprintf('  Processing probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    
    % Load data
    data = ctl.load_formatted_data(probe_ids{pid});
    
    % Get all trials for the specified trial group
    all_trials = {};
    for kk = 1:length(trial_group_labels)
        trials = data.get_trials_with_trial_group_label(trial_group_labels{kk});
        all_trials = [all_trials; trials];
    end
    
    n_trials = length(all_trials);
    
    % Storage for all lick events across trials
    lick_times_per_trial = cell(n_trials, 1);
    lick_positions_per_trial = cell(n_trials, 1);
    lick_ttg_per_trial = cell(n_trials, 1);
    max_time = 0;
    max_position = 0;
    
    % Process each trial
    for trial_idx = 1:n_trials
        trial = all_trials{trial_idx};
        
        % Get lick signal and time vector
        lick_signal = trial.lick;
        trial_time = trial.probe_t;
        trial_time_rel = trial_time - trial_time(1);  % Relative to trial start
        
        % Detect lick events (threshold crossings)
        lick_events = detect_lick_events(lick_signal, lick_threshold, trial_time, min_lick_interval);
        
        if isempty(lick_events)
            continue;  % Skip trials with no licks
        end
        
        % Convert lick times to trial-relative times
        lick_times = lick_events - trial_time(1);
        lick_times_per_trial{trial_idx} = lick_times;
        
        % Compute position at lick times
        % Get motion mask and compute position
        motion_mask = trial.motion_mask();
        position = trial.position();  % This integrates velocity
        
        % Interpolate position at lick times
        lick_positions = interp1(trial_time_rel, position, lick_times, 'linear', nan);
        lick_positions_per_trial{trial_idx} = lick_positions;
        
        % Compute time-to-goal (TTG) at lick times
        % TTG is computed from motion samples only
        motion_times = trial.probe_t(motion_mask);
        n_motion_samples = length(motion_times);
        
        if n_motion_samples >= 2
            % Create continuous time from motion samples
            sampling_rate = 10000;  % Hz
            trial_duration = n_motion_samples / sampling_rate;
            continuous_time = linspace(0, trial_duration, n_motion_samples)';
            
            % Compute TTG (time remaining to end of trial)
            ttg_abs = trial_duration - continuous_time;
            
            % Normalize TTG to 0-100%
            ttg_norm = (ttg_abs / trial_duration) * 100;
            
            % Map lick times to motion times and get TTG
            lick_ttg = zeros(length(lick_times), 1);
            for lick_idx = 1:length(lick_times)
                % Find nearest motion sample to this lick time
                [~, nearest_idx] = min(abs(trial_time_rel - lick_times(lick_idx)));
                
                % Check if this sample is in motion
                if motion_mask(nearest_idx)
                    % Find index in motion samples
                    motion_idx = sum(motion_mask(1:nearest_idx));
                    if motion_idx > 0 && motion_idx <= length(ttg_norm)
                        lick_ttg(lick_idx) = ttg_norm(motion_idx);
                    else
                        lick_ttg(lick_idx) = nan;
                    end
                else
                    lick_ttg(lick_idx) = nan;
                end
            end
            
            lick_ttg_per_trial{trial_idx} = lick_ttg;
        else
            lick_ttg_per_trial{trial_idx} = nan(length(lick_times), 1);
        end
        
        % Track maximum values for axis limits
        if ~isempty(lick_times)
            max_time = max(max_time, max(lick_times));
            global_max_time = max(global_max_time, max(lick_times));
        end
        if ~isempty(lick_positions)
            max_position = max(max_position, max(lick_positions(~isnan(lick_positions))));
            global_max_position = max(global_max_position, max(lick_positions(~isnan(lick_positions))));
        end
    end
    
    % Plot 1: Time (column 1, row = pid)
    subplot(n_probes, 3, (pid-1)*3 + 1);
    hold on;
    for trial_idx = 1:n_trials
        lick_times = lick_times_per_trial{trial_idx};
        if ~isempty(lick_times)
            plot(lick_times, trial_idx * ones(size(lick_times)), 'k.', 'MarkerSize', 4);
        end
    end
    if pid == n_probes
        xlabel('Time from trial start (s)');
    end
    ylabel(sprintf('%s\nTrial #', probe_ids{pid}), 'Interpreter', 'none');
    if pid == 1
        title('Lick Events vs. Time');
    end
    ylim([0.5, n_trials + 0.5]);
    xlim([0, global_max_time * 1.05]);
    grid on;
    hold off;
    
    % Plot 2: Position (column 2, row = pid)
    subplot(n_probes, 3, (pid-1)*3 + 2);
    hold on;
    for trial_idx = 1:n_trials
        lick_positions = lick_positions_per_trial{trial_idx};
        if ~isempty(lick_positions)
            valid_idx = ~isnan(lick_positions);
            plot(lick_positions(valid_idx), trial_idx * ones(sum(valid_idx), 1), 'k.', 'MarkerSize', 4);
        end
    end
    if pid == n_probes
        xlabel('Position (cm)');
    end
    ylabel('Trial #');
    if pid == 1
        title('Lick Events vs. Position');
    end
    ylim([0.5, n_trials + 0.5]);
    xlim([0, global_max_position * 1.05]);
    grid on;
    hold off;
    
    % Plot 3: Time-to-goal (TTG, column 3, row = pid)
    subplot(n_probes, 3, (pid-1)*3 + 3);
    hold on;
    for trial_idx = 1:n_trials
        lick_ttg = lick_ttg_per_trial{trial_idx};
        if ~isempty(lick_ttg)
            valid_idx = ~isnan(lick_ttg);
            plot(lick_ttg(valid_idx), trial_idx * ones(sum(valid_idx), 1), 'k.', 'MarkerSize', 4);
        end
    end
    if pid == n_probes
        xlabel('Time-to-goal (%)');
    end
    ylabel('Trial #');
    if pid == 1
        title('Lick Events vs. Time-to-Goal');
    end
    ylim([0.5, n_trials + 0.5]);
    xlim([0, 100]);
    grid on;
    hold off;
end

% Add figure title
FigureTitle(fig_combined, sprintf('Lick Events - All Probes (%s)', strjoin(experiment_groups, ', ')));

% Save combined figure
if save_figs
    ctl.figs.save_fig_to_join();
end

% Join figures if requested
if save_figs
    fprintf('\nJoining figures into single PDF...\n');
    fname = sprintf('lick_events_%s.pdf', strjoin(experiment_groups, '_'));
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end

fprintf('\n=== All probes processed ===\n');

%% Helper function to detect lick events
function lick_events = detect_lick_events(lick_signal, threshold, time_vector, min_interval)
    % Detect lick events as threshold crossings
    %
    % Inputs:
    %   lick_signal   - Voltage signal from lick sensor
    %   threshold     - Voltage threshold for detection (V)
    %   time_vector   - Time vector corresponding to lick_signal
    %   min_interval  - Minimum time between consecutive licks (s)
    %
    % Outputs:
    %   lick_events   - Times of detected lick events
    
    % Find rising edge crossings (signal goes above threshold)
    above_threshold = lick_signal > threshold;
    rising_edges = find(diff([0; above_threshold]) == 1);
    
    if isempty(rising_edges)
        lick_events = [];
        return;
    end
    
    % Get times of rising edges
    lick_events = time_vector(rising_edges);
    
    % Remove events that are too close together (debouncing)
    if length(lick_events) > 1
        intervals = diff(lick_events);
        keep_idx = [true; intervals >= min_interval];
        lick_events = lick_events(keep_idx);
    end
end
