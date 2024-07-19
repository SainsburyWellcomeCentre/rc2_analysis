% Statistics on the the motion and baseline durations

% initialize controller and variables
ctl = RC2Analysis();
probe_ids  = ctl.get_probe_ids('visual_flow', ...
                               'mismatch_nov20', ...
                               'mismatch_jul21', ...
                               'mismatch_darkness_oct21', ...
                               'darkness', ...
                               'passive_same_luminance', ...
                               'passive_same_luminance_muscimol');
trial_types                 = {'RT', 'RT_gain_up', 'R', 'T_bank', 'T_RT', 'T_R', 'T', 'RVT', 'RVT_gain_up','RV', 'RV_gain_up', 'VT_RVT', 'VT_RV', 'V_RVT', 'V_RV', 'T_Vstatic', 'V', 'VT'};

stationary_times = [];
motion_times = [];

% Collect stationary times for every probe
for probe_i = 1 : length(probe_ids)
    data   = ctl.load_formatted_data(probe_ids{probe_i});

    % for every condition
    for type_i = 1 : length(trial_types)
        trials = data.get_trials_with_trial_group_label(trial_types{type_i});

        % for every trial
        for trial_i = 1 : length(trials)
            trial  = trials{trial_i}.to_aligned;
            
            % Exclude artifacts coming from jittering
            if trial.stationary_time < 3
                stationary_times(end+1) = 3;
            else
                stationary_times(end+1) = trial.stationary_time;
            end

            motion_times(end+1) = trial.motion_time;
        end
    end
end


% Calculate and print all relevant statistics
stationary_times_mean = mean(stationary_times);
stationary_times_median = median(stationary_times);
stationary_times_std = std(stationary_times);
stationary_times_max = max(stationary_times);
stationary_times_min = min(stationary_times);
stationary_times_q1 = prctile(stationary_times, 25);
stationary_times_q3 = prctile(stationary_times, 75);


motion_times_mean = mean(motion_times);
motion_times_median = median(motion_times);
motion_times_std = std(motion_times);
motion_times_max = max(motion_times);
motion_times_min = min(motion_times);
motion_times_q1 = prctile(motion_times, 25);
motion_times_q3 = prctile(motion_times, 75);

[p] = signrank(stationary_times, motion_times);

sprintf('Stationary times. Mean: %.2f, std: %.2f, max: %.2f, min: %.2f.', stationary_times_mean, stationary_times_std, stationary_times_max, stationary_times_min)
sprintf('Median: %.2f, q1: %.2f, q3: %.2f', stationary_times_median, stationary_times_q1, stationary_times_q3)
sprintf('Motion times. Mean: %.2f, std: %.2f, max: %.2f, min: %.2f.', motion_times_mean, motion_times_std, motion_times_max, motion_times_min)
sprintf('Median: %.2f, q1: %.2f, q3: %.2f', motion_times_median, motion_times_q1, motion_times_q3)
sprintf('p: %.5f', p)



