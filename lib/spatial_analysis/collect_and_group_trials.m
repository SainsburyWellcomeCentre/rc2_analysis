function [trial_groups, group_names, group_labels] = collect_and_group_trials(sessions, threshold_cm)
% COLLECT_AND_GROUP_TRIALS Collect trials from sessions and group by track length
%
%   [trial_groups, group_names, group_labels] = collect_and_group_trials(sessions)
%   [trial_groups, group_names, group_labels] = collect_and_group_trials(sessions, threshold_cm)
%
%   Collects all trials from the provided sessions and groups them into
%   'long' and 'short' categories based on the maximum in-motion position.
%
%   Inputs:
%       sessions     - Cell array of session objects
%       threshold_cm - (optional) Position threshold for long/short classification
%                      Default: 90 cm (trials > 90 cm are 'long')
%
%   Outputs:
%       trial_groups  - Struct with fields 'long' and 'short', each containing
%                       an array of structs with fields:
%                           .trial       - The aligned trial object
%                           .session_idx - Index of the source session
%                           .trial_idx   - Index within the session
%                           .max_pos     - Maximum in-motion position (cm)
%       group_names   - Cell array {'long', 'short'}
%       group_labels  - Cell array with descriptive labels for plotting
%
%   Example:
%       [groups, names, labels] = collect_and_group_trials(sessions);
%       long_trials = groups.long;
%       short_trials = groups.short;

    if nargin < 2
        threshold_cm = 90;
    end
    
    % Initialize structure to collect all trials
    all_trials_info = struct('trial', {}, 'session_idx', {}, 'trial_idx', {}, 'max_pos', {}, 'global_trial_idx', {});
    
    % Collect all trials with their metadata
    global_trial_counter = 0;
    for s = 1:length(sessions)
        session = sessions{s};
        trials_in_session = session.trials;
        
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            
            if ~isempty(pos)
                global_trial_counter = global_trial_counter + 1;
                max_pos = max(pos);
                all_trials_info(end+1) = struct(...
                    'trial', trial, ...
                    'session_idx', s, ...
                    'trial_idx', t, ...
                    'max_pos', max_pos, ...
                    'global_trial_idx', global_trial_counter); %#ok<AGROW>
            end
        end
    end
    
    % Group trials by max position
    long_trials = all_trials_info([all_trials_info.max_pos] > threshold_cm);
    short_trials = all_trials_info([all_trials_info.max_pos] <= threshold_cm);
    
    % Package results
    trial_groups = struct('long', long_trials, 'short', short_trials);
    group_names = {'long', 'short'};
    group_labels = {'Long (0-120 cm)', 'Short (60-120 cm)'};
end
