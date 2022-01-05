classdef LocoVestLocoVestSession < RVTSession
% LocoVestLocoVestSession Class for handling details of a
% locovest_loco_vest session.
%
%   LocoVestLocoVestSession Properties:
%       trial_group_ids         - integer IDs for each trial group
%       trial_group_labels      - string IDs for each trial group
%
%   LocoVestLocoVestSession Methods:
%       get_trial_group_id      - Given a trial, return the trial group ID
%       remove_trials_without_motion - removes some trials without any motion
%
%   See also: RVTSession   
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 5
        trial_group_labels    = {'RT', 'R', 'T_bank', 'T_RT', 'T_R'};
    end
    
    
    
    methods
        
        function obj = LocoVestLocoVestSession(session)
        % LocoVestLocoVestSession
        %
        %   LocoVestLocoVestSession(SESSION) takes a session structure,
        %   SESSION, saved in the formatted data file and creates an object
        %   to handle the data of a 'locovest_loco_vest' session.
        
            obj = obj@RVTSession(session);
            
            for ii = 1 : session.n_trials
                
                % create trial object (may have to subclass this at one
                % point)
                obj.trials{ii} = Trial(session.trials(ii), obj);
                
                % get which group the trial belongs to
                trial_group_id = obj.get_trial_group_id(obj.trials{ii});
                
                % set the group id and label
                obj.trials{ii}.trial_group_id = trial_group_id;
                obj.trials{ii}.trial_group_label = obj.get_group_label_from_id(trial_group_id);
            end
            
            obj.remove_trials_without_motion();
        end
        
        
        
        function group_id = get_trial_group_id(~, trial)
        %%get_trial_group_id Given a trial, return the trial group ID
        %
        %   GROUP_ID = get_trial_group_id(TRIAL)
        %   get which trial group ID the given trial in TRIAL, of class
        %   Trial, belongs to. 
        
            if strcmp(trial.protocol, 'Coupled')
                group_id = 1;
            elseif strcmp(trial.protocol, 'EncoderOnly')
                group_id = 2;
            elseif strcmp(trial.protocol, 'StageOnly')
                if strcmp(trial.replay_of, 'Bank')
                    group_id = 3;
                elseif strcmp(trial.replay_of, 'Coupled')
                    group_id = 4;
                elseif strcmp(trial.replay_of, 'EncoderOnly')
                    group_id = 5;
                end
            end
        end
        
        
        
        function remove_trials_without_motion(obj, varargin)
        %%some RT trials didn't have motion... remove them so they don't
        %%influence later analysis
            
            trial_group_labels = cellfun(@(x)(x.trial_group_label), obj.trials, 'uniformoutput', false); %#ok<*PROPLC>
            motion_time = cellfun(@(x)(x.motion_time), obj.trials);
            
            idx = strcmp(trial_group_labels, 'RT') & motion_time < 0.0001;
            
            obj.trials(idx) = [];
        end
    end
end