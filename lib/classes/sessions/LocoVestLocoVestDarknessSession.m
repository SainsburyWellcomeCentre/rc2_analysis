classdef LocoVestLocoVestDarknessSession < RVTSession
% LocoVestLocoVestDarknessSession Class for handling details of a locovest_loco_vest_darkness
% session.
%
%   LocoVestLocoVestDarknessSession Properties:
%       trial_group_ids         - integer IDs for each trial group
%       trial_group_labels      - string IDs for each trial group
%
%   LocoVestLocoVestDarknessSession Methods:
%       get_trial_group_id      - Given a trial, return the trial group ID
%
%   See also: RVTSession    
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 4
        trial_group_labels    = {'RT', 'R', 'T_RT', 'T_R'}
    end
    
    
    
    methods
        
        function obj = LocoVestLocoVestDarknessSession(session)
        % LocoVestLocoVestDarknessSession
        %
        %   LocoVestLocoVestDarknessSession(SESSION) takes a session structure,
        %   SESSION, saved in the formatted data file and creates an object
        %   to handle the data of a 'locovest_loco_vest_darkness' session.
        
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
                if strcmp(trial.replay_of, 'Coupled')
                    group_id = 3;
                elseif strcmp(trial.replay_of, 'EncoderOnly')
                    group_id = 4;
                end
            end
        end
    end
end
