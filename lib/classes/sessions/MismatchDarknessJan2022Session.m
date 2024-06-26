classdef MismatchDarknessJan2022Session < RVTSession
% MismatchDarknessJan2022Session Class for handling details of a
% mismatch_darkness_jan2022 session.
%
%   MismatchDarknessJan2022Session Properties:
%       trial_group_ids         - integer IDs for each trial group
%       trial_group_labels      - string IDs for each trial group
%
%   MismatchDarknessJan2022Session Methods:
%       get_trial_group_id      - Given a trial, return the trial group ID
%
%   See also: RVTSession

    properties (Constant = true)
        
        trial_group_ids      = 1 : 2
        trial_group_labels    = {'T_RT', 'RT_gain_up'};
    end
    
    
    
    methods
        
        function obj = MismatchDarknessJan2022Session(session)
        % MismatchDarknessJan2022Session
        %
        %   MismatchDarknessJan2022Session(SESSION) takes a session structure,
        %   SESSION, saved in the formatted data file and creates an object
        %   to handle the data of a 'mismatch_darkness_jan2021' session.
        
            obj = obj@RVTSession(session);
            
            for ii = 1 : session.n_trials
                
                % get which group the trial belongs to
                trial_group_id = obj.get_trial_group_id(session.trials(ii));
                
                % create trial object (may have to subclass this at one
                % point)
                obj.trials{ii} = Trial(session.trials(ii), obj);
                
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
        
            if strcmp(trial.protocol, 'StageOnly')
                group_id = 1;
            elseif strcmp(trial.protocol, 'CoupledMismatch')
                group_id = 2;
            end
        end
        
        
        
        function trials = get_trials_with_trial_group_label(obj, varargin)
        %%get_trials_with_trial_group_label Return cell array of Trial
        %%objects which are part of a trial group label
        %
        %   TRIALS = get_trials_with_trial_group_label(VARARGIN)
        %   reimplementation of get_trials_with_trial_group_label to remove
        %   some trials which were not true repeats
        %
        %   See also: RVTSession.get_trials_with_trial_group_label
        
            trials = get_trials_with_trial_group_label@RVTSession(obj, varargin{:});
            
            if strcmp(obj.probe_id, 'CAA-1116247_rec1')
                idx = cellfun(@(x)(x.trial_id), trials) == 2;
                if sum(idx) > 0
                    warning('Trial with ID 2 removed from data');
                    trials(idx) = [];
                end
            end
        end
    end
end
