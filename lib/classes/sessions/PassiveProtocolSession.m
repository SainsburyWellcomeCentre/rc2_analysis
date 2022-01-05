classdef PassiveProtocolSession < RVTSession
% PassiveProtocolSession Class for handling details of a
% passive_protocol session.
%
%   PassiveProtocolSession Properties:
%       trial_group_ids         - integer IDs for each trial group
%       trial_group_labels      - string IDs for each trial group
%
%   PassiveProtocolSession Methods:
%       get_trial_group_id      - Given a trial, return the trial group ID
%
%   See also: RVTSession   
    
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 3
        trial_group_labels    = {'VT', 'V', 'T'};
    end
    
    
    
    methods
        
        function obj = PassiveProtocolSession(session)
        % PassiveProtocolSession
        %
        %   PassiveProtocolSession(SESSION) takes a session structure,
        %   SESSION, saved in the formatted data file and creates an object
        %   to handle the data of a 'passive_protocol' session.
        
            obj = obj@RVTSession(session);
        end
        
        
        
        function group_id = get_trial_group_id(~, trial)
        %%get_trial_group_id Given a trial, return the trial group ID
        %
        %   GROUP_ID = get_trial_group_id(TRIAL)
        %   get which trial group ID the given trial in TRIAL, of class
        %   Trial, belongs to. 
        
            if strcmp(trial.protocol, 'StageOnly')
                if trial.config.enable_vis_stim
                    group_id = 1;
                else
                    group_id = 3;
                end
            elseif strcmp(trial.protocol, 'ReplayOnly')
                group_id = 2;
            end
        end
    end
end