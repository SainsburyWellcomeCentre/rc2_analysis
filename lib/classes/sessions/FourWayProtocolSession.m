classdef FourWayProtocolSession < RVTSession
% FourWayProtocolSession Class for handling details of a four_way_protocol
% session.
%
%   FourWayProtocolSession Properties:
%       trial_group_ids         - integer IDs for each trial group
%       trial_group_labels      - string IDs for each trial group
%
%   FourWayProtocolSession Methods:
%       get_trial_group_id      - Given a trial, return the trial group ID
%
%   See also: RVTSession

    properties (Constant = true)
        
        trial_group_ids      = 1 : 6
        trial_group_labels    = {'RVT', 'RV', 'VT_RVT', 'VT_RV', 'V_RVT', 'V_RV'};
    end
    
    
    
    methods
        
        function obj = FourWayProtocolSession(session)
        % FourWayProtocolSession
        %
        %   FourWayProtocolSession(SESSION) takes a session structure,
        %   SESSION, saved in the formatted data file and creates an object
        %   to handle the data of a 'four_way_protocol' session.
        
            obj = obj@RVTSession(session);
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
                else
                    group_id = 4;
                end
            elseif strcmp(trial.protocol, 'ReplayOnly')
                if strcmp(trial.replay_of, 'Coupled')
                    group_id = 5;
                else
                    group_id = 6;
                end
            end
        end
    end
end