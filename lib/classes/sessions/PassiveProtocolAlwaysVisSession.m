classdef PassiveProtocolAlwaysVisSession < RVTSession
    % PassiveProtocolAlwaysVisSession Class for managing passive_protocol sessions
    % with a static visual stimulus for all trials.
    %
    %   This class handles the specific details of a 'passive_protocol'
    %   session in which visual stimuli are always present, although they
    %   may be static or dynamic based on the trial configuration.
    %
    %   PassiveProtocolAlwaysVisSession Properties:
    %       trial_group_ids         - Integer IDs representing each trial group
    %       trial_group_labels      - String labels for each trial group, including
    %                                 'VT' for visual and tactile stimuli, 'V' for
    %                                 visual-only, and 'T_Vstatic' for tactile
    %                                 with static visual stimuli.
    %
    %   PassiveProtocolAlwaysVisSession Methods:
    %       get_trial_group_id      - Determines the trial group ID for a given trial
    %       PassiveProtocolAlwaysVisSession - Constructor that initializes the
    %                                         session object and assigns trial
    %                                         group IDs and labels.
    %
    %   See also: RVTSession
    
        properties (Constant = true)
            
            trial_group_ids      = 1 : 3
            trial_group_labels    = {'VT', 'V', 'T_Vstatic'};
        end
        
        methods
            
            function obj = PassiveProtocolAlwaysVisSession(session)
            % PassiveProtocolAlwaysVisSession Constructor
            %
            %   obj = PassiveProtocolAlwaysVisSession(SESSION) initializes a
            %   session structure, SESSION, to create an object for handling
            %   'passive_protocol' sessions where visual stimuli are always enabled.
            %   This constructor assigns group IDs and labels to each trial in the session.
                obj = obj@RVTSession(session);
                
                for ii = 1 : session.n_trials
                    
                    % Get the group ID of the trial
                    trial_group_id = obj.get_trial_group_id(session.trials(ii));
                    
                    % Create a Trial object for each session trial
                    obj.trials{ii} = Trial(session.trials(ii), obj);
                    
                    % Set the group ID and label for each trial
                    obj.trials{ii}.trial_group_id = trial_group_id;
                    obj.trials{ii}.trial_group_label = obj.get_group_label_from_id(trial_group_id);
                end
            end
            
            function group_id = get_trial_group_id(~, trial)
            %%get_trial_group_id Determines the trial group ID for a given trial
            %
            %   group_id = get_trial_group_id(trial) returns the group ID for the
            %   specified trial, which is based on the trial's sequence configuration
            %   within the protocol. This ID corresponds to visual and tactile
            %   conditions defined in the session.
                group_id = trial.config.trial_sequence;
            end
        end
    end
    