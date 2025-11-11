classdef LocoVestVariableDistanceSession < RVTSession
    % LocoVestVariableDistanceSession Class for managing running with translation sessions
    % with one recordings over 2 different distances.
    %
    %   This class handles the specific details of a 'locovest_variable_distance'
    %   session.
    %
    %   LocoVestVariableDistanceSession Properties:
    %       trial_group_ids         - Integer IDs representing each trial group
    %       trial_group_labels      - String labels for trial group is always RT (in ambient light).
    %
    %   LocoVestVariableDistanceSession Methods:
    %       get_trial_group_id      - Determines the trial group ID for a given trial
    %       LocoVestVariableDistanceSession - Constructor that initializes the
    %                                         session object and assigns trial
    %                                         group IDs and labels.
    %
    
        properties (Constant = true)
            
            trial_group_ids      = 1
            trial_group_labels    = {'RT'};
        end
        
        methods
            
            function obj = LocoVestVariableDistanceSession(session)
            % LocoVestVariableDistanceSession Constructor
            %
            %   obj = LocoVestVariableDistanceSession(SESSION) initializes a
            %   session structure, SESSION, to create an object for handling
            %   'locovest_variable_distance' sessions.
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
                group_id = 1;
            end
        end
    end
    
    