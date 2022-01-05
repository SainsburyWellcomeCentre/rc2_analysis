classdef MotionBout < handle
% MotionBout Class for handling data about motion bouts during a trial
%
%   MotionBout Properties:
%       start_idx       - sample point of the trial on which motion starts
%       end_idx         - sample point of the trial on which motion ends
%       trial           - the Trial object with which the motion is associated
%       prepad          - time in s before motion onset to include in velocity trace for bout
%       postpad         - time in s after motion offset to include in velocity trace for bout
%       duration        - duration in s of the bout from start to end
%       start_time      - start of the motion bout in seconds (in probe time)
%       end_time        - end of the motion bout in seconds (in probe time)
%       velocity        - get the velocity of the motion bout from the Trial object
%
%   See also: Trial

    properties
        
        start_idx
        end_idx
        trial
        
        prepad = 0
        postpad = 0
    end
    
    properties (Dependent = true)
        duration
        start_time
        end_time
        velocity
    end
    
    
    methods
        
        function obj = MotionBout(start_idx, end_idx, trial)
        % MotionBout
        %
        %   MotionBout(START_SAMPLE, END_SAMPLE, TRIAL)
        %   creates the object by taking sample points START_SAMPLE and
        %   END_SAMPLE in the trial TRIAL (object of class Trial) which
        %   reference a motion period.
        
            obj.start_idx = start_idx;
            obj.end_idx = end_idx;
            
            % useful info
            obj.trial = trial;
        end
        
        
        function val = get.duration(obj)
        %%return the duration of the motion bout
            val = (obj.end_idx - obj.start_idx)/obj.trial.fs;
        end
        
        
        function val = get.start_time(obj)
        %%return the start time of the motion bout in probe time
            val = obj.trial.probe_t(obj.start_idx);
        end
        
        
        function val = get.end_time(obj)
        %%return the end time of the motion bout in probe time
            val = obj.trial.probe_t(obj.end_idx);
        end
        
        
        function val = get.velocity(obj)
        %%return the velocty trace of the motion bout between `prepad`
        %%seconds before onset and `postpad` seconds after offset
        
            start_idx = obj.start_idx - obj.trial.fs*obj.prepad; %#ok<PROP>
            end_idx = obj.end_idx + obj.trial.fs*obj.postpad; %#ok<PROP>
            val = obj.trial.velocity(start_idx:end_idx); %#ok<PROP>
        end
    end
end