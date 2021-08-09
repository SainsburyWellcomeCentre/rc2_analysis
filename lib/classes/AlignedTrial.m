classdef AlignedTrial < handle
    
    properties
        
        id
        protocol
        fs
        trial
        trial_to_align_to
        offset
        n_points
        is_replay
        replay_of
        replayed_trial_id
        min_bout_duration = 2
    end
    
    
    methods
        
        function obj = AlignedTrial(trial, trial_to_align_to, offset)
            
            obj.trial = trial;
            obj.trial_to_align_to = trial_to_align_to;
            
            obj.fs = obj.trial.fs;
            obj.offset = offset;
            obj.n_points = min(length(trial.probe_t) - offset + 1, ...
                               length(trial_to_align_to.probe_t));
            
        end
        
        
        
        function val = get.id(obj)
            
            val = obj.trial.id;
        end
        
        
        
        function val = get.protocol(obj)
            
            val = obj.trial.protocol;
        end
        
        
        
        function val = get.is_replay(obj)
            
            val = obj.trial.is_replay;
        end
        
        
        
        function val = get.replay_of(obj)
            
            val = obj.trial.replay_of;
        end
        
        
        
        function val = get.replayed_trial_id(obj)
            
            val = obj.trial.replayed_trial_id;
        end
        
        
        
        function val = config(obj)
            
            val = obj.trial.config;
        
        end
        
        
        
        function val = probe_t(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.probe_t(obj.offset + (0:obj.n_points-1));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = velocity(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.velocity(obj.offset + (0:obj.n_points-1));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = filtered_teensy(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.filtered_teensy(obj.offset + (0:obj.n_points-1));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = stage(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.stage(obj.offset + (0:obj.n_points-1));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = multiplexer_output(obj, idx)
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.multiplexer_output)
                val = obj.trial.multiplexer_output(obj.offset + (0:obj.n_points-1));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
                return
            end
        end
        
        
        
        function val = camera1(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.camera1(obj.offset + (0:obj.n_points-1));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = solenoid(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.solenoid(obj.offset + (0:obj.n_points-1));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function mask = analysis_window(obj)
        
            mask = obj.trial_to_align_to.analysis_window();
            mask = mask(1:obj.n_points);
        end
        
        
        
        function mask = treadmill_motion_mask(obj, vel_thresh, acc_thresh, min_dur)
            
            VariableDefault('vel_thresh', 1);
            VariableDefault('acc_thresh', 0.5);
            VariableDefault('min_dur', 0.2);
        
            mask = obj.trial.treadmill_motion_mask(vel_thresh, acc_thresh, min_dur);
            mask = mask(obj.offset + (0:obj.n_points-1));
        end
        
        
        function mask = stationary_mask(obj)
            
            mask = obj.trial_to_align_to.stationary_mask;
            mask = mask(1:obj.n_points);
%             mask = mask(obj.offset + (0:obj.n_points-1));
%             mask = mask & obj.analysis_window();
        end
        
        
        
        function mask = motion_mask(obj)
            
            mask = obj.trial_to_align_to.motion_mask();
            mask = mask(1:obj.n_points);
%             mask = mask & ~obj.stationary_mask();
        end
        
        
        
        function val = stationary_time(obj)
            
            val = sum(obj.stationary_mask()) / obj.fs;
        end
        
        
        
        function val = motion_time(obj)
            
            val = sum(obj.motion_mask()) / obj.fs;
        end
        
        
        
        function val = analysis_window_time(obj)
            
            val = sum(obj.analysis_window()) / obj.fs;
        end
        
        
        
        function bouts = motion_bouts(obj, include_200ms)
            
            VariableDefault('include_200ms', false);
            
            bouts_ori = obj.trial_to_align_to.motion_bouts(include_200ms);
            
            if isempty(bouts_ori)
                bouts = [];
                return
            end
            
            for i = 1 : length(bouts_ori)
                s = bouts_ori(i).start_idx;
                e = bouts_ori(i).end_idx;
                bouts(i) = MotionBout(s, e, obj);
            end
        end
    end
end
