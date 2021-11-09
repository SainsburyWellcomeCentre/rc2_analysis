classdef AlignedTrial < handle
    
    properties (SetAccess = private)
        
        probe_id
        session_id
        trial_id
        trial_group_label
        original_trial_id
        
        protocol
        fs
        
        trial
        original_trial
        
        offset
        n_points
        is_replay
        replay_of
        replayed_trial_id
    end
    
    
    
    methods
        
        function obj = AlignedTrial(trial, original_trial)
            
            obj.trial = trial;
            obj.original_trial = original_trial;
            
            obj.fs = obj.trial.fs;
            obj.offset = obj.trial.replay_offset;
            obj.n_points = min(length(trial.probe_t) - obj.offset, length(original_trial.probe_t));
        end
        
        
        
        function val = get.probe_id(obj)
            
            val = obj.trial.probe_id;
        end
        
        
        
        function val = get.session_id(obj)
            
            val = obj.trial.session_id;
        end
        
        
        
        function val = get.trial_id(obj)
            
            val = obj.trial.trial_id;
        end
        
        
        
        function val = get.original_trial_id(obj)
            
            val = obj.trial.original_trial_id;
        end
        
        
        
        function val = get.trial_group_label(obj)
            
            val = obj.trial.trial_group_label;
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
            
            val = obj.trial.probe_t(obj.offset + (1:obj.n_points));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = velocity(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.velocity(obj.offset + (1:obj.n_points));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = treadmill_speed(obj, idx)
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.treadmill_speed)
                val = obj.trial.treadmill_speed(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = visual_speed(obj, idx)
            
            VariableDefault('idx', []);
            val = obj.get_aligned_data('visual_speed', idx);
        end
        
        
        
        function val = filtered_teensy(obj, idx)
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.filtered_teensy)
                val = obj.trial.filtered_teensy(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = filtered_teensy_2(obj, idx)
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.filtered_teensy_2)
                val = obj.trial.filtered_teensy_2(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = stage(obj, idx)
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.stage)
                val = obj.trial.stage(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = multiplexer_output(obj, idx)
            
            VariableDefault('idx', []);
            val = obj.get_aligned_data('multiplexer_output', idx);
        end
        
        
        
        function val = teensy_gain(obj, idx)
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.teensy_gain)
                val = obj.trial.teensy_gain(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = camera1(obj, idx)
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.camera1)
                val = obj.trial.camera1(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = solenoid(obj, idx)
            
            VariableDefault('idx', []);
            
            val = obj.trial.solenoid(obj.offset + (1:obj.n_points));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function mask = analysis_window(obj)
        
            mask = obj.original_trial.analysis_window();
            mask = mask(1:obj.n_points);
        end
        
        
        
        function mask = treadmill_motion_mask(obj, vel_thresh, acc_thresh, min_dur)
            
            VariableDefault('vel_thresh', 1);
            VariableDefault('acc_thresh', 0.5);
            VariableDefault('min_dur', 0.2);
        
            mask = obj.trial.treadmill_motion_mask(vel_thresh, acc_thresh, min_dur);
            mask = mask(obj.offset + (1:obj.n_points));
        end
        
        
        function mask = stationary_mask(obj)
            
            mask = obj.original_trial.stationary_mask;
            mask = mask(1:obj.n_points);
%             mask = mask(obj.offset + (1:obj.n_points));
%             mask = mask & obj.analysis_window();
        end
        
        
        
        function mask = motion_mask(obj)
            
            mask = obj.original_trial.motion_mask();
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
        
        
        
        function bouts = motion_bouts(obj, include_200ms, true_motion_start)
            
            VariableDefault('include_200ms', false);
            VariableDefault('true_motion_start', true);
            
            bouts_ori = obj.original_trial.motion_bouts(include_200ms, true_motion_start);
            
            bouts = {};
            for ii = 1 : length(bouts_ori)
                s = bouts_ori{ii}.start_idx;
                e = bouts_ori{ii}.end_idx;
                bouts{ii} = MotionBout(s, e, obj);
            end
        end
        
        
        
        function val = get_aligned_data(obj, type, idx)
            
            val = [];
            if ~isempty(obj.trial.(type))
                val = obj.trial.(type)(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            end
        end
    end
end
