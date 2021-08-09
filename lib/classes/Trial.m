classdef Trial < handle
    
    properties
        
        id
        trial
        session
        config
        
        protocol
        start_idx
        end_idx
        
        is_replay = false
        replay_of
        replayed_trial_id
        
        min_bout_duration = 2;
        use_camera_for_stationary = false;
    end
    
    properties (Dependent = true)
        
        fs
        rc2_t
        probe_t
        
        treadmill_speed
        treadmill_acceleration
        stage_speed
        stage_acceleration
        multiplexer_speed
        multiplexer_acceleration
        
        velocity
        acceleration
        
        solenoid
        stage
        gain_teensy
        filtered_teensy
        filtered_teensy_2
        raw_teensy
        multiplexer_output
        pump
        lick
        teensy_gain
        camera0
        camera1
        camera_motion_mask
    end
    
    
    methods
        
        function obj = Trial(trial, session)
        
            obj.trial = trial;
            obj.id = trial.id;
            obj.session = session;
            obj.config = trial.config;
            
            obj.protocol = trial.protocol;
            obj.start_idx = trial.start_idx;
            obj.end_idx = trial.end_idx;
            
            if any(strcmp(obj.protocol, {'StageOnly', 'ReplayOnly'}))
                
                replayed_trial = obj.find_replayed_trial();
                               
                if ~isempty(replayed_trial)
                    if strcmp(obj.protocol, 'ReplayOnly')
                        obj.config.forward_limit = obj.config.start_pos - ...
                            (max(replayed_trial.config.stage_pos, replayed_trial.config.start_pos) - replayed_trial.config.forward_limit);
                    end
                    
                    obj.is_replay = true;
                    obj.replay_of = replayed_trial.protocol;
                    obj.replayed_trial_id = replayed_trial.id;
                else
                    
                    obj.is_replay = true;
                    obj.replay_of = 'Bank';
                    obj.config.forward_limit = 250; %%%% HMMMM, needs to be changed
                    obj.replayed_trial_id = [];
                end
            end
            
            
            % correct for missing data
            obj = correct_at_trial_object_creation(obj);
        end
        
        
        function val = get.fs(obj)
            val = obj.session.fs;
        end
        
        
        
        function val = get.rc2_t(obj)
            
            val = obj.session.rc2_t(obj.start_idx:obj.end_idx);
        end
        
        
        function val = get.probe_t(obj)
            
            val = obj.session.probe_t(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.solenoid(obj)
            
            val = obj.session.solenoid(obj.start_idx:obj.end_idx);
        end
        
        
        function val = get.lick(obj)
            
            val = obj.session.lick(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.pump(obj)
            
            val = obj.session.pump(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.camera0(obj)
            val = obj.session.camera0(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.camera1(obj)
            val = [];
            if ~isempty(obj.session.camera1)
                val = obj.session.camera1(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.camera_motion_mask(obj)
            val = obj.session.camera_motion_mask(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.teensy_gain(obj)
            val = obj.session.teensy_gain(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.filtered_teensy_2(obj)
            val = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.raw_teensy(obj)
            val = obj.session.raw_teensy(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.multiplexer_output(obj)
            val = [];
            if ~isempty(obj.session.multiplexer_output)
                val = obj.session.multiplexer_output(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.gain_teensy(obj)
             val = obj.session.gain_teensy(obj.start_idx:obj.end_idx);
        end
        
        
        function val = get.treadmill_speed(obj)
            
            switch obj.protocol
                case {'EncoderOnly', 'Coupled', 'StageOnly', 'ReplayOnly'}
                    trace = obj.session.filtered_teensy(obj.start_idx:obj.end_idx);
                case {'EncoderOnlyMismatch', 'CoupledMismatch'}
                    trace = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
            end
            
            val = filter_trace(trace);
        end
        
        
        function val = get.treadmill_acceleration(obj)
        
            acc = 100*diff(obj.treadmill_speed); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
            
        end
        
        
        function val = get.stage_speed(obj)
            
            val = filter_trace(obj.stage);
        end
        
        
        function val = get.stage_acceleration(obj)
        
            acc = 100*diff(obj.stage_speed); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
            
        end
        
        
        
        function val = get.multiplexer_speed(obj)
            
            trace = obj.session.multiplexer_output(obj.start_idx:obj.end_idx);
            val = filter_trace(trace);
            
        end
        
        
        
        function val = get.multiplexer_acceleration(obj)
        
            acc = 100 * diff(obj.multiplexer_speed); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
            
        end
        
        
        
        function val = get.velocity(obj)
            
            switch obj.protocol
                case {'EncoderOnly', 'Coupled'}
                    trace = obj.session.filtered_teensy(obj.start_idx:obj.end_idx);
                case {'StageOnly'}
                    trace = obj.session.stage(obj.start_idx:obj.end_idx);
                case {'ReplayOnly'}
                    trace = obj.session.multiplexer_output(obj.start_idx:obj.end_idx);
                case {'EncoderOnlyMismatch', 'CoupledMismatch'}
                    trace = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
            end
            
            val = filter_trace(trace);
        end
        
        
        function val = get.acceleration(obj)
            acc = 100*diff(obj.velocity); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
        end
        
        
        function val = get.stage(obj)
             val = obj.session.stage(obj.start_idx:obj.end_idx);
        end
        
        
        function val = get.filtered_teensy(obj)
            val = [];
            if ~isempty(obj.session.filtered_teensy)
                val = obj.session.filtered_teensy(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function trial = find_replayed_trial(obj)
            c = [obj.session.trials(:).config];
            if isempty(c)
                trial = [];
            else
                idx = strcmp({c(:).log_fname}, obj.config.wave_fname);
                trial = obj.session.trials(idx);
            end
        end
        
        
        function baseline_mask = baseline_mask(obj)
            % original name was not good
            baseline_mask = obj.baseline_window();
            
        end
        
        function idx = baseline_window(obj)
             
            switch obj.protocol
                case {'EncoderOnly', 'Coupled', 'CoupledMismatch', 'EncoderOnlyMismatch'}
                    sol_down_i = find(diff(obj.solenoid > 2.5) == -1, 1);
                    
                    idx = obj.rc2_t > obj.rc2_t(sol_down_i) - 3 & ...
                        obj.rc2_t < obj.rc2_t(sol_down_i) - 1;
                    
                case {'StageOnly', 'ReplayOnly'}
                    sol_up_i = find(diff(obj.solenoid > 2.5) == 1, 1);
                    idx = obj.rc2_t > obj.rc2_t(sol_up_i) + 2 & ...
                        obj.rc2_t < obj.rc2_t(sol_up_i) + 4;
            end
        end
        
        

        function val = position(obj, deadband, mask)
            
            VariableDefault('deadband', 0.1);
            VariableDefault('mask', false(size(obj.velocity)));
            
            % set optional mask positions to zero
            vel = obj.velocity;
            vel(mask) = 0;
            
            % set velocity within limits to zero
            vel(abs(vel) < deadband) = 0;
            
            val = cumsum(vel)/obj.fs;
        end
        
        
        
        function motion_mask = treadmill_motion_mask(obj, vel_thresh, acc_thresh, min_dur)
            
            VariableDefault('vel_thresh', 1);
            VariableDefault('acc_thresh', 0.5);
            VariableDefault('min_dur', 0.2);
            
            % perform thresholding of velocity and acceleration
            stationary_mask = obj.treadmill_speed < vel_thresh & ...
                abs(obj.treadmill_acceleration) < acc_thresh;
            
            cc = bwconncomp(stationary_mask);
            s = cellfun(@(x)(x(1)), cc.PixelIdxList);
            e = cellfun(@(x)(x(end)), cc.PixelIdxList);
            
            dur = (e - s)/obj.fs;
            
            rm_idx = dur < min_dur;
            
            stationary_mask(cat(1, cc.PixelIdxList{rm_idx})) = false;
            motion_mask = ~stationary_mask;
        end
        
        
        
        function motion_mask = stage_motion_mask(obj, vel_thresh, acc_thresh, min_dur)
            
            VariableDefault('vel_thresh', 1);
            VariableDefault('acc_thresh', 0.5);
            VariableDefault('min_dur', 0.2);
            
            % perform thresholding of velocity and acceleration
            stationary_mask = obj.stage_speed < vel_thresh & ...
                abs(obj.stage_acceleration) < acc_thresh;
            
            cc = bwconncomp(stationary_mask);
            s = cellfun(@(x)(x(1)), cc.PixelIdxList);
            e = cellfun(@(x)(x(end)), cc.PixelIdxList);
            
            dur = (e - s)/obj.fs;
            
            rm_idx = dur < min_dur;
            
            stationary_mask(cat(1, cc.PixelIdxList{rm_idx})) = false;
            motion_mask = ~stationary_mask;
            
        end
        
        
        
        function motion_mask = multiplexer_motion_mask(obj, vel_thresh, acc_thresh, min_dur)
            
            VariableDefault('vel_thresh', 1);
            VariableDefault('acc_thresh', 0.5);
            VariableDefault('min_dur', 0.2);
            
            % perform thresholding of velocity and acceleration
            stationary_mask = obj.multiplexer_speed < vel_thresh & ...
                abs(obj.multiplexer_acceleration) < acc_thresh;
            
            cc = bwconncomp(stationary_mask);
            s = cellfun(@(x)(x(1)), cc.PixelIdxList);
            e = cellfun(@(x)(x(end)), cc.PixelIdxList);
            
            dur = (e - s)/obj.fs;
            
            rm_idx = dur < min_dur;
            
            stationary_mask(cat(1, cc.PixelIdxList{rm_idx})) = false;
            motion_mask = ~stationary_mask;
            
        end
        
        
        
        function mask = analysis_window(obj, remove_after_solenoid_start)
            
            VariableDefault('remove_after_solenoid_start', 0.2);
            
            remove_at_end = 0.55;
            add_before_solenoid_low = 2;
            
            if any(strcmp(obj.protocol, {'Coupled', 'EncoderOnly'}))
                
                % where is solenoid low
                idx = find(obj.solenoid < 2.5) + 1; idx(end) = [];
                
                assert(length(unique(diff(idx))) == 1, 'solenoid not down contiguously');
                
                % add in 2s before the solenoid goes low, 
                % remove 200ms after solenoid goes low
                % remove 500ms before solenoid goes high
                idx_new = (idx(1) - add_before_solenoid_low * obj.fs) : (idx(1)-1);
                idx_new = [idx_new, (idx(1) + remove_after_solenoid_start * obj.fs) : (idx(end) - remove_at_end * obj.fs)];
            
            elseif any(strcmp(obj.protocol, {'CoupledMismatch', 'EncoderOnlyMismatch'}))
                
                after_mm_to_remove = 3;
                
                % where is solenoid low
                idx = find(obj.solenoid < 2.5) + 1; idx(end) = [];
                
                assert(length(unique(diff(idx))) == 1, 'solenoid not down contiguously');
                
                % add in 2s before the solenoid goes low, 
                % remove 200ms after solenoid goes low
                % remove 500ms before solenoid goes high
                idx_new = (idx(1) - add_before_solenoid_low * obj.fs) : (idx(1)-1);
                idx_new = [idx_new, (idx(1) + remove_after_solenoid_start * obj.fs) : (idx(end) - remove_at_end * obj.fs)];
                
                % find mismatch onset
                mm_onset_idx = find(diff(obj.teensy_gain > 2.5) == 1) + 1;
                
                % remove 3s after mismatch onset
                idx_new = [idx_new(1):mm_onset_idx, (mm_onset_idx + after_mm_to_remove*obj.fs):idx_new(end)];
                
            elseif any(strcmp(obj.protocol, {'StageOnly', 'ReplayOnly'}))
                
                % find where the solenoid goes high (mouse initiated trial)
                start_idx = find(obj.solenoid > 2.5, 1); %#ok<PROPLC>
                
                % remove first 200ms after solenoid high
                start_idx = start_idx + remove_after_solenoid_start * obj.fs; %#ok<PROPLC>
                
                % where trace reaches final position
                end_idx = find(obj.position > 0.98*(obj.config.start_pos - obj.config.forward_limit)/10, 1); %#ok<PROPLC>
                
                % full index
                idx_new = start_idx : end_idx; %#ok<PROPLC>
            end
            
            % create the mask
            mask = false(length(obj.solenoid), 1);
            mask(idx_new) = true;
            
            % include the baseline period in the analysis window
            mask = mask | obj.baseline_mask();
        end
        
        
        
        function mask = stationary_mask(obj)
            
            % in analysis window, treadmill not in motion, camera not
            % in motion, stage not in motion
            mask = obj.analysis_window() & ...
                   ~ obj.treadmill_motion_mask() & ...
                   ~ obj.stage_motion_mask;
               % ~ obj.camera_motion_mask & ...
            
            if obj.use_camera_for_stationary
                mask = mask & ~ obj.stage_motion_mask;
            end
               
               
            if strcmp(obj.protocol, 'ReplayOnly')
                mask = mask & ~ obj.multiplexer_motion_mask();
            end
        end
        
        
        
        function mask = motion_mask(obj, remove_after_solenoid_start)
            
            VariableDefault('remove_after_solenoid_start', 0.2);
            
            switch obj.protocol
                case {'Coupled', 'EncoderOnly', 'CoupledMismatch', 'EncoderOnlyMismatch'}
                    
                    mask = obj.analysis_window(remove_after_solenoid_start) & ...
                       obj.treadmill_motion_mask() & ...
                       obj.solenoid < 2.5;
                    
                case 'StageOnly'
                   
                    mask = obj.analysis_window(remove_after_solenoid_start) & ...
                           obj.stage_motion_mask();
                       
                case 'ReplayOnly'
                    
                    mask = obj.analysis_window(remove_after_solenoid_start) & ...
                           obj.multiplexer_motion_mask();
            end
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
        %% include_200ms - whether to include the first 200ms of the trial after the solenoid goes low
        %       by default this is false
        
            VariableDefault('include_200ms', false);
            VariableDefault('true_motion_start', true);
            
            % split the motion mask into running bouts
            if include_200ms
                cc = bwconncomp(obj.motion_mask(0));
            else
                cc = bwconncomp(obj.motion_mask());
            end
            s = cellfun(@(x)(x(1)), cc.PixelIdxList);
            e = cellfun(@(x)(x(end)), cc.PixelIdxList);
            
%             if true_motion_start
%                 % make sure the sample before is genuinely stationary and not
%                 % just the start of an analysis window
%                 stat_mask = obj.stationary_mask();
%                 true_motion_start = stat_mask(s-1);
%                 
%                 s = s(true_motion_start);
%                 e = e(true_motion_start);
%             end
            
            if true_motion_start
                
                if include_200ms
                    analysis_window = obj.analysis_window(0);
                else
                    analysis_window = obj.analysis_window();
                end
                
                starts_in_analysis_window = analysis_window(s-1);
                s = s(starts_in_analysis_window);
                e = e(starts_in_analysis_window);
            end

            
            if isempty(s)
                bouts = [];
                return
            end
            
            for i = 1 : length(s)
                bouts(i) = MotionBout(s(i), e(i), obj);
            end
            
%             bouts = bouts([bouts(:).duration] > obj.min_bout_duration);
        end
        
        
        
        function t = mismatch_onset_t(obj)
            
            mm_onset_idx = find(diff(obj.teensy_gain > 2.5) == 1) + 1;
            t = obj.probe_t(mm_onset_idx);
        end
        
        
        
        function t = mismatch_offset_t(obj)
            
            mm_offset_idx = find(diff(obj.teensy_gain > 2.5) == -1) + 1;
            t = obj.probe_t(mm_offset_idx);
        end
        
        
        
        function mask = mismatch_window(obj)
            
            mm_onset_idx = find(diff(obj.teensy_gain > 2.5) == 1) + 1;
            mm_offset_idx = find(diff(obj.teensy_gain > 2.5) == -1) + 1;
            
            mask = false(length(obj.rc2_t), 1);
            mask(mm_onset_idx : (mm_offset_idx + 0.05*obj.fs)) = true;
        end
    end
end
