classdef Trial < handle
% Trial Class for handling details of a trial (1 lap of the linear track,
% linear corridor etc.)
%
%   Trial Properties:
%       probe_id            - string with the probe recording ID the trial was part of
%       session_id          - string with the session ID the trial was part of
%       trial_id            - integer with the trial ID
%       trial_group_id      - integer ID of the trial group, see RVTSession
%       trial_group_label   - string ID of the trial group, see RVTSession
%       trial               - the original structure from the formatted data file
%       session             - handle to the parent Session object
%       config              - structure containing the configuration for
%                             this trial (loaded from .cfg file associated with session .bin file)
%       protocol            - 'Coupled', 'EncoderOnly', 'ReplayOnly', 'StageOnly', 'CoupledMismatch', 'EncoderOnlyMismatch'
%       start_idx           - sample point in the session on which the trial starts
%       end_idx             - sample point in the session on which the trial ends
%       is_replay           - true if trial is a replay of another trial
%       replay_of           - if trial is a replay, this is a string containing the protocol
%                             name of the original trial.
%       original_trial_id   - integer ID of original trial
%       replay_offset       - if the trial is a replay, the offset applied
%                             to the replay trial in order to match it to the original
%       use_camera_for_stationary - whether to include the motion energy from the camera in computing the stationary periods
%       
%       fs                  - sample rate of the NIDAQ session recording
%       rc2_t               - timebase of the trial in the NIDAQ session recording (0s on first sample point, 1/fs on second etc.)
%       probe_t             - timebase of the trial synced to the probe recording
%       camera_t            - timebase of the camera frames synced to the probe recording
%       camera_idx          - index of the session camera frames appearing in this trial
%       treadmill_speed     - speed of the treadmill
%       treadmill_acceleration - acceleration of the treadmill
%       stage_speed         - speed of the stage
%       stage_acceleration  - acceleration of the stage
%       multiplexer_speed   - speed command coming from the multiplexer
%       multiplexer_acceleration - acceleration of multiplexer_speed
%       visual_speed        - command speed to the visual corridor during this trial
%       velocity            - speed for a particular trial (will select
%                             from `treadmill_speed` or `stage_speed` as appropriate)
%       acceleration        - acceleration for a particular trial (will select
%                             from `treadmill_speed` or `stage_speed` as appropriate)
%       solenoid            - solenoid channel
%       stage               - speed of the stage (c.f. stage_speed)
%       gain_teensy         - gain_teensy channel
%       filtered_teensy     - filtered_teensy channel
%       filtered_teensy_2   - filtered_teensy_2 channel
%       raw_teensy          - raw_teensy channel
%       multiplexer_output  - multiplexer_output channel
%       pump                - pump channel
%       lick                - lick channel
%       teensy_gain         - teensy_gain channel
%       camera0             - motion energy for camera0
%       camera1             - motion energy for camera1
%
%   Trial Methods:
%       find_original_trial - for replay trials, finds the original trial in the set of all trials
%       baseline_mask       - same as `baseline_window`
%       baseline_window     - get mask defining a 'baseline' period where the animal should not be in motion
%       position            - integrate `velocity` trace to give estimate of position along track
%       treadmill_motion_mask - mask defining where periods of treadmill is in motion
%       stage_motion_mask   - mask defining where periods of linear stage is in motion
%       multiplexer_motion_mask - mnask defining where periods of multiplexer command is in motion
%       analysis_window     - mask defining periods where the trial is to be analyzed
%       stationary_mask     - mask defining the stationary periods for the trial
%       motion_mask         - mask defining the motion periods for the trial
%       stationary_time     - amount of time the trial is stationary
%       motion_time         - amount of time the trial is in motion
%       analysis_window_time - amount of time in the analysis window for the trial
%       motion_bouts        - calculate and return motion bouts in the trial
%       mismatch_onset_t    - time (in probe recording time) on which a mismatch event starts
%       mismatch_offset_t   - time (in probe recording time) on which a mismatch event begins to ramp down
%       mismatch_onset_sample - sample point (of the trial) on which a mismatch event starts
%       mismatch_offset_sample - sample point (of the trial) on which a mismatch event end (gain is ramped_down)
%       mismatch_saliency   - return the saliency of the mismatch trial
%       is_mismatch_trial   - boolean, whether this trial is a "mismatch" trial
%       mismatch_duration   - duration of the mismatch between starting to ramp up and starting to ramp down
%       to_aligned          - create an object of type AlignedTrial for this trial
%       
%
%   See also: RVTSession, AlignedTrial
%
%   TODO:   1. this class is too big and needs to be refactored
%           2. properties for each channel should become a Map object?
%           3. masks are too complicated
%           4. make all replay trials aligned by default

    properties
        
        probe_id
        session_id
        trial_id
        
        trial_group_id
        trial_group_label
        trial
        session
        config
        
        protocol
        start_idx
        end_idx
        
        is_replay = false
        replay_of
        original_trial_id
        replay_offset = 0  % samples
        
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
        visual_speed
        
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
    end
    
    properties (Dependent = true, Hidden = true)
        
        camera0_
        camera1_
        camera_t
        camera_idx
    end
    
    
    methods
        
        function obj = Trial(trial, session)
        % Trial
        %
        %   Trial(TRIAL, SESSION) creates an object handling details of a
        %   trial. Takes as arguments TRIAL, which is a MATLAB structure
        %   created on the formatting step and saved to the formatted data
        %   structure, and SESSION an object of class RVTSession or one of
        %   its subclasses.
        %
        %   Currently, this is a monolithic class for a trial
        %   this needs to be split up
        
            obj.trial = trial;
            obj.trial_id = trial.trial_id;
            obj.probe_id = trial.probe_id;
            obj.session_id = trial.session_id;
            obj.session = session;
            obj.config = trial.config;
            
            obj.protocol = trial.protocol;
            obj.start_idx = trial.start_idx;
            obj.end_idx = trial.end_idx;
            
            if any(strcmp(obj.protocol, {'StageOnly', 'ReplayOnly'}))
                
                original_trial = obj.find_original_trial();
                               
                if ~isempty(original_trial)
                    if strcmp(obj.protocol, 'ReplayOnly')
                        obj.config.forward_limit = obj.config.start_pos - ...
                            (max(original_trial.config.stage_pos, original_trial.config.start_pos) - original_trial.config.forward_limit);
                    end
                    
                    obj.is_replay = true;
                    obj.replay_of = original_trial.protocol;
                    obj.original_trial_id = original_trial.trial_id;
                else
                    
                    obj.is_replay = true;
                    obj.replay_of = 'Bank';
                    obj.config.forward_limit = 250; %%%% HMMMM, needs to be changed
                    obj.original_trial_id = [];
                end
            end
            
            
            % correct for missing data
            obj = correct_at_trial_object_creation(obj);
        end
        
        
        function val = get.fs(obj)
        %%sampling frequency of the NIDAQ data in the trial
            val = obj.session.fs;
        end
        
        
        
        function val = get.rc2_t(obj)
        %%timebase of the trial in the session
            val = obj.session.rc2_t(obj.start_idx:obj.end_idx);
        end
        
        
        function val = get.probe_t(obj)
        %%timebase of the trial synced to the probe recording
            val = obj.session.probe_t(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.camera_idx(obj)
        %%index of camera frames in the session associated with this trial
            start_t = obj.session.probe_t(obj.start_idx);
            end_t = obj.session.probe_t(obj.end_idx);
            val = obj.session.camera_t >= start_t & obj.session.camera_t <= end_t;
        end
        
        
        function val = get.camera_t(obj)
        %%timebase of the camera frames in this trial synced to the probe recording    
            val = obj.session.camera_t(obj.camera_idx);
        end
        
        
        
        function val = get.solenoid(obj)
        %%solenoid channel during this trial    
            val = obj.session.solenoid(obj.start_idx:obj.end_idx);
        end
        
        
        function val = get.lick(obj)
        %%lick channel during this trial    
            val = obj.session.lick(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.pump(obj)
        %%pump channel during this trial    
            val = obj.session.pump(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.camera0_(obj)
        %%motion energy for camera0, camera timebase
            val = [];
            if ~isempty(obj.session.camera0_)
                val = obj.session.camera0_(obj.camera_idx);
            end
        end
        
        
        
        function val = get.camera1_(obj)
        %%motion energy for camera1, camera timebase
            val = [];
            if ~isempty(obj.session.camera1_)
               val = obj.session.camera1_(obj.camera_idx);
            end
        end
        
        
        function val = get.camera0(obj)
        %%motion energy for camera0, interpolated to session
            val = [];
            if ~isempty(obj.session.camera0)
                val = obj.session.camera0(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.camera1(obj)
        %%motion energy for camera1, interpolated to session
            val = [];
            if ~isempty(obj.session.camera1)
               val = obj.session.camera1(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.teensy_gain(obj)
        %%teensy_gain channel during this trial  
            val = [];
            if ~isempty(obj.session.teensy_gain)
                val = obj.session.teensy_gain(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.filtered_teensy_2(obj)
        %%filtered_teensy_2 channel during this trial
            val = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.raw_teensy(obj)
        %%raw_teensy channel during this trial
            val = obj.session.raw_teensy(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.multiplexer_output(obj)
        %%multiplexer_output channel during this trial
            val = [];
            if ~isempty(obj.session.multiplexer_output)
                val = obj.session.multiplexer_output(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.gain_teensy(obj)
        %%gain_teensy channel during this trial
             val = obj.session.gain_teensy(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.treadmill_speed(obj)
        %%speed of the treadmill during this trial
        
            switch obj.protocol
                case {'EncoderOnly', 'Coupled', 'StageOnly', 'ReplayOnly'}
                    if ~isempty(obj.session.filtered_teensy)
                        trace = obj.session.filtered_teensy(obj.start_idx:obj.end_idx);
                    else
                        trace = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
                    end
                case {'EncoderOnlyMismatch', 'CoupledMismatch'}
                    trace = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
            end
            
            val = filter_trace(trace);
        end
        
        
        
        function val = get.treadmill_acceleration(obj)
        %%acceleration of the treadmill during this trial
            acc = 100*diff(obj.treadmill_speed); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
            
        end
        
        
        function val = get.stage_speed(obj)
        %%speed of the stage during this trial
        
            val = filter_trace(obj.stage);
        end
        
        
        function val = get.stage_acceleration(obj)
        %%acceleration of the stage during this trial
        
            acc = 100*diff(obj.stage_speed); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
            
        end
        
        
        
        function val = get.visual_speed(obj)
        %%command speed to the visual corridor during this trial
        %   = multiplexer_speed if visual stimulus is enabled
        %   = 0 if visual stimulus is not enabled
        
            if obj.config.enable_vis_stim
                val = obj.multiplexer_speed;
            else
                val = zeros(size(obj.multiplexer_speed));
            end
        end
        
        
        
        function val = get.multiplexer_speed(obj)
        %%speed command coming from the multiplexer
        
            trace = obj.session.multiplexer_output(obj.start_idx:obj.end_idx);
            val = filter_trace(trace);  
        end
        
        
        
        function val = get.multiplexer_acceleration(obj)
        %%acceleration of multiplexer_speed
            acc = 100 * diff(obj.multiplexer_speed); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
            
        end
        
        
        
        function val = get.velocity(obj)
        %%speed for a particular trial (will select from `treadmill_speed` or `stage_speed` as appropriate)
        
            switch obj.protocol
                case {'EncoderOnly', 'Coupled'}
                    if ~isempty(obj.session.filtered_teensy)
                        trace = obj.session.filtered_teensy(obj.start_idx:obj.end_idx);
                    else
                        trace = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
                    end
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
        %%acceleration for a particular trial (will select from `treadmill_speed` or `stage_speed` as appropriate)
        
            acc = 100*diff(obj.velocity); % cm/s to m/s^2
            val = [acc(1); (acc(1:end-1)+acc(2:end))/2; acc(end)];
        end
        
        
        function val = get.stage(obj)
        %%unfiltered speed of the stage (c.f. stage_speed)
             val = obj.session.stage(obj.start_idx:obj.end_idx);
        end
        
        
        
        function val = get.filtered_teensy(obj)
        %%filtered_teensy channel during this trial
            val = [];
            if ~isempty(obj.session.filtered_teensy)
                val = obj.session.filtered_teensy(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function trial = find_original_trial(obj)
        %%find_original_trial For replay trials, finds the original trial in the set of all trials
        %
        %   TRIAL = find_original_trial()
        %   returns the Trial object associated with the original trial
        
            trial_array = [obj.session.trials{:}];
            config_array = [trial_array(:).config];
            if isempty(config_array)
                trial = [];
            else
                idx = strcmp({config_array(:).log_fname}, obj.config.wave_fname);
                if sum(idx) == 1
                    trial = obj.session.trials{idx};
                else
                    trial = [];
                end
            end
        end
        
        
        function baseline_mask = baseline_mask(obj)
        %renames `baseline_window`
            baseline_mask = obj.baseline_window();
        end
        
        
        
        function idx = baseline_window(obj)
        %%baseline_window Get mask defining a 'baseline' period where the
        %%animal should not be in motion
        
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
        %%position Integrate `velocity` trace to give estimate of position along track
        %
        %   POSITION = position(DEADBANK, MASK)
        %   integrates the `velocity` trace to give a POSITION trace. If
        %   DEADBAND is supplied it determines the velocity below which no
        %   movement is included (i.e. velocity set to 0, rather than
        %   integrating a residual noise which might have a slight offset).
        %   Default = 0.1cm/s.
        %
        %   If MASK is supplied it should be a boolean vector of length
        %   equal to `velocity` defining any other periods of the velocity
        %   trace which should be set to 0 and therefore not included. By
        %   default it is false everywhere.
        
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
        %%treadmill_motion_mask Mask defining where periods of treadmill is in motion
        %
        %   MASK = treadmill_motion_mask(VELOCITY_THRESHOLD, ACCELERATION_THRESHOLD, MIN_STATIONARY_DURATION)
        %   returns a boolean mask, MASK, of length equal to `velocity`
        %   defining periods of the trial where the treadmill is in
        %   motion. First, stationary periods are defined as those where
        %   the velocity < VELOCITY_THRESHOLD, the acceleration is <
        %   ACCELERATION_THRESHOLD and the period lasts for >
        %   MIN_STATIONARY_DURATION. All are optional and defaults are:
        %       VELOCITY_THRESHOLD = 1 cm/s
        %       ACCELERATION_THRESHOLD = 0.5 m/s^2
        %       MIN_STATIONARY_DURATION = 0.2 s
        %
        %   The motion mask returned is then the logical negation of this stationary
        %   mask.
        %
        %   Note that no `analysis_window` is applied in this method so
        %   mask applies to the entire trial.
        
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
        %%stage_motion_mask Mask defining where periods of linear stage is in motion
        %
        %   MASK = stage_motion_mask(VELOCITY_THRESHOLD, ACCELERATION_THRESHOLD, MIN_STATIONARY_DURATION)
        %   returns a boolean mask, MASK, of length equal to `velocity`
        %   defining periods of the trial where the linear stage is in
        %   motion. First, stationary periods are defined as those where
        %   the velocity < VELOCITY_THRESHOLD, the acceleration is <
        %   ACCELERATION_THRESHOLD and the period lasts for >
        %   MIN_STATIONARY_DURATION. All are optional and defaults are:
        %       VELOCITY_THRESHOLD = 1 cm/s
        %       ACCELERATION_THRESHOLD = 0.5 m/s^2
        %       MIN_STATIONARY_DURATION = 0.2 s
        %
        %   The motion mask returned is then the logical negation of this stationary
        %   mask.
        %
        %   Note that no `analysis_window` is applied in this method so
        %   mask applies to the entire trial.
        
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
        %%multiplexer_motion_mask Mask defining where periods of
        %%multiplexer command is in motion
        %
        %   MASK = multiplexer_motion_mask(VELOCITY_THRESHOLD, ACCELERATION_THRESHOLD, MIN_STATIONARY_DURATION)
        %   returns a boolean mask, MASK, of length equal to `velocity`
        %   defining periods of the trial where the multiplexer output is in
        %   motion. First, stationary periods are defined as those where
        %   the velocity < VELOCITY_THRESHOLD, the acceleration is <
        %   ACCELERATION_THRESHOLD and the period lasts for >
        %   MIN_STATIONARY_DURATION. All are optional and defaults are:
        %       VELOCITY_THRESHOLD = 1 cm/s
        %       ACCELERATION_THRESHOLD = 0.5 m/s^2
        %       MIN_STATIONARY_DURATION = 0.2 s
        %
        %   The motion mask returned is then the logical negation of this stationary
        %   mask.
        %
        %   Note that no `analysis_window` is applied in this method so
        %   mask applies to the entire trial.
        
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
        %%analysis_window Mask defining periods where the trial is to be analyzed
        %
        %   MASK = analysis_window(REMOVE_AFTER_SOLENOID_START)
        %   creates a boolean mask, MASK, of length equal to `velocity`
        %   defining periods of the trial where traces are to be analyzed.
        %
        %   If REMOVE_AFTER_SOLENOID_START is supplied, it defines a
        %   period, in s, after the solenoid trace goes low to remove from
        %   the analysis window. If it is not supplied or empty, the
        %   default is set to 0.2 s. (This is removed because there is a
        %   200ms period where the gain is ramped up from 0 to 1 to avoid
        %   very large accelerations when the solenoid goes low.
        %
        %   The analysis window depends on the trial type.
        
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
        %%stationary_mask Mask defining the stationary periods
        %
        %   MASK = stationary_mask()
        %   returns a mask, MASK, defining periods of the trial in the
        %   analysis window, which are defined as being stationary (no
        %   motion of any kind).
        
            % in analysis window, treadmill not in motion, camera not
            % in motion, stage not in motion
            mask = obj.analysis_window() & ...
                   ~ obj.treadmill_motion_mask() & ...
                   ~ obj.stage_motion_mask;
               
            if strcmp(obj.protocol, 'ReplayOnly')
                mask = mask & ~ obj.multiplexer_motion_mask();
            end
        end
        
        
        
        function mask = motion_mask(obj, remove_after_solenoid_start)
         %%stationary_mask Mask defining the motion periods for the trial
        %
        %   MASK = motion_mask(REMOVE_AFTER_SOLENOID_START)
        %   returns a mask, MASK, defining periods of the trial in the
        %   analysis window, which are defined as being in motion (no
        %   motion of any kind). If REMOVE_AFTER_SOLENOID_START is
        %   supplied, it defines the amount, in s, of the analysis window to
        %   remove after the solenoid goes low. By default, it is 0.2s.
        %
        %   See also: analysis_window
        
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
        %%stationary_time Amount of time the trial is stationary
        %
        %   TIME = stationary_time()
        
            val = sum(obj.stationary_mask()) / obj.fs;
        end
        
        
        
        function val = motion_time(obj)
        %%motion_time Amount of time the trial is in motion
        %
        %   TIME = motion_time()
            
            val = sum(obj.motion_mask()) / obj.fs;
        end
        
        
        
        function val = analysis_window_time(obj)
        %%analysis_window_time Amount of time in the analysis window for
        %%the trial
        %
        %   TIME = analysis_window_time()
        
            val = sum(obj.analysis_window()) / obj.fs;
        end
        
        
        
        function bouts = motion_bouts(obj, include_200ms, true_motion_start)
        %%motion_bouts Calculate and return motion bouts in the trial
        %
        %   BOUTS_STRUCT = motion_bouts(INCLUDE_200MS, TRUE_MOTION_START)
        %   calculates and returns the motion bouts occuring in the trial.
        %   BOUTS_STRUCT is a structure array containing information about
        %   each motion bout in the trial. See MotionBout for more details.
        %
        %   If supplied, INCLUDE_200MS should be a boolean value declaring
        %   whether to include the first 200ms after the solenoid goes low
        %   to look for motion bouts. If true, the 200ms is included and
        %   any motion straight after the solenoid goes low is included in
        %   the first bout. If not supplied or empty, INCLUDE_200MS is set
        %   to false.
        %
        %   TRUE_MOTION_START should not be an option, but remains for
        %   legacy reasons. Since the analysis window in a mismatch trial
        %   begins again 3s after the mismatch onset, if there is motion
        %   at the re-start of this analysis window, it used to be deemed
        %   as the start of a bout. However, it is not necessarily the case
        %   that it is the start of a motion bout. Therefore, this option
        %   set to true ensures that it really is a motion bout by looking
        %   at the sample point before the motion is first detected. By
        %   default it is true, and it should always be set to true.
        %
        % See also: MotionBout
        
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
                bouts{i} = MotionBout(s(i), e(i), obj);
            end
            
%             bouts = bouts([bouts(:).duration] > obj.min_bout_duration);
        end
        
        
        
        function t = mismatch_onset_t(obj)
        %%mismatch_onset_t Time (in probe recording time) on which a
        %%mismatch event starts
        %
        %   TIME = mismatch_onset_t()
        
            idx = obj.mismatch_onset_sample();
            t = obj.probe_t(idx);
        end
        
        
        
        function t = mismatch_offset_t(obj)
        %%mismatch_offset_t Time (in probe recording time) on which a
        %%mismatch event begins to ramp down
        %
        %   TIME = mismatch_offset_t()
        
            idx = obj.mismatch_offset_sample();
            t = obj.probe_t(idx);
        end
        
        
        
        function idx = mismatch_onset_sample(obj)
        %%mismatch_onset_sample Sample point (of the trial) on which a
        %%mismatch event starts
        %
        %   SAMPLE = mismatch_onset_sample()
        
            idx = find(diff(obj.teensy_gain > 2.5) == 1) + 1;
        end
        
        
        
        function idx = mismatch_offset_sample(obj)
        %%mismatch_offset_sample Sample point (of the trial) on which a
        %%mismatch event end (begins ramp down of gain)
        %
        %   SAMPLE = mismatch_offset_sample()
        
            idx = find(diff(obj.teensy_gain > 2.5) == -1) + 1;
        end
        
        
        
        function [delta_speed, idx] = mismatch_saliency(obj)
        %%mismatch_saliency Saliency of the mismatch trial
        %
        %   [DELTA_SPEED, IDX] = mismatch_saliency() returns the maximum of
        %   the difference between the running speed and the command speed
        %   to the stage or visual stimulus during a mismatch event.
        %
        %   DELTA_SPEED is in cm/s and IDX is the sample point of the trial
        %   at which the maximum difference occurs.
        
            % if not a mismatch trial return NaN
            if ~obj.is_mismatch_trial()
                delta_speed = nan;
                return
            end
            
            onset_sample = obj.mismatch_onset_sample();
            offset_sample = obj.mismatch_offset_sample();
            
            speed_difference = obj.gain_teensy - obj.filtered_teensy_2;
            
            [delta_speed, idx] = max(speed_difference(onset_sample:offset_sample));
            idx = idx + onset_sample - 1;
        end
        
        
        
        function val = is_mismatch_trial(obj)
        %%is_mismatch_trial Boolean, whether this trial is a "mismatch"
        %%trial
        %
        %   BOOL = is_mismatch_trial()
        
            val = any(strcmp(obj.protocol, {'CoupledMismatch', 'EncoderOnlyMismatch'}));
        end
        
        
        
        function val = mismatch_duration(obj)
        %%is_mismatch_trial Duration of the mismatch between starting to ramp up and starting
        %%to ramp down
        %
        %   DURATION = mismatch_duration()
        %
        %   Note: this is not the full duration of mismatch. It is just
        %   when the command is sent to the Teensy to ramp up or ramp down
        %   the gain.
        
            val = obj.mismatch_offset_t - obj.mismatch_onset_t;
        end
        
        
        function aligned_trial = to_aligned(obj)
        %%to_aligned Create an object of type AlignedTrial for this trial
        %
        %   ALIGNED_TRIAL = to_aligned()
        %   gets the trial to which this trial should be aligned... makes
        %   call to Session and back, which is a bit ugly but makes life
        %   easy
        
            aligned_trial = obj.session.to_aligned(obj);
        end
    end
end
