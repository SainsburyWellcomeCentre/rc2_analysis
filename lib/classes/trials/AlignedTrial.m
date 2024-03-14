classdef AlignedTrial < handle
% AlignedTrial Class for handling details of a replay trial which has been 
%
%   AlignedTrial Properties:
%       probe_id            - string with the probe recording ID the trial was part of
%       session_id          - string with the session ID the trial was part of
%       trial_id            - integer with the trial ID
%       trial_group_label   - string ID of the trial group, see RVTSession
%       original_trial_id   - integer ID of original trial
%       protocol            - 'Coupled', 'EncoderOnly', 'ReplayOnly', 'StageOnly', 'CoupledMismatch', 'EncoderOnlyMismatch'
%       fs                  - sample rate of the NIDAQ session recording
%       trial               - replay Trial object
%       original_trial      - original trial Trial object
%       offset              - sample offset
%       n_points            - number of sample points from `offset` to take for the traces
%       is_replay           - true if trial is a replay of another trial
%       replay_of           - if trial is a replay, this is a string containing the protocol
%                             name of the original trial.
%
%   AlignedTrial Methods:
%       config                  - return structure containing the configuration for
%                                 this trial (loaded from .cfg file associated with session .bin file)
%       probe_t                 - return timebase of the aligned trial synced to the probe recording
%       velocity                - return the speed for the aligned trial
%       treadmill_speed         - return the treadmill speed for the aligned trial
%       visual_speed            - return the command speed to the visual corridor during this trial
%       filtered_teensy         - return the filtered_teensy channel for the aligned trial
%       filtered_teensy_2       - return the filtered_teensy_2 channel for the aligned trial
%       stage                   - return the stage channel for the aligned trial
%       multiplexer_output      - return the multiplexer_output channel for the aligned trial
%       teensy_gain             - return the teensy_gain channel for the aligned trial
%       camera1                 - return the motion energy for camera1 channel for the aligned trial
%       solenoid                - return the solenoid channel for the aligned trial
%       analysis_window         - return the analysis window for the aligned trial
%       treadmill_motion_mask   - return the motion mask of the treadmill for the aligned trial
%       stationary_mask         - return the stationary mask for the aligned trial
%       motion_mask             - return the motion mask for the aligned trial
%       stationary_time         - amount of time spent stationary for the aligned trial
%       motion_time             - amount of time spent in motion for the aligned trial
%       analysis_window_time    - amount of time in the analysis window for the aligned trial
%       motion_bouts            - calculate and return motion bouts in the aligned trial
%       get_aligned_data        - shared function for getting traces from the original trial but aligned
%   
%   AlignedTrial is a class for handling replay trials and aligning them to
%   their original trials. Any trial can be 'cast' to AlignedTrial, even
%   trials which are not replays.
%
%   See also: Trial
%
%   TODO: refactor this and Trial so that replay trials are always aligned
%           if possible
    
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
        % AlignedTrial
        %
        %   AlignedTrial(REPLAY_TRIAL, ORIGINAL_TRIAL) creates an object to
        %   handle replay trials which are aligned to the trial from which
        %   they are replayed. REPLAY_TRIAL is an object of class Trial
        %   associated with a replay trial, and ORIGINAL_TRIAL is also of
        %   class trial associated with the original trial.
        
            obj.trial = trial;
            obj.original_trial = original_trial;
            
            obj.fs = obj.trial.fs;
            obj.offset = obj.trial.replay_offset;
            obj.n_points = min(length(trial.probe_t) - obj.offset, length(original_trial.probe_t));
        end
        
        
        
        function val = get.probe_id(obj)
        %%string with the probe recording ID the trial was part of 
            val = obj.trial.probe_id;
        end
        
        
        
        function val = get.session_id(obj)
        %%string with the session ID the trial was part of
            val = obj.trial.session_id;
        end
        
        
        
        function val = get.trial_id(obj)
        %%integer with the trial ID    
            val = obj.trial.trial_id;
        end
        
        
        
        function val = get.original_trial_id(obj)
        %%integer ID of original trial    
            val = obj.trial.original_trial_id;
        end
        
        
        
        function val = get.trial_group_label(obj)
        %%string ID of the trial group, see RVTSession    
            val = obj.trial.trial_group_label;
        end
        
        
        
        function val = get.protocol(obj)
        %%'Coupled', 'EncoderOnly', 'ReplayOnly', 'StageOnly', 'CoupledMismatch', 'EncoderOnlyMismatch'
            val = obj.trial.protocol;
        end
        
        
        
        function val = get.is_replay(obj)
        %%true if trial is a replay of another trial
            val = obj.trial.is_replay;
        end
        
        
        
        function val = get.replay_of(obj)
        %%if trial is a replay, this is a string containing the protocol name of the original trial
            val = obj.trial.replay_of;
        end
        
        
        
        function val = get.replayed_trial_id(obj)
        %%TODO: this method is not needed   
            val = obj.trial.replayed_trial_id;
        end
        
        
        
        function val = config(obj)
        %%config Return structure containing the configuration for this trial
        %
        %   CONFIG = config()
        %
        %   See also: Trial
        
            val = obj.trial.config;
        end
        
        
        
        function val = probe_t(obj, idx)
        %%probe_t Return timebase of the aligned trial synced to the probe recording
        %
        %   TIMEBASE = probe_t(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            val = obj.trial.probe_t(obj.offset + (1:obj.n_points));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function val = velocity(obj, idx)
        %%velocity Return the speed for the aligned trial
        %
        %   SPEED = velocity(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            val = obj.trial.velocity(obj.offset + (1:obj.n_points));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        function val = acceleration(obj, idx)
        %%velocity Return the speed for the aligned trial
        %
        %   SPEED = velocity(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            val = obj.trial.acceleration(obj.offset + (1:obj.n_points));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        function val = treadmill_speed(obj, idx)
        %%treadmill_speed Return the treadmill speed for the aligned trial
        %
        %   SPEED = treadmill_speed(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
            
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.treadmill_speed)
                val = obj.trial.treadmill_speed(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = visual_speed(obj, idx)
        %%visual_speed Return the command speed to the visual corridor during this trial
        %
        %   SPEED = visual_speed(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            val = obj.get_aligned_data('visual_speed', idx);
        end
        
        
        
        function val = filtered_teensy(obj, idx)
        %%filtered_teensy Return the filtered_teensy channel for the aligned trial
        %
        %   SPEED = filtered_teensy(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.filtered_teensy)
                val = obj.trial.filtered_teensy(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = filtered_teensy_2(obj, idx)
        %%filtered_teensy_2 Return the filtered_teensy_2 channel for the aligned trial
        %
        %   SPEED = filtered_teensy_2(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.filtered_teensy_2)
                val = obj.trial.filtered_teensy_2(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = stage(obj, idx)
        %%stage Return the stage channel for the aligned trial
        %
        %   SPEED = stage(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.stage)
                val = obj.trial.stage(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = multiplexer_output(obj, idx)
        %%multiplexer_output Return the multiplexer_output channel for the aligned trial
        %
        %   SPEED = multiplexer_output(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            val = obj.get_aligned_data('multiplexer_output', idx);
        end
        
        
        
        function val = teensy_gain(obj, idx)
        %%teensy_gain Return the teensy_gain channel for the aligned trial
        %
        %   TRACE = teensy_gain(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.teensy_gain)
                val = obj.trial.teensy_gain(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        
        function val = camera1(obj, idx)
        %%camera1 Return the motion energy for camera1 channel for the aligned trial
        %
        %   MOTION_ENERGY = camera1(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.camera1)
                val = obj.trial.camera1(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        function val = camera0(obj, idx)
        %%camera1 Return the motion energy for camera1 channel for the aligned trial
        %
        %   MOTION_ENERGY = camera1(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.camera0)
                val = obj.trial.camera0(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end

        function val = pupil_diameter(obj, idx)
        %%TODO
        
            VariableDefault('idx', []);
            
            if ~isempty(obj.trial.pupil_diameter)
                val = obj.trial.pupil_diameter(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            else
                val = [];
            end
        end
        
        
        function val = solenoid(obj, idx)
        %%solenoid Return the solenoid channel for the aligned trial
        %
        %   SOLENOID = solenoid(INDEX)
        %   If INDEX is supplied it should be a boolean vector of length
        %   `n_points` indicating which part of the trace to return. If not
        %   supplied or empty, the entire trace is returned.
        %
        %   See also: Trial
        
            VariableDefault('idx', []);
            
            val = obj.trial.solenoid(obj.offset + (1:obj.n_points));
            
            if ~isempty(idx), val = val(idx); end
        end
        
        
        
        function mask = analysis_window(obj)
        %%analysis_window Return the analysis window for the aligned trial
        %
        %   MASK = analysis_window()
        %
        %   See also: Trial
        
            mask = obj.original_trial.analysis_window();
            mask = mask(1:obj.n_points);
        end
        
        
        
        function mask = treadmill_motion_mask(obj, vel_thresh, acc_thresh, min_dur)
        %%treadmill_motion_mask Return the motion mask of the treadmill for the aligned trial
        %
        %   MASK = treadmill_motion_mask(VELOCITY_THRESHOLD, ACCELERATION_THRESHOLD, MIN_STATIONARY_DURATION)
        %
        %   For definition of arguments see Trial.treadmill_motion_mask
        %
        %   See also: Trial.treadmill_motion_mask
        
            VariableDefault('vel_thresh', 1);
            VariableDefault('acc_thresh', 0.5);
            VariableDefault('min_dur', 0.2);
        
            mask = obj.trial.treadmill_motion_mask(vel_thresh, acc_thresh, min_dur);
            mask = mask(obj.offset + (1:obj.n_points));
        end
        
        
        function mask = stationary_mask(obj)
        %%stationary_mask Return the stationary mask for the aligned trial
        %
        %   MASK = stationary_mask()
        %
        %   See also: Trial
        
            mask = obj.original_trial.stationary_mask;
            mask = mask(1:obj.n_points);
%             mask = mask(obj.offset + (1:obj.n_points));
%             mask = mask & obj.analysis_window();
        end
        
        
        
        function mask = motion_mask(obj)
        %%motion_mask Return the motion mask for the aligned trial
        %
        %   MASK = motion_mask()
        %
        %   See also: Trial
        
            mask = obj.original_trial.motion_mask();
            mask = mask(1:obj.n_points);
%             mask = mask & ~obj.stationary_mask();
        end
        
        
        
        function val = stationary_time(obj)
        %%stationary_time Amount of time spent stationary for the aligned trial
        %
        %   TIME = stationary_time()
        %
        %   See also: Trial
        
            val = sum(obj.stationary_mask()) / obj.fs;
        end
        
        
        
        function val = motion_time(obj)
        %%motion_time Amount of time spent in motion for the aligned trial
        %
        %   TIME = motion_time()
        %
        %   See also: Trial
        
            val = sum(obj.motion_mask()) / obj.fs;
        end
        
        
        
        function val = analysis_window_time(obj)
        %%analysis_window_time Amount of time in the analysis window for the aligned trial
        %
        %   TIME = analysis_window_time()
        %
        %   See also: Trial
        
            val = sum(obj.analysis_window()) / obj.fs;
        end
        
        
        
        function bouts = motion_bouts(obj, include_200ms, true_motion_start)
        %%motion_bouts Calculate and return motion bouts in the trial
        %
        %   BOUTS_STRUCT = motion_bouts(INCLUDE_200MS, TRUE_MOTION_START)
        %
        %   For definition of arguments and defaults see
        %   Trial.motion_bouts.
        %
        %   See also: Trial.motion_bouts
        
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
        %%get_aligned_data Shared function for getting traces from the original trial but aligned
        %
        %   TRACE = get_aligned_data(STRING, INDEX) shared function for
        %   getting traces for the aligned trial. STRING refers to a
        %   property name (e.g. 'visual', 'multiplexer_output'). INDEX
        %   should be supplied should be a boolean vector of length  `n_points`
        %   indicating which part of the trace to return.
        %
            val = [];
            if ~isempty(obj.trial.(type))
                val = obj.trial.(type)(obj.offset + (1:obj.n_points));
                if ~isempty(idx), val = val(idx); end
            end
        end
    end
end
