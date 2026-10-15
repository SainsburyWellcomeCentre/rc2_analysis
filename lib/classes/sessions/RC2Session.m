classdef RC2Session < handle
% RC2Session Class for handling details of a session
%
%   RC2Session Properties:
%       probe_id            - string with the probe recording ID the session was part of
%       session_id          - string with the session ID
%       fs                  - sample rate of the NIDAQ session recording
%       config              - structure containing the configuration for
%                             this trial (loaded from .cfg file associated with session .bin file)
%
%       gain_teensy         - gain_teensy channel
%       filtered_teensy     - filtered_teensy channel
%       filtered_teensy_2   - filtered_teensy_2 channel
%       raw_teensy          - raw_teensy channel
%       stage               - speed of the stage (c.f. stage_speed)
%       lick                - lick channel
%       pump                - pump channel
%       solenoid            - solenoid channel
%       photodiode          - photodiode channel
%       gain_change         - gain_change channel
%       minidaq_ao0         - minidaq_ao0 channel
%       multiplexer_output  - multiplexer_output channel
%       teensy_gain         - teensy_gain channel
%       camera0             - motion energy for camera0
%       camera1             - motion energy for camera1
%       pupil_diameter      - pupil diameter (pixels) from the eye camera
%       pupil_position      - pupil x-y position (pixels, N x 2) from the eye camera
%
%       rc2_t               - timebase of the trial in the NIDAQ session recording (0s on first sample point, 1/fs on second etc.)
%       probe_t             - timebase of the trial synced to the probe recording
%       camera_t            - timebase of the camera frames synced to the probe recording
%       eyecam_t            - timebase of the eye camera frames (eyecam_0 pulses) synced to the probe recording
%
%   RC2Session Methods:
%       camera0_interp         - interpolate the motion energy data on camera0 so it has the same timebase as the other channel.
%       camera1_interp         - interpolate the motion energy data on camera1 so it has the same timebase as the other channel.
%       pupil_diameter_interp  - interpolate the pupil diameter data so it has the same timebase as the other channel.
%       pupil_position_interp  - interpolate the pupil x-y position data so it has the same timebase as the other channel.
%       eyecam_t_from_pulses   - timebase of the eye camera frames from the eyecam_0 pulses
%       align_to_eyecam        - match eye camera frames to the eyecam_0 pulses
%
%   See also: RVTSession

    properties (SetAccess = private)
        
        probe_id
        session_id
        fs
        config
        
        gain_teensy
        filtered_teensy
        filtered_teensy_2
        raw_teensy
        stage
        lick
        pump
        solenoid
        photodiode
        gain_change
        minidaq_ao0
        multiplexer_output
        teensy_gain
        camera0
        camera1
        pupil_diameter
        pupil_position
        
        rc2_t
        probe_t
    end
    
    properties (SetAccess = private, Hidden = true)
        
        camera0_
        camera1_
        camera_t
        eyecam_t
        pupil_diameter_
        pupil_position_
    end
    
    properties (Constant = true, Hidden = true)
        
        eyecam_threshold = 1.5  % volts, eyecam_0 pulses are ~0-3.3V
    end
    
    
    
    methods
        
        function obj = RC2Session(session)
        % RC2Session
        %
        %   RC2Session(SESSION) takes a session structure saved in the
        %   formatted data and creates an object to handle the data.
            
            chan_names = {'session_id', ...
                          'probe_id', ...
                          'fs', ...
                          'config', ...
                          'filtered_teensy', ...
                          'raw_teensy', ...
                          'stage', ...
                          'lick', ...
                          'pump', ...
                          'solenoid', ...
                          'photodiode', ...
                          'gain_change', ...
                          'multiplexer_output', ...
                          'minidaq_ao0', ...
                          'teensy_gain', ...
                          'gain_teensy', ...
                          'filtered_teensy_2', ...
                          'camera0', ...
                          'camera1', ...
                          'pupil_diameter', ...
                          'pupil_position', ...
                          'rc2_t', ...
                          'probe_t', ...
                          'camera_t'};
            
            % don't include camera data if something wrong with recording
            if length(session.camera0) ~= length(session.camera_t)
                idx = cellfun(@(x)(any(ismember({'camera0', 'camera1', 'camera_t'}, x))), chan_names);
                chan_names(idx) = [];
            end
                      
            for chan_i = 1 : length(chan_names)
                if isfield(session, chan_names{chan_i})
                    obj.(chan_names{chan_i}) = session.(chan_names{chan_i});
                end
            end
            
            if ~isempty(obj.camera0)
                obj.camera0_ = obj.camera0;
                obj.camera1_ = obj.camera1;
                obj.camera0 = obj.camera0_interp();
                obj.camera1 = obj.camera1_interp();
            end

            % the pupil data come from the eye camera, a separate system
            % from camera0/camera1, which emits one pulse per frame on the
            % eyecam_0 channel
            obj.eyecam_t = obj.eyecam_t_from_pulses(session);
            obj.pupil_diameter_ = obj.align_to_eyecam(obj.pupil_diameter, 'pupil_diameter');
            obj.pupil_position_ = obj.align_to_eyecam(obj.pupil_position, 'pupil_position');
            obj.pupil_diameter = obj.pupil_diameter_interp();
            obj.pupil_position = obj.pupil_position_interp();
        end



        function val = eyecam_t_from_pulses(obj, session)
        %%eyecam_t_from_pulses Timebase of the eye camera frames
        %
        %   EYECAM_T = eyecam_t_from_pulses(SESSION)
        %   returns the time, synced to the probe recording, of each rising
        %   edge on the eyecam_0 channel of the session structure SESSION.
        %   Each edge marks one eye camera frame. Empty if there is no
        %   eyecam_0 channel.

            val = [];
            if ~isfield(session, 'eyecam_0') || isempty(session.eyecam_0) || isempty(obj.probe_t)
                return
            end

            rising_idx = find(diff(session.eyecam_0 > obj.eyecam_threshold) == 1) + 1;
            val = obj.probe_t(rising_idx);
            val = val(:);
        end



        function val = align_to_eyecam(obj, val, name)
        %%align_to_eyecam Match eye camera frames to the eyecam_0 pulses
        %
        %   ALIGNED = align_to_eyecam(VAL, NAME)
        %   takes VAL, an N x M matrix with one row per eye camera frame,
        %   and returns the first rows matching the number of pulses in
        %   `eyecam_t`. The eye camera keeps recording for a few frames
        %   after the NIDAQ stops, so extra frames are at the end and are
        %   dropped. If there are fewer frames than pulses, or no pulses,
        %   the frames cannot be matched and ALIGNED is empty. NAME is used
        %   in the printed messages.

            if isempty(val)
                return
            end

            n_frames = size(val, 1);
            n_pulses = length(obj.eyecam_t);

            if n_pulses == 0
                fprintf('[pupil] %s: no eyecam_0 pulses, %s left empty\n', obj.session_id, name);
                val = [];
            elseif n_frames < n_pulses
                warning('RC2Session:align_to_eyecam:tooFewFrames', ...
                    '%s: %s has %i frames but there are %i eyecam_0 pulses, cannot align, left empty', ...
                    obj.session_id, name, n_frames, n_pulses);
                val = [];
            else
                if n_frames > n_pulses
                    fprintf('[pupil] %s: %s has %i frames, %i eyecam_0 pulses, dropping last %i frames\n', ...
                        obj.session_id, name, n_frames, n_pulses, n_frames - n_pulses);
                end
                val = val(1:n_pulses, :);
            end
        end
        
        
        
        function val = camera0_interp(obj)
        %%camera0_interp Interpolate the motion energy data so it has the same
        %%timebase as the other channel.
        %
        %   MOTION_ENERGY = camera0_interp()
        %   increases the resolution of the motion energy trace from the
        %   camera data so it has the same resolution as the other NIDAQ
        %   traces.
        
            val = [];
            if ~isempty(obj.camera0)
                val = interp1(obj.camera_t, obj.camera0, obj.probe_t);
            end
        end
        
        
        
        function val = camera1_interp(obj)
        %%camera1_interp Interpolate the motion energy data so it has the same
        %%timebase as the other channel.
        %
        %   MOTION_ENERGY = camera1_interp()
        %   increases the resolution of the motion energy trace from the
        %   camera data so it has the same resolution as the other NIDAQ
        %   traces.
        
            val = [];
            if ~isempty(obj.camera1)
                val = interp1(obj.camera_t, obj.camera1, obj.probe_t);
            end
        end
        
        function val = pupil_diameter_interp(obj)
        %%pupil_diameter_interp Interpolate the pupil diameter data so it has the same
        %%timebase as the other channel.
        %
        %   DIAMETER = pupil_diameter_interp()
        %   increases the resolution of the pre-computed pupil diameter trace 
        %   so it has the same resolution as the other NIDAQ
        %   traces.
        
            val = [];
            if ~isempty(obj.pupil_diameter_)
                val = interp1(obj.eyecam_t, obj.pupil_diameter_, obj.probe_t);
            end
        end



        function val = pupil_position_interp(obj)
        %%pupil_position_interp Interpolate the pupil x-y position data so it has the same
        %%timebase as the other channel.
        %
        %   POSITION = pupil_position_interp()
        %   increases the resolution of the pre-computed pupil x-y position
        %   (N x 2) so it has the same resolution as the other NIDAQ
        %   traces.

            val = [];
            if ~isempty(obj.pupil_position_)
                val = interp1(obj.eyecam_t, obj.pupil_position_, obj.probe_t(:));
            end
        end
    end
end
