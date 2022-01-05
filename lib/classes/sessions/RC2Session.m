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
%
%       rc2_t               - timebase of the trial in the NIDAQ session recording (0s on first sample point, 1/fs on second etc.)
%       probe_t             - timebase of the trial synced to the probe recording
%       camera_t            - timebase of the camera frames synced to the probe recording
%
%   RC2Session Methods:
%       camera0_interp      - interpolate the motion energy data on camera0 so it has the same timebase as the other channel.
%       camera1_interp      - interpolate the motion energy data on camera1 so it has the same timebase as the other channel.
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
        
        rc2_t
        probe_t
    end
    
    properties (SetAccess = private, Hidden = true)
        
        camera0_
        camera1_
        camera_t
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
    end
end
