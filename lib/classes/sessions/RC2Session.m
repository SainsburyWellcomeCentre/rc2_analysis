classdef RC2Session < handle
    
    properties
        
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
        
        camera_motion_mask
        
        rc2_t
        probe_t
        camera_t
    end
    
    
    
    methods
        
        function obj = RC2Session(session)
            % takes a session structure and converts it to a Session object
            
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
                obj.camera0 = obj.camera0_interp();
                obj.camera1 = obj.camera1_interp();
            end
        end
        
        
        
        function val = camera0_interp(obj)
            
            val = [];
            if ~isempty(obj.camera0)
                val = interp1(obj.camera_t, obj.camera0, obj.probe_t);
            end
        end
        
        
        
        function val = camera1_interp(obj)
            
            val = [];
            if ~isempty(obj.camera1)
                val = interp1(obj.camera_t, obj.camera1, obj.probe_t);
            end
        end
        
        
        
        function compute_camera_motion_mask(obj)
            obj.camera_motion_mask = compute_camera_motion_mask(obj.camera1_interp, false);
        end
    end
end