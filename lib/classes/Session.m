classdef Session < handle
    
    properties
        
        id
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
        
        trials = Trial.empty()
    end
    
    methods
        
        function obj = Session(session, t_sync, cameras, t_trig)
            % takes a session structure and converts it to a Session object
            
            VariableDefault('cameras', []);
            VariableDefault('t_trig', []);
            
            obj.id = session.id;
            obj.fs = session.fs;
            obj.config = session.config;
            
            chan_names = {'filtered_teensy', 'raw_teensy', ...
                          'stage', 'lick', 'pump', 'solenoid', ...
                          'photodiode', 'gain_change', 'multiplexer_output', 'minidaq_ao0', ...
                          'teensy_gain', 'gain_teensy', 'filtered_teensy_2'};
            
            for chan_i = 1 : length(chan_names)
                if isfield(session, chan_names{chan_i})
                    obj.(chan_names{chan_i}) = session.(chan_names{chan_i});
                end
            end
            
            obj.rc2_t = (0:length(obj.solenoid)-1)'*(1/obj.fs);
            
            % handle case where t_sync is missing
            if isempty(t_sync)
                obj.probe_t = nan(length(obj.rc2_t), 1);
            else
                obj.probe_t = t_sync;
            end
            
            for trial_i = 1 : length(session.trials)
                obj.trials(trial_i) = Trial(session.trials(trial_i), obj);
            end
            
            if ~isempty(cameras)
                obj.camera0 = interp1(t_trig, cameras.camera0, obj.probe_t);
                obj.camera1 = interp1(t_trig, cameras.camera1, obj.probe_t);
                obj.camera_motion_mask = compute_camera_motion_mask(obj.camera1, false);
            end
        end
        
        
        function trials = trials_by_protocol(obj, protocol)
            
            protocol_list = {obj.trials(:).protocol};
            
            idx = strcmp(protocol_list, protocol);
            
            trials = obj.trials(idx);
        end
    end
end