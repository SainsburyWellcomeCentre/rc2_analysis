classdef RC2Trial < handle
% Currently unused

    properties (SetAccess = private)
        
        id
        session_id
        probe_id
        config
    end
    
    properties (Dependent = true)
        
        fs
        
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
    
    
    properties (SetAccess = private)
        
        start_idx
        end_idx
    end
    
    
    properties (Hidden = true)
        
        trial
        session
    end
    
    
    
    methods
        
        function obj = RC2Trial(trial, session)
        
            obj.trial = trial;
            obj.session = session;
            
            obj.id = trial.id;
            obj.session_id = trial.session_id;
            obj.probe_id = session.probe_id;
            obj.config = trial.config;
            
            obj.start_idx = trial.start_idx;
            obj.end_idx = trial.end_idx;
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
            val = [];
            if ~isempty(obj.session.camera0_interp)
                val = obj.session.camera0_interp(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.camera1(obj)
            val = [];
            if ~isempty(obj.session.camera1_interp)
                val = obj.session.camera1_interp(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.camera_motion_mask(obj)
            val = [];
            if ~isempty(obj.session.camera_motion_mask)
                val = obj.session.camera_motion_mask(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.teensy_gain(obj)
            val = [];
            if ~isempty(obj.session.teensy_gain)
                val = obj.session.teensy_gain(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.filtered_teensy_2(obj)
            val = [];
            if ~isempty(obj.session.filtered_teensy_2)
                val = obj.session.filtered_teensy_2(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.raw_teensy(obj)
            val = [];
            if ~isempty(obj.session.filtered_teensy_2)
                val = obj.session.raw_teensy(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.multiplexer_output(obj)
            val = [];
            if ~isempty(obj.session.multiplexer_output)
                val = obj.session.multiplexer_output(obj.start_idx:obj.end_idx);
            end
        end
        
        
        
        function val = get.gain_teensy(obj)
            val = [];
            if ~isempty(obj.session.filtered_teensy_2)
                val = obj.session.gain_teensy(obj.start_idx:obj.end_idx);
            end
        end
    end
end
