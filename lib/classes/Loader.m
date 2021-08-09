classdef Loader < handle
    
    properties
        
        file_manager
    end
    
    
    methods
        
        function obj = Loader(config)
            
            obj.file_manager = FileManager(config);
        end
        
        
        
        function data = formatted_data(obj, recording_id, variables)
            
            VariableDefault('variables', []);
            
            fname = obj.file_manager.formatted_data(recording_id);
            
            if isempty(variables)
                
                data = load_data(fname);
                
                data.probe_recording = recording_id;
                
                s_list = obj.session_list();
                
                for i = length(data.sessions):-1:1
                    
                    idx = strcmp(s_list.session_id, data.sessions(i).id);
                    
                    if sum(idx) == 0
                        exp_type{i} = 'unknown';
                    else
                        exp_type{i} = s_list.short_name{idx};
                    end
                end
                
                data.experiment_type = exp_type;
                data = DataController(data);
                
            else
                
                data = load(fname, variables{:});
            end
        end
        
        
        
        function s_list = session_list(obj)
            
            fname = obj.file_manager.session_list();
            s_list = readtable(fname);
        end
    end
end
