classdef Loader < handle
    
    properties
        
        file_manager
    end
    
    
    methods
        
        function obj = Loader(file_manager)
            
            obj.file_manager = file_manager;
        end
        
        
        
        function data = formatted_data(obj, recording_id, variables)
            
            VariableDefault('variables', []);
            
            fname = obj.file_manager.formatted_data(recording_id);
            
            if isempty(variables)
                
                data = load_data(fname);
                
                data.probe_recording = recording_id;
                3
                s_list = obj.session_list();
                
                for i = length(data.sessions):-1:1
                    
                    idx = strcmp(s_list.session_id, data.sessions(i).id);
                    
                    if sum(idx) == 0
                        exp_type{i} = 'unknown';
                    else
                        exp_type{i} = s_list.short_name{idx};
                    end
                end
                
                data.experiment_group = exp_type;
                data = DataController(data);
                
            else
                
                data = load(fname, variables{:});
            end
        end
        
        
        
        function s_list = session_list(obj)
            
            fname = obj.file_manager.session_list();
            s_list = readtable(fname);
        end
        
        
        function e_list = experiment_list(obj)
            
            [fname, exists] = obj.file_manager.experiment_list();
            if exists
                e_list = readtable(fname);
            else
                e_list = [];
            end
        end
    end
end
