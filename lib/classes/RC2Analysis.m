classdef RC2Analysis < handle
    
    properties
        
        path_config
        
        file
        load
    end
    
    
    methods
        
        function obj = RC2Analysis()
        %%main class for accessing data from RC2 setup
            obj.path_config     = path_config();
            
            obj.file            = FileManager(obj.path_config);
            obj.load            = Loader(obj.file);
        end
        
        
        
        function val = get_probe_ids(obj, experiment_group)
        %%get recording ids for a named group of experiments
            experiment_list = obj.load.experiment_list();
            idx = strcmp(experiment_list.experiment_group, experiment_group);
            val = experiment_list.probe_id(idx);
        end
        
        
        
        function data = load_formatted_data(obj, probe_id)
        %%loads the formatted data for a probe recording
            data = obj.load.formatted_data(probe_id);
        end
    end
end
