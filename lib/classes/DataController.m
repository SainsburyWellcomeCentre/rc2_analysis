classdef DataController < handle
    
    properties
        
        experiment_type
        data
    end
    
    methods
        
        function obj = DataController(data, ~)
            
            obj.data = data;
            obj.experiment_type = data.experiment_type{1};
        end
        
        
        function selected_clusters = selected_clusters(obj)
            
            selected_clusters = obj.data.clusters;
            idx = ismember([obj.data.clusters(:).id], obj.data.selected_clusters);
            selected_clusters(~idx) = [];
            
        end
        
        
        function visp_clusters = VISp_clusters(obj, is_selected, spiking_class)
            
            VariableDefault('is_selected', true);
            VariableDefault('spiking_class', 'any');
            
            if is_selected
                selected_clusters = obj.selected_clusters();
            else
                selected_clusters = obj.data.clusters;
            end
            
            idx = regexp({selected_clusters(:).region_str}, 'VISp[\dX]');
            idx = ~cellfun(@isempty, idx);
            
            visp_clusters = selected_clusters(idx);
            
            if strcmp(spiking_class, 'RS')
                idx = [visp_clusters(:).duration] >= 0.45;
                visp_clusters = visp_clusters(idx);
            elseif strcmp(spiking_class, 'FS')
                idx = [visp_clusters(:).duration] < 0.45;
                visp_clusters = visp_clusters(idx);
            end
            
        end
        
        
        function trials = mismatch_trials(obj)
            
            % get mismatch trials
            if strcmp(obj.data.probe_recording, 'CAA-1112872_rec1_rec1b_rec2_rec3')
                trials = [obj.data.sessions(1).trials, obj.data.sessions(2).trials];
            else
                trials = [obj.data.sessions(1).trials];
            end
            
        end
        
        
        function fs = sample_rate(obj)
            
            fs = obj.data.sessions(1).fs;
            
        end
        
    end
    
end