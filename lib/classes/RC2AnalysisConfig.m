classdef RC2AnalysisConfig < handle
    
    properties
        
        figure_dir = 'C:\Users\Lee\Documents\mvelez\figures'
        local_data_dir = 'D:\mvelez';
        formatted_data_dir
        summary_data_dir
    end
    
    methods
        
        function obj = RC2AnalysisConfig()
        end
        
        
        function val = get.formatted_data_dir(obj)
            
            val = fullfile(obj.local_data_dir, 'formatted_data');
        end
        
        
        function val = get.summary_data_dir(obj)
            
            val = fullfile(obj.local_data_dir, 'summary_data');
        end
        
        
        
        
        
    end
end