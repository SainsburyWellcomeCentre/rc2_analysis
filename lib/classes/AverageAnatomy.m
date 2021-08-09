classdef AverageAnatomy < handle
    
    properties
        
        anatomy_array
        
    end
    
    properties (Constant = true)
        
        VISp_layers = {'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
        
    end
    
    methods
        
        function obj =  AverageAnatomy(anatomy_array)
            
            
            obj.anatomy_array = anatomy_array;
            
        end
        
        
        function boundaries = average_VISp_boundaries(obj)
            
            for i = 1 : length(obj.anatomy_array)
                
                b(i) = obj.anatomy_array(i).VISp_boundaries;
                
                assert(isequal(b(i).region, obj.VISp_layers'));
            end
            
            boundaries.region = obj.VISp_layers;
            boundaries.upper = mean(cat(2, b(:).upper), 2);
            boundaries.lower = mean(cat(2, b(:).lower), 2);            
            
        end
        
        
        function [boundary_position, cluster_position] = mi_vs_depth_positions(obj, rel_pos, layer)
            
            boundaries = obj.average_VISp_boundaries();
            
            cluster_position = nan(length(rel_pos), 1);
            
            for i = 1 : length(rel_pos)
                
                idx = find(strcmp(layer{i}, boundaries.region));
                cluster_position(i) = boundaries.upper(idx) - (boundaries.upper(idx) - boundaries.lower(idx))*rel_pos(i);
            end
            
            boundary_position = [boundaries.upper(:); boundaries.lower(end)];
            
        end
        
    end
end