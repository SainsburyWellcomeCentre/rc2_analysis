classdef AverageAnatomy < handle
    
    properties
        
        anatomy_array
    end
    
    properties (Constant = true)
        
        VISp_layers = {'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
    end
    
    
    
    methods
        
        function obj =  AverageAnatomy(anatomy_array)
        %%class for combining a set of anatomies to get averaged boundaries
            obj.anatomy_array = anatomy_array;
        end
        
        
        
        function boundaries = average_VISp_boundaries_from_pia_struct(obj)
        %%computes the average distance of the boundaries of each VISp layer when averaged across several anatomies    
            
            for ii = 1 : length(obj.anatomy_array)
                
                b(ii) = obj.anatomy_array{ii}.VISp_boundaries_from_pia;
                
                % make sure that none of the anatomies are 'strange' (i.e.
                % have multiple of the same layers)
                assert(isequal(b(ii).region, obj.VISp_layers'), 'anatomy likely has multiple versions of the same layer');
            end
            
            boundaries.region = obj.VISp_layers;
            boundaries.upper = mean(cat(2, b(:).upper), 2);
            boundaries.lower = mean(cat(2, b(:).lower), 2);
        end
        
        
        
        function boundary_position = average_VISp_boundaries_from_pia(obj)
            
            boundaries = obj.average_VISp_boundaries_from_pia_struct();
            boundary_position = [boundaries.upper(:); boundaries.lower(end)];
        end
        
%         function [boundary_position, cluster_position] = mi_vs_depth_positions(obj, rel_pos, layer)
%             
%             boundaries = obj.average_VISp_boundaries();
%             
%             cluster_position = nan(length(rel_pos), 1);
%             
%             for i = 1 : length(rel_pos)
%                 
%                 idx = find(strcmp(layer{i}, boundaries.region));
%                 cluster_position(i) = boundaries.upper(idx) - (boundaries.upper(idx) - boundaries.lower(idx))*rel_pos(i);
%             end
%             
%             boundary_position = [boundaries.upper(:); boundaries.lower(end)];
%             
%         end
        
        
        
        
    
        function cluster_position = from_pia_using_relative_position(obj, relative_position, layer)
        %%takes a set of relative positional points (each a point within a VISp layer, relative to the upper boundary), along
        %%with the layer, to compute the position of the point in a map of
        %%boudaries averaged across several anatomies
            
            % get averaged boundaries
            boundaries = obj.average_VISp_boundaries_from_pia_struct();
            
            cluster_position = nan(length(relative_position), 1);
            
            for i = 1 : length(relative_position)
                
                % get region
                idx = find(strcmp(layer{i}, boundaries.region));
                cluster_position(i) = boundaries.upper(idx) + (boundaries.lower(idx) - boundaries.upper(idx))*relative_position(i);
            end
            
            boundary_position = [boundaries.upper(:); boundaries.lower(end)];
        end
    end
end