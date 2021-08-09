classdef Anatomy < handle
    
    properties
        
        probe_fname
        anatomy
        
    end
    
    properties (Constant = true)
        
        VISp_layers = {'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
        
    end
    
    
    methods
        
        function obj = Anatomy(data)
            
            obj.probe_fname = data.data.probe_recording;
            obj.anatomy = data.data.anatomy;
            
        end
        
        
        function boundaries = VISp_boundaries(obj)
            
            r = 0;
            
            for i = 1 : length(obj.anatomy.region_str)
                
                if ~any(strcmp(obj.VISp_layers, obj.anatomy.region_str{i}))
                    continue
                end
                
                r = r + 1;
                
                this_instance = 1;
                
                if r > 1
                    
                    idx = strcmp(boundaries.region(1:r-1), obj.anatomy.region_str{i});
                    
                    if sum(idx) > 0
                        this_instance = sum(idx) + 1;
                    end
                end
                
                boundaries.region{r, 1} = obj.anatomy.region_str{i};
                boundaries.upper(r, 1) = obj.anatomy.region_boundaries(i);
                boundaries.lower(r, 1) = obj.anatomy.region_boundaries(i+1);
                boundaries.instance(r, 1) = this_instance;
                
            end
            
        end
        
        
        function [relative_depth, layer] = VISp_layer_relative_depth(obj, from_tip)
            
            boundaries = obj.VISp_boundaries();
            
            idx = find(from_tip < boundaries.upper & from_tip > boundaries.lower);
            
            % if no region found
            if isempty(idx)
                layer = '';
                relative_depth = nan;
                return
            end
            
            relative_depth = 1 - ((from_tip - boundaries.lower(idx)) / ...
                (boundaries.upper(idx) - boundaries.lower(idx)));
            layer = boundaries.region{idx};
        end
    end
end