classdef Anatomy < handle
    
    properties (SetAccess = private)
        
        anatomy
    end
    
    properties (Dependent = true)
        
        probe_id
    end
    
    properties (Constant = true)
        
        VISp_layers = {'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
    end
    
    
    
    methods
        
        function obj = Anatomy(anatomy)
            
            obj.anatomy = anatomy;
        end
        
        
        
        function val = get.probe_id(obj)
            val = obj.anatomy.probe_id;
        end
        
        
        
        function from_pia = from_tip_to_from_pia(obj, from_tip)
        %%take a value in um from the probe tip and convert to um from the
        %%pial surface
            % find pial surface from tip
            visp1_cmp =  regexp(obj.anatomy.region_str, 'VISp*\w1');
            visp1_idx = find(cellfun(@(x)(~isempty(x)), visp1_cmp));
            pial_surface_from_tip = obj.anatomy.region_boundaries(visp1_idx);
            from_pia = pial_surface_from_tip - from_tip;
        end
        
        
        
        function boundaries = VISp_boundaries(obj)
        %%returns a structure containing information about the boundaries
        %%of VISp layers for a single anatomy (probe track)
        %   boundaries.region   - the VISp layer
        %             .upper    - the upper boundary of the layer (from probe_tip)
        %             .lower    - the lower boundary of the layer (from probe_tip)
        %             .instance - if the probe goes in and out of VISp (can
        %                       happen) returns which instance of this VISp
        %                       layer we have come across (working downward
        %                       from the pia)
        
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
                
                boundaries.region{r, 1}     = obj.anatomy.region_str{i};
                boundaries.upper(r, 1)      = obj.anatomy.region_boundaries(i);
                boundaries.lower(r, 1)      = obj.anatomy.region_boundaries(i+1);
                boundaries.instance(r, 1)   = this_instance;
            end
        end
        
        
        
        function boundaries = VISp_boundaries_from_pia(obj)
        %%as for VISp_boundaries() but boundary values are *from pial
        %%surface* and not from probe_tip
            boundaries = obj.VISp_boundaries();
            boundaries.upper = obj.from_tip_to_from_pia(boundaries.upper);
            boundaries.lower = obj.from_tip_to_from_pia(boundaries.lower);
        end
        
        
        
        function boundaries = VISp_boundaries_from_pia_flat(obj)
        %%as for VISp_boundaries() but boundary values are *from pial
        %%surface* and not from probe_tip
            b = obj.VISp_boundaries_from_pia();
            boundaries = [b.upper; b.lower(end)];
        end
        
        
        
        function [relative_depth, layer] = VISp_layer_relative_depth(obj, from_tip)
        %%get relative depth of a point in a VISp layer
        %   from_tip = um of point from tip
        %  the relative value is FROM THE UPPER BOUNDARY
        %  if point not in VISp returns [nan, '']
        
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