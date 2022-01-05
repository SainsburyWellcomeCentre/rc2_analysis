classdef Anatomy < handle
% Anatomy Class for handling the anatomy information (after loading anatomy
% structure from a formatted data file)
%
%   Anatomy Properties:
%       anatomy         - the `anatomy` structure from the formatted data file
%       probe_id        - string containing the probe recording ID
%       shank_id        - integer of the shank ID (zero-indexed)
%       VISp_layers     - cell array of strings containing acronyms for VISp layers
%
%   Anatomy Methods:
%       from_tip_to_from_pia    - convert a distance from the probe tip to a distance from the pia (in um)
%       VISp_boundaries         - return structure with boundary distances for the VISp layers from probe tip (in um)
%       VISp_boundaries_from_pia - return structure with boundary distances for the VISp layers from pia (in um)
%       VISp_boundaries_from_pia_flat - return vector with boundary distances for the VISp layers from pia (in um)
%       VISp_layer_relative_depth - return relative depth in a layer, for a point given from probe tip (in um)
%
%   See also: ProbeTrack

    properties (SetAccess = private)
        
        anatomy
    end
    
    properties (Dependent = true)
        
        probe_id
        shank_id
    end
    
    properties (Constant = true)
        
        VISp_layers = {'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
    end
    
    
    
    methods
        
        function obj = Anatomy(anatomy)
        % Anatomy
        %
        %   Anatomy(ANATOMY) creates an object to handle anatomy
        %   information for a probe recording. Takes as argument ANATOMY, a
        %   structure saved in the formatted data file.
        %
        %   See also: ProbeTrack
        
            obj.anatomy = anatomy;
        end
        
        
        
        function val = get.probe_id(obj)
        %%string of probe recording ID
            val = obj.anatomy.probe_id;
        end
        
        
        
        function val = get.shank_id(obj)
        %%integer of shank ID (zero-indexed)
            val = obj.anatomy.shank_id;
        end
        
        
        
        function from_pia = from_tip_to_from_pia(obj, from_tip)
        %%from_tip_to_from_pia Convert a distance from the probe tip to a distance from the pia (in um)
        %
        %   FROM_PIA = from_tip_to_from_pia(FROM_TIP)
        %   take a value in um from the probe tip, FROM_TIP, and convert to
        %   um from the pial surface, FROM_PIA.
        
            % find pial surface from tip
            visp1_cmp =  regexp(obj.anatomy.region_str, 'VISp*\w1');
            visp1_idx = find(cellfun(@(x)(~isempty(x)), visp1_cmp));
            pial_surface_from_tip = obj.anatomy.region_boundaries(visp1_idx);
            from_pia = pial_surface_from_tip - from_tip;
        end
        
        
        
        function boundaries = VISp_boundaries(obj)
        %%VISp_boundaries Return structure with boundary distances for the VISp layers from probe tip (in um)
        %
        %   BOUNDARIES = VISp_boundaries()
        %   returns a structure, BOUNDARIES, containing information about the boundaries
        %   of VISp layers for a single anatomy (probe track). Fields are:
        %       region   - string with the VISp layer
        %       upper    - the upper boundary of the layer (from probe_tip)
        %       lower    - the lower boundary of the layer (from probe_tip)
        %       instance - if the probe goes in and out of a VISp layber (can
        %                  happen at VISp edges) returns which instance of this VISp
        %                  layer we have come across (working downward
        %                   from the pia)
        
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
        %%VISp_boundaries_from_pia Return structure with boundary distances for the VISp layers from pia (in um)
        %
        %   BOUNDARIES = VISp_boundaries_from_pia()
        %   returns a structure, BOUNDARIES, containing information about the boundaries
        %   of VISp layers for a single anatomy (probe track). Fields are:
        %       region   - string with the VISp layer
        %       upper    - the upper boundary of the layer (from pia)
        %       lower    - the lower boundary of the layer (from pia)
        %       instance - if the probe goes in and out of a VISp layber (can
        %                  happen at VISp edges) returns which instance of this VISp
        %                  layer we have come across (working downward
        %                   from the pia)
        
            boundaries = obj.VISp_boundaries();
            boundaries.upper = obj.from_tip_to_from_pia(boundaries.upper);
            boundaries.lower = obj.from_tip_to_from_pia(boundaries.lower);
        end
        
        
        
        function boundaries = VISp_boundaries_from_pia_flat(obj)
        %%VISp_boundaries_from_pia_flat Return vector with boundary distances for the VISp layers from pia (in um)
        %
        %   BOUNDARIES = VISp_boundaries_from_pia_flat()
        %   returns a vector, BOUNDARIES, containing the boundaries of each
        %   layer, starting with the upper boundary of VISp1 and ending with
        %   the lower boundary of VISp6b.
        
            b = obj.VISp_boundaries_from_pia();
            boundaries = [b.upper; b.lower(end)];
        end
        
        
        
        function [relative_depth, layer] = VISp_layer_relative_depth(obj, from_tip)
        %%VISp_layer_relative_depth Return relative depth in a layer, for a point given from probe tip (in um)
        %
        %   [RELATIVE_DEPTH, LAYER] = VISp_layer_relative_depth(FROM_TIP)
        %   get relative depth of a point in a VISp layer. FROM_TIP is a
        %   distance from the probe tip in um. RELATIVE_DEPTH is a value
        %   between 0 and 1 and specifies the relative depth in the layer
        %   from the *upper* boundary of the layer. LAYER is the string
        %   acronym of the layer. If point specified is not in VISp returns
        %   returned values are RELATIVE_DEPTH = NaN, LAYER = ''. 
        
            boundaries = obj.VISp_boundaries();
            
            idx = find(from_tip < boundaries.upper & from_tip >= boundaries.lower);
            
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