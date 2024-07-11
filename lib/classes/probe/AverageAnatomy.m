classdef AverageAnatomy < handle
% AverageAnatomy Class for handling averaging the information across
% several Anatomy objects
%
%   AverageAnatomy Properties:
%       anatomy_array       - array of Anatomy objects
%       VISp_layers         - cell array of strings containing acronyms for VISp layers
%
%   AverageAnatomy Methods:
%       average_VISp_boundaries_from_pia_struct     - computes the average distance of the boundaries of each
%                                                     VISp layer when averaged across several anatomies (returns structure)
%       average_VISp_boundaries_from_pia            - computes the average distance of the boundaries of each VISp 
%                                                     layer when averaged across several anatomies (returns vector)
%       from_pia_using_relative_position            - computes depth of a point from the pia in the averaged anatomy 
%                                                     from its relative position of the point in its layer 
%
%   See also: Anatomy, ProbeTrack

    properties
        
        anatomy_array
    end
    
    properties (Constant = true)
        
        VISp_layers = {'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
    end
    
    
    
    methods
        
        function obj =  AverageAnatomy(anatomy_array)
        %%AverageAnatomy
        %
        %   AverageAnatomy(ANATOMY_ARRAY) class for combining a set of
        %   anatomies to get averaged boundaries.
        
            obj.anatomy_array = [anatomy_array{:}];
        end
        
        
        
        function boundaries = average_VISp_boundaries_from_pia_struct(obj)
        %%average_VISp_boundaries_from_pia_struct Computes the average
        %%distance of the boundaries of each VISp layer when averaged
        %%across several anatomies
        %
        %   BOUNDARIES = average_VISp_boundaries_from_pia_struct() returns
        %   an averaged boundary structure, BOUNDARIES, which simply
        %   averages the upper and lower boundaries of each layer across
        %   all the Anatomy objects in `anatomy_array`.
        %
        %   See also:   Anatomy.VISp_boundaries_from_pia
            
            for ii = 1 : length(obj.anatomy_array)
                
                b(ii) = obj.anatomy_array{ii}.VISp_boundaries_from_pia;
                
                % make sure that none of the anatomies are 'strange' (i.e.
                % have multiple of the same layers)
%                 assert(isequal(b(ii).region, obj.VISp_layers'), 'anatomy likely has multiple versions of the same layer');
            end
            
            for ii = 1 : length(obj.anatomy_array)
                if size(b(ii).upper) < length(obj.VISp_layers)
                    b(ii) = obj.correct_for_missing_layers(b(ii));
                end
            end
            
            boundaries.region = obj.VISp_layers;
            boundaries.upper = mean(cat(2, b(:).upper), 2);
            boundaries.lower = mean(cat(2, b(:).lower), 2);
        end
        
            
        
        
        function boundary_position = average_VISp_boundaries_from_pia(obj)
        %%average_VISp_boundaries_from_pia Computes the average
        %%distance of the boundaries of each VISp layer when averaged
        %%across several anatomies
        %
        %   BOUNDARIES = average_VISp_boundaries_from_pia() returns
        %   an averaged vector, BOUNDARIES, which simply
        %   averages the upper and lower boundaries of each layer across
        %   all the Anatomy objects in `anatomy_array`.  Vector starts with
        %   the upper boundary of VISp1 and ends with the lower boundary of
        %   VISp6b. 
        %
        %   See also:   Anatomy.VISp_boundaries_from_pia_flat
        
            boundaries = obj.average_VISp_boundaries_from_pia_struct();
            boundary_position = [boundaries.upper(:); boundaries.lower(end)];
        end
        
        

        function cluster_position = from_pia_using_relative_position(obj, relative_position, layer)
        %%from_pia_using_relative_position Computes depth of a point from
        %%the pia in  the averaged anatomy from its relative position in
        %%its layer 
        %
        %   DEPTH = from_pia_using_relative_position(RELATIVE_POSITION, LAYER)
        %   takes a vector of relative positional points in RELATIVE
        %   POSITION (each a value between 0 and 1, representing the
        %   distance of a point from the upper boundary of a VISp layer)
        %   and a cell array of strings specifying the VISp layer of each
        %   point in LAYER. From this, it computes the position of each
        %   point in the averaged VISp "cortex" from the pia (in um),
        %   returned in DEPTH. RELATIVE_POSITON is a #cluster x 1 vector,
        %   LAYER is a #cluster x 1 cell array of strings, and DEPTH will
        %   be a #cluster x 1 vector.
        %   
        %   See also:   Anatomy.VISp_layer_relative_depth
        
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
    
    
    methods (Static = true)
        function boundary = correct_for_missing_layers(boundary)
            dummy_layer_6a_length = 235;
            dummy_layer_6b_length = 80;
            
            if size(boundary.region) == [4 1]
                % add layer 6a info
                boundary.region = [boundary.region; 'VISp6a'];
                boundary.upper = [boundary.upper; boundary.lower(end)];
                boundary.lower = [boundary.lower; boundary.lower(end) + dummy_layer_6a_length];
                boundary.instance = [boundary.instance; 0];
            end

            if size(boundary.region) == [5 1]
                % add layer 6b info
                boundary.region = [boundary.region; 'VISp6b'];
                boundary.upper = [boundary.upper; boundary.lower(end)];
                boundary.lower = [boundary.lower; boundary.lower(end) + dummy_layer_6b_length];
                boundary.instance = [boundary.instance; 0];
            end
        end 
    end
end