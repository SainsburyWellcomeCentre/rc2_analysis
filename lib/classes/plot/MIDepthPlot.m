classdef MIDepthPlot < RC2Axis
% MIDepthPlot Class for plotting the depth of a cluster against 
%
%  MIDepthPlot Properties:
%       marker_style    - style of dot to plot for each cluster
%   
%  MIDepthPlot Methods:
%       plot            - make the plot
%       print_layers    - print anatomical region labels at the side
%
%   TODO:   1. this class has too many required arguments.
%           2. make p-value threshold an adjustable property

    properties
        
        marker_style = 'o'
    end
    
    properties (SetAccess = private)
        
        h_lines
        h_txt
        h_dot
        
        mi
        direction
        depths
        boundaries
        regions
        p
    end
    
    
    methods
        
        function obj = MIDepthPlot(mi, p, direction, depths, boundaries, regions, h_ax)
        %%MIDepthPlot
        %
        %   MIDepthPlot(MODULATION_INDEX, P_VALUE, DIRECTION, DEPTHS, BOUNDARIES, REGIONS, AXIS_HANDLE)
        %   creates an object in preparation for plotting the depth vs.
        %   modulation index of each cluster.
        %
        %   Args:
        %       MODULATION_INDEX - # clusters x 1 array containing the modulation index of each cluster
        %       P_VALUE          - # clusters x 1 array containing a p-value for each cluster
        %       DIRECTION        - # clusters x 1 array containing 1, 0 or
        %                          -1 to indicate whether a response was +ve, no change or negative
        %       DEPTHS           - # clusters x 1 array containing the depth of each cluster
        %       BOUNDARIES       - (# regions + 1) x 1 array containing the depth of the anatomical region boundaries
        %       REGIONS          - # regions x 1 cell array of strings with region labels
        
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
            
            obj.mi = mi;
            obj.direction = direction;
            obj.depths = depths;
            obj.boundaries = boundaries;
            obj.regions = regions;
            obj.p = p;
        end
        
        
        
        function plot(obj)
        %%plot Make the depth vs. MI plot
        %
        %   plot() scatters the depth (`depth`) of each cluster against its modulation
        %   index (`mi`). A cluster is grey if the p-value (`p`) is > 0.05,
        %   red if p-value is < 0.05 and `direction` is 1, and blue if
        %   p-value is < 0.05 and `direction` is -1.
        %
        %   Anatomical boundaries are also plotted. To print region labels use
        %   `print_layers` method.
        
            n_clust = length(obj.p);
            
            % warn user if there are NaN p-values
            if sum(isnan(obj.p)) > 0
                warning('There are %i NaN p-values', sum(isnan(obj.p)))
            end
            
            obj.h_dot = cell(n_clust, 1);
            
            % loop through clusters and put non-significant dots in the
            % background
            for clust_i = 1 : n_clust
                if obj.p(clust_i) >= 0.05
                    obj.h_dot{clust_i} = scatter(obj.mi(clust_i), obj.depths(clust_i), scatterball_size(0.8), obj.black, 'fill', obj.marker_style);
                end
            end
            
            % loop through clusters again and plot the significant ones
            for clust_i = 1 : n_clust
                if obj.p(clust_i) < 0.05
                    if obj.direction(clust_i) == 1
                        col = obj.red;
                    elseif obj.direction(clust_i) == -1
                        col = obj.blue;
                    else
                        error('p-value significant but no direction');
                    end
                    obj.h_dot{clust_i} = scatter(obj.mi(clust_i), obj.depths(clust_i), scatterball_size(1.2), col, obj.marker_style);
                end
            end
            
            
            
            for b_i = 1 : length(obj.boundaries)
                line([-1, 1], obj.boundaries(b_i)*[1, 1], 'color', 'k', 'linestyle', '--');
            end
            
            set(obj.h_ax, 'xlim', [-1, 1], 'ylim', [obj.boundaries(1)-10, obj.boundaries(end)+10], ...
                'xtick', [-1, 0, 1], 'ydir', 'reverse');
            
        end
        
        
        
        function print_layers(obj, position)
        %%print_layers  Add anatomical region labels at the side
        %
        %   print_layers(POSITION) adds the labels in `regions` to the
        %   plot. POSITION is either 'right' (default) or 'left' and
        %   indicates which side of the axis to place the text.
        
            VariableDefault('position', 'right');
            
            % format
            for b_i = 1 : length(obj.boundaries)
    
                if b_i < length(obj.boundaries)-1
                    if strcmp(position, 'right')
                        text(1, sum(obj.boundaries(b_i) + obj.boundaries(b_i+1))/2, obj.regions{b_i}, ...
                            'horizontalalignment', 'left', 'verticalalignment', 'middle')
                    else
                        text(-1, sum(obj.boundaries(b_i) + obj.boundaries(b_i+1))/2, obj.regions{b_i}, ...
                            'horizontalalignment', 'right', 'verticalalignment', 'middle')
                    end
                end
            end
        end
    end
end