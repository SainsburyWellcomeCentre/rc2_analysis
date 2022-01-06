classdef RasterDisplaySection < handle
% RasterDisplaySection Class for controlling a section of a
% RasterDisplayFigure.
%
%   RasterDisplaySection Properties:
%       h_ax            - array of axes in this section
%
%   RasterDisplaySection Methods:
%       set_position    - Sets the position of the axes in the figure
%       
%   TDDO:   1. make three axes more general

    properties
        
        h_ax = gobjects(0);
    end
    
    
    
    methods
        
        function obj = RasterDisplaySection(h_fig)
        %%RasterDisplaySection
        %
        %   RasterDisplaySection(FIGURE_HANDLE) associates a section with a
        %   figure referenced by FIGURE_HANDLE. Currently a section is
        %   composed of three axes which display associated data for raster
        %   plots.
        
            if ~exist('h_fig', 'var')
                return
            end
            
            figure(h_fig);
            
            for i = 1 : 3
                obj.h_ax(i) = axes('units', 'centimeters');
                axis normal
            end
        end
        
        
        
        function set_position(obj, pos)
        %%set_position Sets the position of the axes in the section
        %
        %   set_position(POSITION) puts the axes of the section in the box
        %   defined by POSITION. POSITION is of the form [from left,
        %   from_right, width, height] and defines an area on the figure.
        
            x_spacing = pos(3)/10;
            y_spacing = pos(4)/10;
            
            ax_width = pos(3) - 2*x_spacing;
            ax_height = (pos(4) - 4*y_spacing)/3;
            
            for i = 1 : 3
                obj.h_ax(i).Position = [pos(1) + x_spacing, ...
                    pos(2) + (4-i)*y_spacing + (3-i)*ax_height, ...
                    ax_width, ...
                    ax_height];
            end
        end
    end
end
            