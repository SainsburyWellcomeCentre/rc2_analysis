classdef RasterDisplaySection < handle
    
    properties
        h_ax = gobjects(0);
    end
    
    methods
        
        function obj = RasterDisplaySection(h_fig)
            
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
            