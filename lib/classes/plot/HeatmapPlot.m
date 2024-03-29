classdef HeatmapPlot < RC2Axis
% HeatmapPlot Class for plotting heatmaps
%
%  HeatmapPlot Properties:
%    h_im       - handle to the heatmap image
%    h_ax_col   - handle to the axis for the colourbar
%    h_im_col   - handle to the image in the colourbar
%    fr_limits  - colour limits of the heatmap and colourbar
%
%  HeatmapPlot Methods:
%    plot       - plots the heatmap
%
%   Creates a heatmap image along with a colorbar in the top right

    properties (SetAccess = private)
        
        h_im
        h_ax_col
        h_im_col
    end
    
    properties
        
        fr_limits = [-8, 8]
    end
    
    
    
    methods
        
        function obj = HeatmapPlot(h_ax)
        %%HeatmapPlot
        %
        %   HeatmapPlot(AXIS_HANDLE) sets up an object for handling
        %   heatmaps. AXIS_HANDLE is optional, if supplied it should be a
        %   handle to an axis object. Otherwise, an axis will be created.
        
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        
        function plot(obj, t, heatmap)
        %%plot Plots the heatmap
        %
        %   plot(TIME, HEATMAP) plots the 2D matrix HEATMAP plotted against TIME on
        %   the x-axis. TIME should have same length as the second
        %   dimension of HEATMAP.
        
            load('r2bcmap.mat', 'map');
            
            n_clusters = size(heatmap, 1);
            
            % plot the data
            obj.h_im  = imagesc(obj.h_ax, heatmap);
            
            % 0 time line
            line(obj.h_ax, [0, 0], [0.5, n_clusters+0.5], 'linestyle', '--', 'color', 'k');
            
            colormap(obj.h_ax, map);
            
            set(obj.h_im, 'xdata', t);
            
            set(obj.h_ax, 'clim', obj.fr_limits, ...
                          'ytick', [1, n_clusters], ...
                          'xlim', t([1, end]), ...
                          'ylim', [0.5, n_clusters+0.5], ...
                          'xcolor', 'none', ...
                          'ycolor', 'none', ...
                          'yticklabel', [n_clusters, 1]);
                      
            text(obj.h_ax, t(1)-(t(2)-t(1))*0.1, n_clusters, sprintf('%i', 1), ...
                'horizontalalignment', 'right', ...
                'verticalalignment', 'middle', ...
                'fontsize', 8);
            text(obj.h_ax, t(1)-(t(2)-t(1))*0.1, 1, sprintf('%i', n_clusters), ...
                'horizontalalignment', 'right', ...
                'verticalalignment', 'middle', ...
                'fontsize', 8);
            text(obj.h_ax, t(1)-(t(2)-t(1))*0.2, (n_clusters + 1) / 2, 'cell #', ...
                'horizontalalignment', 'center', ...
                'verticalalignment', 'bottom', ...
                'fontsize', 8, ...
                'rotation', 90);
            
            % colorbar
            axis_position = get(obj.h_ax, 'position');
            axis_width = axis_position(3);
            axis_height = axis_position(4);
            
            colorbar_width = 0.05 * axis_width;
            colorbar_height = 0.35 * axis_height;
            
            colorbar_position = [(axis_position(1) + axis_width) + colorbar_width, ...
                                 (axis_position(2) + axis_height) - colorbar_height, ...
                                 colorbar_width, colorbar_height];
            
            obj.h_ax_col = axes();
            
            set(obj.h_ax_col, ...
                'color',                    'w', ...
                'xcolor',                   'k', ...
                'ycolor',                   'k', ...
                'units',                    'normalized', ...
                'positionconstraint',       'innerposition', ...
                'innerposition',            colorbar_position, ...
                'plotboxaspectratiomode',   'auto');
            hold on;
            
            % plot the colormap
            cmap_size = [3, 64];
            dummy_map = repmat(linspace(0, 1, cmap_size(2))', 1, cmap_size(1));
            obj.h_im_col = imagesc(obj.h_ax_col, dummy_map, [0, 1]);
            colormap(obj.h_ax_col, map);
            
            set(obj.h_ax_col, 'xlim', [0.5, cmap_size(1)+0.5], 'ylim', [0.5, cmap_size(2)+0.5])
            
            axis off
            
            text(obj.h_ax_col, cmap_size(1)+0.5, -1, sprintf('<%i', obj.fr_limits(1)), ...
                'color', 'k', ...
                'horizontalalignment', 'center', ...
                'verticalalignment', 'top', ...
                'fontsize', 5);
            text(obj.h_ax_col, cmap_size(1)+0.5, cmap_size(2)+1, sprintf('<%i', obj.fr_limits(1)), ...
                'color', 'k', ...
                'horizontalalignment', 'center', ...
                'verticalalignment', 'bottom', ...
                'fontsize', 5);
            text(obj.h_ax_col, cmap_size(1)+0.6, (cmap_size(2)+1)/2, '\DeltaFR (Hz)', ...
                'color', 'k', ...
                'horizontalalignment', 'center', ...
                'verticalalignment', 'bottom', ...
                'fontsize', 5, ...
                'rotation', 270);
        end
    end
end
