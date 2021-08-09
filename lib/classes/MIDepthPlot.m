classdef MIDepthPlot < RC2Axis
    
    properties
        
        h_lines
        h_txt
        h_dot
        
        x
        y
        y_gt_x
        depths
        boundaries
        regions
        p
        mi
        
        marker_style = 'o'
        
    end
    
    
    methods
        
        function obj = MIDepthPlot(x, y, p, y_gt_x, depths, boundaries, regions, h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
            
            if ~exist('x', 'var')
                return
            end
            
            obj.x = x;
            obj.y = y;
            obj.y_gt_x = y_gt_x;
            obj.depths = depths;
            obj.boundaries = boundaries;
            obj.regions = regions;
            obj.p = p;
            
            obj.mi = (y - x) ./ (y + x);
        end
        
        
        
        function plot(obj)
            
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
                    if obj.y_gt_x(clust_i)
                        col = obj.red;
                    else
                        col = obj.blue;
                    end
                    obj.h_dot{clust_i} = scatter(obj.mi(clust_i), obj.depths(clust_i), scatterball_size(1.2), col, obj.marker_style);
                end
            end
            
            
            
            for b_i = 1 : length(obj.boundaries)
                line([-1, 1], obj.boundaries(b_i)*[1, 1], 'color', 'k', 'linestyle', '--');
            end
            
            set(obj.h_ax, 'xlim', [-1, 1], 'ylim', [obj.boundaries(end)-10, obj.boundaries(1)+10], ...
                'ytick', [], 'xtick', [-1, 0, 1]);
            
        end
        
        
        
        function print_layers(obj, position)
            
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