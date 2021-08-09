classdef UnityPlotPopulation < UnityPlot
    
    properties
        
        p
        y_gt_x
        marker_style = 'o'
        h_dot
    end
    
    methods
        
        function obj = UnityPlotPopulation(x, y, p, y_gt_x, h_ax)
            
            len = [length(x), length(y), length(p), length(y_gt_x)];
            assert(length(unique(len)) == 1, 'inputs do not have equal lengths');
            
            obj = obj@UnityPlot(x, y, h_ax);
            
            obj.p = p;
            obj.y_gt_x = y_gt_x;
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
                    obj.h_dot{clust_i} = scatter(obj.x(clust_i), obj.y(clust_i), scatterball_size(0.8), obj.black, 'fill', obj.marker_style);
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
                    obj.h_dot{clust_i} = scatter(obj.x(clust_i), obj.y(clust_i), scatterball_size(1.2), col, obj.marker_style);
                end
            end
            
            % format
            m = min([get(obj.h_ax, 'xlim'), get(obj.h_ax, 'ylim')]);
            M = max([get(obj.h_ax, 'xlim'), get(obj.h_ax, 'ylim')]);
            
            set(obj.h_ax, 'xlim', [m, M], 'ylim', [m, M]);
            obj.h_line = line(obj.h_ax, [m, M], [m, M], 'linestyle', '--', 'color', 'k');
            
            x_med = nanmedian(obj.x);
            y_med = nanmedian(obj.y);
            x_iqr = prctile(obj.x, [25, 75]);
            y_iqr = prctile(obj.y, [25, 75]);
            p_median = signrank(obj.x, obj.y);
            
            scatter(obj.h_ax, x_med, y_med, scatterball_size(1.3), 'g', 'fill');
            line(obj.h_ax, x_iqr, y_med([1, 1]), 'color', 'g');
            line(obj.h_ax, x_med([1, 1]), y_iqr, 'color', 'g');
            
            n_sig_increase = sum(obj.p < 0.05 & obj.y_gt_x);
            n_sig_decrease = sum(obj.p < 0.05 & ~obj.y_gt_x);
            n_sig_no_change = sum(obj.p >= 0.05);
            
            
            if p_median < 0.01
                format_str = 'p = %.2e \nred = %i (%i%%) \nblue = %i (%i%%)\nblack = %i (%i%%)\ntotal = %i';
            else
                format_str = 'p = %.2f \nred = %i (%i%%) \nblue = %i (%i%%)\nblack = %i (%i%%)\ntotal = %i';
            end
            
            txt_str = sprintf(format_str, ...
                    p_median, ...
                    n_sig_increase, ...
                    round(100*n_sig_increase/n_clust), ...
                    n_sig_decrease, ...
                    round(100*n_sig_decrease/n_clust), ...
                    n_sig_no_change, ...
                    round(100*n_sig_no_change/n_clust), ...
                    n_clust);
            
            obj.h_txt = text(obj.h_ax, M, m, txt_str, 'verticalalignment', 'bottom', 'horizontalalignment', 'right', ...
                'fontsize', 6);
            
        end
        
        
        function highlight(obj, idx)
            
            assert(length(idx) == length(obj.x))
            scatter(obj.h_ax, obj.x(idx), obj.y(idx), scatterball_size(1.3), 'm');
        end
    end
end