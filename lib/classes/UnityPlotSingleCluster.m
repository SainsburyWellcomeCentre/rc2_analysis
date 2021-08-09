classdef UnityPlotSingleCluster < UnityPlot
    
    properties
        
    end
    
    methods
        
        function obj = UnityPlotSingleCluster(x, y, h_ax)
            
            obj = obj@UnityPlot(x, y, h_ax);
            
            scatter(x, y, 10, [0.5, 0.5, 0.5], 'fill');
            
            % format
            m = min([get(obj.h_ax, 'xlim'), get(obj.h_ax, 'ylim')]);
            M = max([get(obj.h_ax, 'xlim'), get(obj.h_ax, 'ylim')]);
            
            set(obj.h_ax, 'xlim', [m, M], 'ylim', [m, M]);
            obj.h_line = line(obj.h_ax, [m, M], [m, M], 'linestyle', '--', 'color', 'k');
            
            x_med = nanmedian(x);
            y_med = nanmedian(y);
            x_iqr = prctile(x, [25, 75]);
            y_iqr = prctile(y, [25, 75]);
            p_median = signrank(x, y);
            
            if p_median < 0.05 && x_med > y_med
                col = [30,144,255]/255;
            elseif p_median < 0.05 && x_med < y_med
                col = 'r';
            elseif p_median < 0.05 && x_med == y_med
                col = 'm';
            else
                col = 'k';
            end
            
            scatter(obj.h_ax, x_med, y_med, 20, col, 'fill');
            line(obj.h_ax, x_iqr, y_med([1, 1]), 'color', col);
            line(obj.h_ax, x_med([1, 1]), y_iqr, 'color', col);
            
            if p_median < 0.01
                txt_str = sprintf('p = %.2e', p_median);
            else
                txt_str = sprintf('p = %.2f', p_median);
            end
            
            obj.h_txt = text(obj.h_ax, M, m, txt_str, 'verticalalignment', 'bottom', 'horizontalalignment', 'right', ...
                'fontsize', 8);
            
        end
        
    end
    
end