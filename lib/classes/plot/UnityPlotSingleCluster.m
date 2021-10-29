classdef UnityPlotSingleCluster < UnityPlot
    
    properties
        
        p
        change_direction
        marker_style = 'o'
        h_dot
        
        dot_size = 10
        
    end
    
    methods
        
        function obj = UnityPlotSingleCluster(x, y, p, change_direction, h_ax)
            
            VariableDefault('h_ax', []);
            
            len = [length(x), length(y)];
            assert(length(unique(len)) == 1, 'inputs do not have equal lengths');
            
            obj = obj@UnityPlot(x, y, h_ax);
            
            obj.p = p;
            obj.change_direction = change_direction;
        end
        
        
        function plot(obj)
            
            n_trials = length(obj.x);
            
            if sum(isnan(obj.x)) > 0 || sum(isnan(obj.y)) > 0
                warning('There are %i x and %i y NaN FR values', sum(isnan(obj.x)), sum(isnan(obj.x)));
            end
            
            obj.h_dot = cell(n_trials, 1);
            
            scatter(obj.h_ax, obj.x, obj.y, obj.dot_size, [0.5, 0.5, 0.5], 'fill');
            
            obj.add_median_dot();
            obj.sync_axes();
            obj.add_unity_line();
            obj.print_stats();
        end
        
        
        function sync_axes(obj)
        %%ensure x and y axis have same range
            m = obj.min;
            M = obj.max;
            obj.xlim([m, M]);
            obj.ylim([m, M]);
        end
        
        
        
        function add_unity_line(obj)
        %%add unity line to plot
            axis_limits = obj.xlim();
            if ~isempty(obj.h_line)
                delete(obj.h_line);
            end
            obj.h_line = line(obj.h_ax, axis_limits, axis_limits, 'linestyle', '--', 'color', 'k');
        end
        
        
        
        function add_median_dot(obj)
        %%adds the median dot to the figure
            x_med = nanmedian(obj.x);
            y_med = nanmedian(obj.y);
            x_iqr = prctile(obj.x, [25, 75]);
            y_iqr = prctile(obj.y, [25, 75]);
            col = obj.get_p_colour();
            
            scatter(obj.h_ax, x_med, y_med, obj.dot_size, col, 'fill')
            line(obj.h_ax, x_iqr, y_med([1, 1]), 'color', col);
            line(obj.h_ax, x_med([1, 1]), y_iqr, 'color', col);
        end
        
        
        
        function col = get_p_colour(obj)
        %%get the colour of the median dot from p-value
        
            if obj.p < 0.05 && obj.change_direction == -1
                col = [30,144,255]/255;
            elseif obj.p < 0.05 && obj.change_direction == 1
                col = 'r';
            elseif obj.p < 0.05 && obj.change_direction == 0
                col = 'm';
            else
                col = 'k';
            end
        end
        
        
        function print_stats(obj)
       %% print the statistics
       
            m = obj.min;
            M = obj.max;
            
            if obj.p < 0.01
                format_str = 'p = %.2e \nx median = %.2f (n=%i) \ny median = %.2f (n=%i)\n# nan (x=%i,y=%i)';
            else
                format_str = 'p = %.2f \nx median = %.2f (n=%i) \ny median = %.2f (n=%i)\n# nan (x=%i,y=%i)';
            end
            
            txt_str = sprintf(format_str, ...
                    obj.p, ...
                    nanmedian(obj.x), ...
                    sum(~isnan(obj.x)), ...
                    nanmedian(obj.y), ...
                    sum(~isnan(obj.y)), ...
                    sum(isnan(obj.x)), ...
                    sum(isnan(obj.y)));
            
            obj.h_txt = text(obj.h_ax, M, m, txt_str, 'verticalalignment', 'bottom', 'horizontalalignment', 'right', ...
                'fontsize', 4);
        end
    end
end
