classdef PairedBoxPlot < RC2Axis
    
    properties
        
        h_txt
        
        data
        
        min
        max
    end
    
    
    
    methods
        
        function obj = PairedBoxPlot(x, y, h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);

            obj.data = [x(:), y(:)];
            
            pos = [2, 3];
            
            for i = 1 : size(obj.data, 1)
                line(pos, [obj.data(i, 1), obj.data(i, 2)], 'color', 'k', 'linewidth', 0.25);
            end
            for i = 1 : 2
                scatter(pos(i)*ones(size(obj.data, 1), 1), obj.data(:, i), scatterball_size(1), 'k', 'fill');
            end
            
            medians = nan(1, 2);
            box = cell(1, 2);
            outliers = cell(2, 2);
            whiskers = nan(2, 2);
            
            for i = 1 : 2
                medians(i) = median(obj.data(:, i));
                box{i} = prctile(obj.data(:, i), [25, 75]);
                for j = 1 : 2
                    if j == 1
                        outliers{i, j} = obj.data(:, i) < box{i}(j) - 1.5 * iqr(obj.data(:, i));
                        whiskers(i, j) = min(obj.data(~outliers{i, j}, i));
                    else
                        outliers{i, j} = obj.data(:, i) > box{i}(j) + 1.5 * iqr(obj.data(:, i));
                        whiskers(i, j) = max(obj.data(~outliers{i, j}, i));
                    end
                end
            end
            
            
            pos = [1, 4];
            width = 0.3;
            for i = 1 : 2
                line(pos(i) + [-width, width], medians([i, i]), 'color', 'k', 'linewidth', 0.5);
                h = patch(pos(i) + [-width, width, width, -width], box{i}([1, 1, 2, 2]), 'k');
                set(h, 'facecolor', 'none', 'edgecolor', 'k');
                
                for j = 1 : 2
                    
                    line(pos([i, i]), [whiskers(i, j), box{i}(j)], 'color', 'k', 'linewidth', 0.25);
                    line(pos(i) + [-0.05, 0.05], whiskers(i, j)*[1, 1], 'color', 'k', 'linewidth', 0.25);
                    scatter(pos(i)*ones(sum(outliers{i, j}), 1), obj.data(outliers{i, j}, i), scatterball_size(1), 'k');
                end
            end
            
            
            set(obj.h_ax, 'xlim', [0, 5], 'xtick', [1, 4]);
            
            % stats
            p = signrank(obj.data(:, 1), obj.data(:, 2));
            if p < 0.01
                txt_str = sprintf('p = %.2e', p);
            else
                txt_str = sprintf('p = %.2f', p);
            end
            
            obj.h_txt = ...
                text(obj.h_ax, 0, max(get(obj.h_ax, 'ylim')), txt_str, ...
                'verticalalignment', 'top', 'horizontalalignment', 'left', ...
                'fontsize', 8);
            
            if p < 0.05
                set(obj.h_txt, 'color', 'r');
            end
        end
        
        
        
        function val = get.min(obj)
            
            val = min(get(obj.h_ax, 'ylim')); %#ok<*CPROP>
        end
        
        
        
        function val = get.max(obj)
            
            val = max(get(obj.h_ax, 'ylim'));
        end
        
        
        
        function xticklabel(obj, str1, str2)
            
            set(obj.h_ax, 'xticklabel', {str1, str2});
        end
        
        
        
        function val = ylim(obj, val)
            
            val = obj.get_set_limits(val, 'ylim');
            set(obj.h_txt, 'position', [0, val(2), 0]);
        end
    end
end