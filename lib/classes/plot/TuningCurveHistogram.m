classdef TuningCurveHistogram < RC2Axis
    
    properties
        
        h_bar
        
        xmin
        xmax
        ymin
        ymax
        
        dot_col = [1, 0, 0]
        dot_size = 10
        plot_labels = true;
    end
    
    
    methods
        
        function obj = TuningCurveHistogram(h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        function plot(obj, tuning)
            
            [n, c] = histcounts(tuning.shuffled.r_shuff, 50);
            obj.h_bar = barh(obj.h_ax, (c(1:end-1)+c(2:end))/2, n, 'facecolor', [0.6, 0.6, 0.6], 'edgecolor', 'none', 'barwidth', 1);
            set(obj.h_ax, 'plotboxaspectratio', [1, 3, 1])
            ylim([-0.5, 0.5]);
            scatter(obj.h_ax, 0, tuning.shuffled.r, obj.dot_size, obj.dot_col, 'fill');
            if obj.plot_labels
                xlabel('count');
                ylabel('r');
            end
        end
        
        
        
        function val = get.xmin(obj)
            
            val = min(get(obj.h_ax, 'xlim')); %#ok<*CPROP>
        end
        
        
        
        function val = get.xmax(obj)
            
            val = max(get(obj.h_ax, 'xlim'));
        end
        
        
        
        function val = get.ymin(obj)
            
            val = min(get(obj.h_ax, 'ylim')); %#ok<*CPROP>
        end
        
        
        
        function val = get.ymax(obj)
            
            val = max(get(obj.h_ax, 'ylim'));
        end
        
        
        
        function color(obj, val)
            
            set(obj.h_line, 'color', val);
            for i = 1 : length(obj.h_errorbars)
                set(obj.h_errorbars(i), 'color', val);
            end
        end
        
        
        
        function val = ylim(obj, val)
            
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'ylim');
            pos = get(obj.h_txt, 'position');
            set(obj.h_txt, 'position', [pos(1), val(1), 0]);
        end
        
        
        
        function val = xlim(obj, val)
            
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'xlim');
            pos = get(obj.h_txt, 'position');
            set(obj.h_txt, 'position', [val(2), pos(2), 0]);
            
            for i = 1 : obj.n_shuffs_to_plot
                
                f = [min(val), max(val)]*obj.shuff.beta_shuff(i, 1) + obj.shuff.beta_shuff(i, 2);
                set(obj.h_fit_shuff(i), 'xdata', [min(val), max(val)], 'ydata', f);
            end
            
            f = [min(val), max(val)]*obj.shuff.beta(1) + obj.shuff.beta(2);
            set(obj.h_fit, 'xdata', [min(val), max(val)], 'ydata', f);
        end
    end
end