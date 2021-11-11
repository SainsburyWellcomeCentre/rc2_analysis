classdef TuningCurvePlot < RC2Axis
    
    properties
        
        h_line
        h_dots
        h_errorbars
        h_fit
        h_txt
        
        h_line_shuff
        h_dots_shuff
        h_errorbars_shuff
        h_fit_shuff
        
        shuff
        
        error_type = 'sem'
        print_stats = true
        dot_size = 10
        line_width = 1
        
        xmin
        xmax
        ymin
        ymax
        
        n_shuffs_to_plot = 4
        main_col = [0, 0, 0]
    end
    
    
    methods
        
        function obj = TuningCurvePlot(h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        function plot(obj, tuning)
        %%fr, sd, n, x, shuff, stat_fr, stat_sd, stat_n, p_signrank, h_ax
            
            obj.shuff = tuning.shuffled;
            fr      = nanmean(tuning.tuning, 2);
            sd      = nanstd(tuning.tuning, [], 2);
            n       = sum(~isnan(tuning.tuning), 2);
            x       = tuning.bin_centers;
        
            for i = 1 : obj.n_shuffs_to_plot
                f = [min(x), max(x)]*tuning.shuffled.beta_shuff(i, 1) + tuning.shuffled.beta_shuff(i, 2);
                obj.h_fit_shuff(i) = line(obj.h_ax, [min(x), max(x)], f, 'color', [0.8, 0.8, 0.8], 'linewidth', obj.line_width);
            end
            
            stat_fr = nanmean(tuning.stationary_fr);
            stat_sd = nanstd(tuning.stationary_fr);
            stat_n  = sum(~isnan(tuning.stationary_fr));
            
            obj.h_dots = scatter(obj.h_ax, x, fr, obj.dot_size, obj.main_col, 'fill');
            for i = 1 : length(fr)
                if strcmp(obj.error_type, 'std')
                    y = fr(i) + sd(i) * [-1, 1];
                else
                    y = fr(i) + (sd(i)/sqrt(n(i))) * [-1, 1];
                end
                
                obj.h_errorbars(i) = line(obj.h_ax, x([i, i]), y, 'color', obj.main_col, 'linewidth', obj.line_width);
            end
            
            f = [min(x), max(x)]*tuning.shuffled.beta(1) + tuning.shuffled.beta(2);
            obj.h_fit = line(obj.h_ax, [min(x), max(x)], f, 'color', obj.main_col, 'linewidth', obj.line_width);
            
            
            scatter(obj.h_ax, 0, stat_fr, obj.dot_size, [0, 0, 0], 'fill')
%             if strcmp(obj.error_type, 'std')
%                 line(obj.h_ax, [0, 0], stat_fr + stat_sd * [-1, 1], 'color', [0.5, 0.5, 0.5])
%             else
%                 line(obj.h_ax, [0, 0], stat_fr + (stat_sd/sqrt(stat_n)) * [-1, 1], 'color', [0.5, 0.5, 0.5])
%             end
            
            if obj.print_stats
                
                str = sprintf('slope = %.2f\n', tuning.shuffled.beta(1));
                str = [str, sprintf('r = %.2f\n', tuning.shuffled.r)];
                str = [str, sprintf('R^2 = %.2f\n', tuning.shuffled.rsq)];
                
                if tuning.shuffled.p < 0.05
                    str = [str, sprintf('p_{shuffled}=%.2e\n', tuning.shuffled.p)];
                    str = [str, sprintf('p_{fitted}=%.2e', tuning.shuffled.p_lm)];
                    col = 'r';
                else
                    str = [str, sprintf('p_{shuffled}=%.2f\n', tuning.shuffled.p)];
                    str = [str, sprintf('p_{fitted}=%.2f', tuning.shuffled.p_lm)];
                    col = 'b';
                end
                
                if tuning.shuffled.beta(1) >= 0 && tuning.shuffled.p < 0.05
                    str = [str, ', +ve'];
                elseif tuning.shuffled.beta(1) < 0 && tuning.shuffled.p < 0.05
                    str = [str, ', -ve'];
                end
                
                
                obj.h_txt = text(obj.h_ax, obj.xmax, obj.ymin, ...
                    str, 'color', col, ...
                    'verticalalignment', 'bottom', 'horizontalalignment', 'right', ...
                    'fontsize', 6);
            end
            
            set(obj.h_ax, 'xlim', [-5, obj.xmax]);
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
            if obj.print_stats
                pos = get(obj.h_txt, 'position');
                set(obj.h_txt, 'position', [pos(1), val(1), 0]);
            end
        end
        
        
        
        function val = xlim(obj, val)
            
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'xlim');
            
            if obj.print_stats
                pos = get(obj.h_txt, 'position');
                set(obj.h_txt, 'position', [val(2), pos(2), 0]);
            end
            
            for i = 1 : obj.n_shuffs_to_plot
                
                f = [min(val), max(val)]*obj.shuff.beta_shuff(i, 1) + obj.shuff.beta_shuff(i, 2);
                set(obj.h_fit_shuff(i), 'xdata', [min(val), max(val)], 'ydata', f);
            end
            
            f = [min(val), max(val)]*obj.shuff.beta(1) + obj.shuff.beta(2);
            set(obj.h_fit, 'xdata', [min(val), max(val)], 'ydata', f);
        end
    end
end