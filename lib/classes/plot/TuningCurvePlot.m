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
        
        fr
        sd
        n
        x
        p_anova
        beta
        shuff
        p_signrank
        stat_fr
        stat_sd
        stat_n
        
        error_type = 'sem'
        
        xmin
        xmax
        ymin
        ymax
        
        n_shuffs_to_plot = 4
        main_col
    end
    
    
    methods
        
        function obj = TuningCurvePlot(fr, sd, n, x, shuff, stat_fr, stat_sd, stat_n, p_signrank, h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
            
            if ~exist('fr', 'var')
                return
            end
            
            obj.x = x;
            obj.fr = fr;
            obj.sd = sd;
            obj.n = n;
%             obj.p_anova = p_anova;
%             obj.beta = beta;
            obj.shuff = shuff;
            obj.p_signrank = p_signrank;
            obj.stat_fr = stat_fr;
            obj.stat_sd = stat_sd;
            obj.stat_n = stat_n;
            
            multicol = lines(4);
            
            if p_signrank < 0.05 && ...
                    mean(obj.fr) >= obj.stat_fr %  p_anova >= 0.05 && ...
                % tonic increase
                main_col = multicol(2, :);
            elseif p_signrank < 0.05 && ...
                    mean(obj.fr) < obj.stat_fr  % p_anova >= 0.05 && ...
                % tonic decrease
                main_col = multicol(1, :);
            else
                main_col = 'k';
            end
            
            for i = 1 : obj.n_shuffs_to_plot
                
%                 obj.h_line_shuff(i) = line(obj.x, obj.shuff.shuff_tuning(:, i), 'color', [0.8, 0.8, 0.8]);
% %                 obj.h_dots_shuff(i) = scatter(obj.x, obj.shuff.shuff_tuning(:, i), [], [0.8, 0.8, 0.8], 'fill');
%                 
%                 for j = 1 : size(obj.shuff.shuff_tuning, 1)
%                     if strcmp(obj.error_type, 'std')     
%                         y = obj.shuff.shuff_tuning(j, i) + obj.shuff.shuff_sd(j, i) * [-1, 1];
%                     else
%                         y = obj.shuff.shuff_tuning(j, i) + (obj.shuff.shuff_sd(j, i) / sqrt(obj.shuff.shuff_n(j, i))) * [-1, 1];
%                     end
%                     
%                     obj.h_errorbars_shuff(i, j) = line(obj.x([j, j]), y, 'color', [0.8, 0.8, 0.8]);
%                 end
                
                f = [min(obj.x), max(obj.x)]*obj.shuff.beta_shuff(i, 1) + obj.shuff.beta_shuff(i, 2);
                obj.h_fit_shuff(i) = line([min(obj.x), max(obj.x)], f, 'color', [0.8, 0.8, 0.8]);
            end
            
            
            obj.h_line = line(obj.x, obj.fr, 'color', main_col);
            obj.h_dots = scatter(obj.x, obj.fr, [], main_col, 'fill');
            
            for i = 1 : length(obj.fr)
                if strcmp(obj.error_type, 'std')
                    
                    y = obj.fr(i) + obj.sd(i) * [-1, 1];
                else
                    y = obj.fr(i) + (obj.sd(i)/sqrt(obj.n(i))) * [-1, 1];
                end
                
                obj.h_errorbars(i) = line(obj.x([i, i]), y, 'color', main_col);
            end
            
            f = [min(obj.x), max(obj.x)]*obj.shuff.beta(1) + obj.shuff.beta(2);
            obj.h_fit = line([min(obj.x), max(obj.x)], f, 'color', main_col);
            
            scatter(0, obj.stat_fr, [], [0.5, 0.5, 0.5], 'fill')
            if strcmp(obj.error_type, 'std')
                line([0, 0], obj.stat_fr + obj.stat_sd * [-1, 1], 'color', [0.5, 0.5, 0.5])
            else
                line([0, 0], obj.stat_fr + (obj.stat_sd/sqrt(obj.stat_n)) * [-1, 1], 'color', [0.5, 0.5, 0.5])
            end
            
            str = sprintf('slope = %.2f\n', obj.shuff.beta(1));
            str = [str, sprintf('r = %.2f\n', obj.shuff.r)];
            str = [str, sprintf('R^2 = %.2f\n', obj.shuff.rsq)];
            
            if obj.shuff.p < 0.05
                str = [str, sprintf('p_{shuffled}=%.2e\n', obj.shuff.p)];
                str = [str, sprintf('p_{fitted}=%.2e', obj.shuff.p_lm)];
                col = 'r';
            else
                str = [str, sprintf('p_{shuffled}=%.2f\n', obj.shuff.p)];
                str = [str, sprintf('p_{fitted}=%.2f', obj.shuff.p_lm)];
                col = 'b';
            end
            
            if obj.shuff.beta(1) >= 0 && obj.shuff.p < 0.05
                str = [str, ', +ve'];
            elseif obj.shuff.beta(1) < 0 && obj.shuff.p < 0.05
                str = [str, ', -ve'];
            end
                
            
            obj.h_txt = text(obj.h_ax, obj.xmin, obj.ymax, ...
                str, 'color', col, ...
                'verticalalignment', 'bottom', 'horizontalalignment', 'right', ...
                'fontsize', 6);
            
            set(obj.h_ax, 'xlim', [-5, obj.xmax]);
            
            obj.main_col = main_col;
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