classdef ExtendedPairedPlot < RC2Axis
    
    properties
        
        h_txt
        
        data
        n_groups
        n_points
        
        min
        max
    end
    
    
    
    methods
        
        function obj = ExtendedPairedPlot(X, h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);

            obj.data = X;
            obj.n_groups = size(X, 2);
            obj.n_points = size(X, 1);
            
            for i = 1 : obj.n_groups
                scatter(i*ones(obj.n_points, 1), obj.data(:, i), scatterball_size(2), 'k', 'fill');
                if i < obj.n_groups
                    for j =  1 : obj.n_points
                        line([i, i+1], [obj.data(j, i), obj.data(j, i+1)], 'color', 'k', 'linewidth', 0.25);
                    end
                end
            end
            
            set(obj.h_ax, 'xlim', [0, obj.n_groups+1], 'plotboxaspectratio', [2, 1, 1], ...
                'xtick', 1:obj.n_groups);
            
            % stats
            p = signrank(obj.data(:, 2), obj.data(:, 4));
            if p < 0.01
                txt_str = sprintf('p = %.2e', p);
            else
                txt_str = sprintf('p = %.2f', p);
            end
            
            line(obj.h_ax, [2, 4], 0.8*max(get(obj.h_ax, 'ylim'))*[1, 1],'color', 'k');
                
            obj.h_txt = ...
                text(obj.h_ax, 3, 0.8*max(get(obj.h_ax, 'ylim')), txt_str, ...
                'verticalalignment', 'bottom', 'horizontalalignment', 'center', ...
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
        
        
        
        function val = ylim(obj, val)
            
            val = obj.get_set_limits(val, 'ylim');
            set(obj.h_txt, 'position', [0, val(2), 0]);
        end
    end
end