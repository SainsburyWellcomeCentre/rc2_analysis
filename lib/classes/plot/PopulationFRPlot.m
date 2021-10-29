classdef PopulationFRPlot < RC2Axis
    
    properties
        h_scale_bar
        scale_bar_s = 0.5
    end
    
    
    
    methods
        
        function obj = PopulationFRPlot(h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        
        function plot(obj, t, avg, sem)
            
            
            plot(obj.h_ax, t, avg, 'color', 'k', 'linewidth', 0.75);
            plot(obj.h_ax, t, avg+sem, 'color', [0.5, 0.5, 0.5], 'linewidth', 0.25);
            plot(obj.h_ax, t, avg-sem, 'color', [0.5, 0.5, 0.5], 'linewidth', 0.25);
            
            set(obj.h_ax, 'xlim', t([1, end]), ...
                          'xcolor', 'none', ...
                          'clipping', 'off', ...
                          'fontsize', 8);
            ylabel(obj.h_ax, '\DeltaFR (Hz)', 'fontsize', 8);
            yl = get(obj.h_ax, 'ylim');
            
            line(obj.h_ax, [0, 0], get(obj.h_ax, 'ylim'), 'color', 'k', 'linewidth', 0.5);
            
            obj.h_scale_bar = line(obj.h_ax, [t(end)-obj.scale_bar_s, t(end)], yl([1, 1]), ...
                               'color', 'k', ...
                               'linewidth', 0.5);
            
            text(obj.h_ax, t(end)-obj.scale_bar_s/2, yl(1)-(yl(2)-yl(1))*0.1, sprintf('%.1fs', obj.scale_bar_s), ...
                'color', 'k', ...
                'horizontalalignment', 'center', 'verticalalignment', 'top', ...
                'fontsize', 8);
            
        end
    end
end
        