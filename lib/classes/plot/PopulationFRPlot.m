classdef PopulationFRPlot < RC2Axis
% PopulationFRPlot Class for plotting the population firing rate with SEM.
%
%  PopulationFRPlot Properties:
%       h_scale_bar     - handle to the scale bar
%       scale_bar_s     - size of the scale bar in seconds
%
%  PopulationFRPlot Methods:
%       plot            - do the plot
%
%   TODO:   1. do we really need this class?

    properties
        
        h_scale_bar
        scale_bar_s = 0.5
    end
    
    
    
    methods
        
        function obj = PopulationFRPlot(h_ax)
        %%PopulationFRPlot
        %
        %   PopulationFRPlot(AXIS_HANDLE). AXIS_HANDLE is optional, if supplied it should be a
        %   handle to an axis object. Otherwise, an axis will be created.
        
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        
        function plot(obj, t, avg, sem)
        %%plot Plots the firing rate trace
        %
        %   plot(TIME, MEAN_TRACE, SEM) plots MEAN_TRACE against TIME and
        %   also plots MEAN_TRACE + SEM against time.
            
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
        