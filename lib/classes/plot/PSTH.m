classdef PSTH < RC2Axis
    
    properties (SetAccess = private)
        
        h_hist
        h_line
    end
    
    properties

        bin_size = 0.02;  % s
    end
    
    
    
    methods
        
        function obj = PSTH(h_ax)
            VariableDefault('h_ax', []);
            obj = obj@RC2Axis(h_ax);
        end
        
        
        function plot(obj, spike_times, t)
            
            % make sure bins are centered around 0
            if t(end) > 0 && t(1) < 0
                n_bins1 = ceil(-t(1) / obj.bin_size);
                n_bins2 = ceil(t(end) / obj.bin_size);
                bin_edges = (-n_bins1:n_bins2)*obj.bin_size;
            else
                n_bins = ceil((t(end) - t(1)) / obj.bin_size);
                bin_edges = linspace(t(1), t(end), n_bins);
                
            end
            
            obj.h_hist = histogram(obj.h_ax, spike_times, bin_edges, 'facecolor', 'k', 'edgecolor', 'none');
            
            set(obj.h_ax, 'plotboxaspectratio', [4, 1, 1])
            
            xlabel(obj.h_ax, 'Time (s)')
                      
            obj.h_line = line(obj.h_ax, [0, 0], get(obj.h_ax, 'ylim'), 'color', 'k', 'linestyle', '--');
        end
    end
end
