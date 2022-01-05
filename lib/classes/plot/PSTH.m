classdef PSTH < RC2Axis
% PSTH Class for plotting PSTHs.
%
%  PSTH Properties:
%       bin_size  - size of the bin for the PSTH
%
%  PSTH Methods:
%       plot      - plot the PSTH


    properties (SetAccess = private)
        
        h_hist
        h_line
    end
    
    properties

        bin_size = 0.02;  % s
    end
    
    
    
    methods
        
        function obj = PSTH(h_ax)
        %%PSTH
        %
        %   PSTH(AXIS_HANDLE) sets up an object for handling
        %   PSTHs. AXIS_HANDLE is optional, if supplied it should be a
        %   handle to an axis object. Otherwise, an axis will be created.
        
            VariableDefault('h_ax', []);
            obj = obj@RC2Axis(h_ax);
        end
        
        
        function plot(obj, spike_times, t)
        %plot Plots the PSTH on the axis
        %
        %   plot(SPIKE_TIMES, LIMITS) where SPIKE_TIMES is a vector with a
        %   list of peri-event spike times, and LIMITS is a 1 x 2 vector
        %   with the limits over which to compute the histogram.
        %   Bin widths are set using the property `bin_size`
        
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
