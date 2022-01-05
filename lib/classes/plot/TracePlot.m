classdef TracePlot < RC2Axis
% TracePlot Class for plotting traces with SD/SEM around a mean, e.g.
% running speed or continuous firing rates.
%
%   TracePlot Properties:
%       h_line          - handle to origin line
%       h_mean_trace    - handle to mean trace line
%       h_traces        - handle to all traces
%       h_sd_traces     - handle to SD traces
%       h_sem_traces    - handle to SEM traces
%       traces          - all traces to plot (# samples x # traces)
%       x               - timebase on x-axiis (# samples x 1)
%
%   TracePlot Methods:
%       add_traces      - add all traces to the plot
%       remove_traces   - remove all traces from the plot
%       add_sd          - add S.D. lines to the plot
%       add_sem         - add S.E.M lines to the plot
%       remove_sem      - remove S.E.M. lines from the plot
%       mean_colour     - change the colour of the mean line
%
%   The constructor plots the mean of the traces provided

    properties
        
        h_line
        h_mean_trace
        h_traces
        h_sd_traces
        h_sem_traces
        
        traces
        x
    end
    
    methods
        
        function obj = TracePlot(traces, x, h_ax)
        %%TracePlot
        %
        %   TracePlot(TRACES, X, AXIS_HANDLE) creates a plot with the mean
        %   trace of TRACES plot, against X on the x-axis. AXIS_HANDLE is optional, if supplied it should be a
        %   handle to an axis object. Otherwise, an axis will be created.
        %   TRACES should be a # samples x # traces 2D matrix and X should
        %   be a vector of length # samples.
        
            VariableDefault('h_ax', []);
            obj = obj@RC2Axis(h_ax);
            
            
            % store the data
            obj.traces = traces;
            obj.x = x;
            
            obj.h_mean_trace = plot(obj.h_ax, obj.x, nanmean(obj.traces, 2), ...
                'color', 'r', 'linewidth', 2);
            
            set(obj.h_ax, 'plotboxaspectratio', [4, 1, 1])
            
            xlabel('Time (s)')
                      
            obj.h_line = line(obj.h_ax, [0, 0], get(obj.h_ax, 'ylim'), ...
                'color', 'k', 'linestyle', '--');
        end
        
        
        function add_traces(obj)
        %%add_traces Adds all traces in property `traces` to the plot.
        %
        %   add_traces() Adds the traces in light gray underneath the mean
        %   trace line.
        
            if ~isempty(obj.h_traces)
                return
            end
            
            if isempty(obj.x) || isempty(obj.traces)
                return
            end
            
            obj.h_traces = plot(obj.h_ax, obj.x, obj.traces, ...
                'color', [0.6, 0.6, 0.6], 'linewidth', 0.5);
            
            % bring mean trace to front
            if ~isempty(obj.h_sd_traces)
                uistack(obj.h_sd_traces, 'top');
            end
            uistack(obj.h_mean_trace, 'top');
        end
        
        
        
        function remove_traces(obj)
        %%remove_traces Removes all the traces from the plot, if they have
        %%been plotted using `add_traces()`
        %
        %   remove_traces() removes the traces if they exist.
        
            if isempty(obj.h_traces)
                return
            end
            
            delete(obj.h_traces);
            obj.h_traces = [];
        end
        
        
        
        function add_sd(obj)
        %%add_sd Adds two traces above and below the mean trace indicating
        %%the mean + or - the standard deviation
        %
        %   add_sd() plots the SD traces
        
            if ~isempty(obj.h_sd_traces)
                return
            end
            
            mean_trace = nanmean(obj.traces, 2);
            sd_trace = nanstd(obj.traces, [], 2);
            
            obj.h_sd_traces = plot(obj.h_ax, obj.x, ...
                [mean_trace + sd_trace, mean_trace - sd_trace], ...
                'color', [0.7, 0, 0], 'linewidth', 0.5);
        end
        
        
        
        function remove_sd(obj)
        %%remove_sd Removes the SD traces plotted with `add_sd()` if they
        %%exist
        %
        %   remove_sd() removes the SD traces if they exist.
        
            if isempty(obj.h_sd_traces)
                return
            end
            
            delete(obj.h_sd_traces);
            
            obj.h_sd_traces = [];
        end
        
        
        
        function add_sem(obj, col)
        %%add_sem Adds two traces above and below the mean trace indicating
        %%the mean + or - the standard error of the mean
        %
        %   add_sem() plots the SD traces
        
            if ~isempty(obj.h_sem_traces)
                return
            end
            
            mean_trace = nanmean(obj.traces, 2);
            sem_trace = nanstd(obj.traces, [], 2)./sqrt(sum(~isnan(obj.traces), 2));
            
            obj.h_sem_traces = plot(obj.h_ax, obj.x, ...
                [mean_trace + sem_trace, mean_trace - sem_trace], ...
                'color', col, 'linewidth', 0.5);
        end
        
        
        
        function remove_sem(obj)
        %%remove_sem Removes the SEM traces plotted with `add_sem()` if they
        %%exist
        %
        %   remove_sem() removes the SEM traces if they exist.
        
            if isempty(obj.h_sem_traces)
                return
            end
            
            delete(obj.h_sem_traces);
            
            obj.h_sem_traces = [];
        end
        
        
        
        function mean_colour(obj, col)
        %%mean_colour Changes the colour of the mean trace
        %
        %   mean_colour(COLOUR) where COLOUR is anything that set(axis,
        %   'color', COLOUR) would work.
        
            set(obj.h_mean_trace, 'color', col);
        end
    end
end