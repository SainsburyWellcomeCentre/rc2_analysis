classdef TracePlot < RC2Axis
    
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
            
            if isempty(obj.h_traces)
                return
            end
            
            delete(obj.h_traces);
            obj.h_traces = [];
        end
        
        
        
        function add_sd(obj)
            
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
            
            if isempty(obj.h_sd_traces)
                return
            end
            
            delete(obj.h_sd_traces);
            
            obj.h_sd_traces = [];
        end
        
        
        
        function add_sem(obj, col)
            
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
            
            if isempty(obj.h_sem_traces)
                return
            end
            
            delete(obj.h_sem_traces);
            
            obj.h_sem_traces = [];
        end
        
        
        
        function mean_colour(obj, col)
            set(obj.h_mean_trace, 'color', col);
        end
    end
end