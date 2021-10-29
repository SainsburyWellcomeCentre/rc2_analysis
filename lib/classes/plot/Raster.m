classdef Raster < RC2Axis
    
    properties
        
        h_line
        h_bars
    end
    
    methods
        
        function obj = Raster(spike_times, h_ax, marker_type)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
            
            n_trials = length(spike_times);
            
            % 
            for trial_i = 1 : n_trials
                
                n_spikes = length(spike_times{trial_i});
                
                if strcmp(marker_type, 'dot')
                    scatter(obj.h_ax, spike_times{trial_i}, ...
                        trial_i*ones(n_spikes, 1), ...
                        10, 'k', 'fill')
                elseif strcmp(marker_type, 'line')
                    for spike_i = 1 : n_spikes
                        line(obj.h_ax, ...
                            [spike_times{trial_i}, spike_times{trial_i}], ...
                            trial_i + [-.5, .5], ...
                            'color', 'k');
                    end
                end
            end
            
            if n_trials <= 1
                yt = 1;
            else
                yt = [1, n_trials];
            end
            
            set(obj.h_ax, 'plotboxaspectratio', [4, 1, 1], ...
                          'ydir', 'reverse', ...
                          'ylim', [0, n_trials+1], ...
                          'ytick', yt)
            
            xlabel(obj.h_ax, 'Time (s)')
            ylabel(obj.h_ax, 'Trial #')
                      
            obj.h_line = line(obj.h_ax, [0, 0], get(obj.h_ax, 'ylim'), 'color', 'k', 'linestyle', '--');
        end
        
        
        
        function add_bars(obj, t_start, t_end)
            
            if ~isempty(obj.h_bars)
                return
            end
            
            obj.h_bars = [];
            
            for i = 1 : length(t_start)
                obj.h_bars(i) = fill(obj.h_ax, [t_start(i), t_end, t_end, t_start(i)], ...
                                     [i-0.5, i-0.5, i+0.5, i+0.5], ...
                                     [0.6, 0.6, 0.6]);
                set(obj.h_bars(i), 'edgecolor', 'none');
            end    
        end
        
        
        function remove_bars(obj)
            
            if isempty(obj.h_bars)
                return
            end
            
            for i = 1 : length(obj.h_bars)
                delete(obj.h_bars(i));
            end
            
            obj.h_bars = [];
        end
        
        function ylabel(obj, str)
            txt = text('string', str);
            set(obj.h_ax, 'ylabel', txt);
        end
        
    end
end