classdef Raster < RC2Axis
% Raster Class for plotting rasters
%
%  Raster Properties:
%
%  Raster Methods:
%       add_bars      - adds grey bars to the end of the raster
%       remove_bars   - removes grey bars added by add_bars
%       ylabel        - adds a y-axis label
%
%   TODO:   1. Make consistent with other plot classes and separate
%              plotting from object creation

    properties
        
        h_line
        h_bars
    end
    
    
    
    methods
        
        function obj = Raster(spike_times, h_ax, marker_type, spike_colors, spike_alphas)
        %%Raster
        %
        %   Raster(SPIKE_TIMES, AXIS_HANDLE, MARKER_TYPE, SPIKE_COLORS, SPIKE_ALPHAS) plots data in
        %   SPIKE_TIMES.  AXIS_HANDLE is optional, if supplied it should be a
        %   handle to an axis object. Otherwise, an axis will be created.
        %   MARKER_TYPE is either 'dot' or 'line' and determines the style
        %   of the spikes.
        %   SPIKE_TIMES should be a cell array with each cell containing
        %   the times of the spikes for a separate row of the raster.
        %   SPIKE_COLORS (optional) is a cell array where each entry is a
        %   # spikes x 3 matrix of RGB colors for each spike.
        %   SPIKE_ALPHAS (optional) is a cell array where each entry is a
        %   # spikes x 1 vector of alpha values for each spike.
         
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
            
            n_trials = length(spike_times);
            
            % Handle optional color/alpha arguments
            use_per_spike_colors = false;
            if nargin >= 4 && ~isempty(spike_colors) && iscell(spike_colors)
                use_per_spike_colors = true;
            end
            
            % 
            for trial_i = 1 : n_trials
                
                n_spikes = length(spike_times{trial_i});
                
                if n_spikes == 0
                    continue;
                end
                
                if use_per_spike_colors && trial_i <= length(spike_colors)
                    % Per-spike colors and alphas - batch by unique color/alpha combinations
                    colors = spike_colors{trial_i};
                    if nargin >= 5 && ~isempty(spike_alphas) && iscell(spike_alphas) && trial_i <= length(spike_alphas)
                        alphas = spike_alphas{trial_i};
                    else
                        alphas = ones(n_spikes, 1);
                    end
                    
                    % Find unique color/alpha combinations
                    color_alpha_key = [colors, alphas];
                    [unique_keys, ~, group_idx] = unique(color_alpha_key, 'rows', 'stable');
                    
                    % Plot each unique color/alpha group together
                    for group = 1:size(unique_keys, 1)
                        mask = group_idx == group;
                        group_spike_times = spike_times{trial_i}(mask);
                        group_color = unique_keys(group, 1:3);
                        group_alpha = unique_keys(group, 4);
                        
                        n_group_spikes = length(group_spike_times);
                        
                        if strcmp(marker_type, 'dot')
                            h = scatter(obj.h_ax, group_spike_times, ...
                                trial_i*ones(n_group_spikes, 1), ...
                                10, group_color, 'fill');
                            h.MarkerFaceAlpha = group_alpha;
                            h.MarkerEdgeAlpha = group_alpha;
                        elseif strcmp(marker_type, 'line')
                            for spike_i = 1:n_group_spikes
                                line(obj.h_ax, ...
                                    [group_spike_times(spike_i), group_spike_times(spike_i)], ...
                                    trial_i + [-.5, .5], ...
                                    'color', [group_color, group_alpha]);
                            end
                        end
                    end
                else
                    % Default: all black
                    if strcmp(marker_type, 'dot')
                        scatter(obj.h_ax, spike_times{trial_i}, ...
                            trial_i*ones(n_spikes, 1), ...
                            10, 'k', 'fill');
                    elseif strcmp(marker_type, 'line')
                        for spike_i = 1 : n_spikes
                            line(obj.h_ax, ...
                                [spike_times{trial_i}(spike_i), spike_times{trial_i}(spike_i)], ...
                                trial_i + [-.5, .5], ...
                                'color', 'k');
                        end
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
        %%add_bars Adds grey bars at the end of each raster line
        %
        %   add_bars(START, END) draw horizontal grey bars at the end of
        %   each raster line to indicate no data in this period. START is a
        %   # raster lines x 1 vector with the start of the grey bar on
        %   each line, and END is a single value indicating where to stop
        %   the bars.
        
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
        %%remove_bars Removes the grey bars added with add_bars, if they exist.
        
            if isempty(obj.h_bars)
                return
            end
            
            for i = 1 : length(obj.h_bars)
                delete(obj.h_bars(i));
            end
            
            obj.h_bars = [];
        end
        
        
        
        function ylabel(obj, str)
        %%ylabel Adds a ylabel text to the axis (what is the usefulness of this?)
        
            txt = text('string', str);
            set(obj.h_ax, 'ylabel', txt);
        end
        
    end
end