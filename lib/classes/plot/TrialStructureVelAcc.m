classdef TrialStructureVelAcc < handle
% 
    properties (SetAccess = private)
        
        h_fig
        h_ax
        h_legend
        
        trial
        clusters
    end
    
    properties (Hidden = true)
        
        mask_scale = 20
        solenoid_scale = 5
        trig_shift = 3
        aw_scale = 22.5
        n_axes = 2
    end
    
    
    
    methods
        
        function obj = TrialStructureVelAcc()
            
            obj.h_fig = figure();
            set(obj.h_fig, 'position', [16, 80, 1440, 890], ...
                           'paperunits', 'inches', ...
                           'papersize', [15, 10], ...
                           'paperposition', [0, 0, 15, 10], ...
                           'renderer', 'painters');
        end
        
        
        function plot_velocity(obj, h_ax)
        %%plot_stage Plots the speed of the stage  
        % Velocity
            obj.plot_general_trace(h_ax, obj.trial.velocity, 'Velocity (cm/s)');
        end
        
        function plot_acceleration(obj, h_ax)
        %%plot_stage Plots the speed of the stage  
        % Velocity
            obj.plot_general_trace(h_ax, obj.trial.acceleration, 'Acceleration (m/s^2)');
        end
        

        
        function plot_general_trace(obj, h_ax, trace, label)
        %%plot_general_trace Shared function for plotting the traces
        
            ylabel(h_ax, label)
            obj.x_limits(h_ax);
            
            if isempty(trace)
                return
            end
            
            highlight = nan(size(trace));
            highlight(obj.trial.analysis_window()) = trace(obj.trial.analysis_window());
            
            % resample at 60hz for saving
            [trace, t] = obj.downsample_trace(trace);
            highlight = obj.downsample_trace(highlight);
            
            plot(h_ax, t, trace, 'color', 'k');
            plot(h_ax, t, highlight, 'color', obj.colours('aw'));
        end
        
        
        
        function [val, downsample_t] = downsample_trace(obj, trace)
        %%downsample_trace Downsamples trace to 60Hz
        
            n_samples = round((obj.trial.probe_t(end)-obj.trial.probe_t(1))*1000);
            downsample_t = linspace(obj.trial.probe_t(1), obj.trial.probe_t(end), n_samples);
            val = interp1(obj.trial.probe_t, trace, downsample_t);
        end
        
        
        
        function x_limits(obj, h_ax)
        %%x_limits Set x-axis limits and plotboxaspectratio
        
            set(h_ax, 'plotboxaspectratio', [10, 1, 1], 'box', 'off');
            set(h_ax, 'xlim', obj.trial.probe_t([1, end]));
        end
        

        
        function overlay_masks(obj, h_ax)
        %%overlay_masks Overlays the various masks on the plot
        %
        %   overlay_masks(AXIS_HANDLE) overlays "motion mask", "stationary
        %   mask", "analysis window", "solenoid trace", "teensy gain"
        
            plot(h_ax, obj.trial.probe_t, obj.mask_scale * obj.trial.motion_mask + obj.trig_shift, 'color', obj.colours('motion'));
            plot(h_ax, obj.trial.probe_t, obj.mask_scale * obj.trial.stationary_mask + obj.trig_shift, 'color', obj.colours('stationary'));
            plot(h_ax, obj.trial.probe_t, obj.solenoid_scale * obj.trial.solenoid + obj.trig_shift, 'color', obj.colours('solenoid'));
            plot(h_ax, obj.trial.probe_t, obj.aw_scale * obj.trial.analysis_window() + obj.trig_shift, 'color', obj.colours('aw'));
            if ~isempty(obj.trial.teensy_gain)
                plot(h_ax, obj.trial.probe_t, obj.mask_scale * obj.trial.teensy_gain + obj.trig_shift, 'color', obj.colours('gain'));
            end
        end
        
        
        
        function plot_on_axes(obj, n, trace_type)
        %%plot_on_axes Plots specific trace on specific axes
        %
        %   plot_on_axes(AXIS_NUMBER, TRACE_TYPE) plots a particular trace
        %   type on a particular axis. AXIS_NUMBER specifies which axis to
        %   plot the trace on (between 1 and `n_axes`), and TRACE_TYPE is
        %   one of 'treadmill', 'visual', 'stage', 'camera', and 'fr'.
        
            obj.h_ax{n} = subplot(obj.n_axes, 1, n); hold on;
            switch trace_type
                case 'velocity'
                    obj.plot_velocity(obj.h_ax{n});
                case 'acceleration'
                    obj.plot_acceleration(obj.h_ax{n});
            end
            obj.overlay_masks(obj.h_ax{n});
        end
        
        
        
        function add_text(obj, h_ax)
        %%add_text Adds time information about stationary and motion periods
        %
        %   add_text(AXIS_HANDLE) adds text information to the axis
        %   specified by AXIS_HANDLE.
        
            text_str = sprintf('Stationary time: %.2f s\nMotion time: %.2f s\nAnalysis window time: %.2f s', ...
                               obj.trial.stationary_time, ...
                               obj.trial.motion_time, ...
                               obj.trial.analysis_window_time);
            
            text(h_ax, obj.trial.probe_t(end), max(get(h_ax, 'ylim')), text_str, ...
                'horizontalalignment', 'right', ...
                'verticalalignment', 'top');
        end
        
        
        
        function plot(obj, trial, clusters)
        %%plot Main plotting function
        %
        %   plot(TRIAL, CLUSTERS) takes the information in the Trial object
        %   TRIAL, and the array of Cluster objects, CLUSTERS, and uses
        %   them to plot the trial structure.
        
            obj.trial = trial;
            obj.clusters = clusters;
            
            obj.plot_on_axes(1, 'velocity');
            obj.plot_on_axes(2, 'acceleration');
            
            obj.add_text(obj.h_ax{1});
            obj.h_legend = legend(obj.h_ax{2}, {'Speed', 'Analysis window', 'In motion', 'Stationary', 'Solenoid', 'Gain'});
            set(obj.h_legend, 'position', [0.835648152294289,0.64842667889739,0.099999998075267,0.11516853611121]);
            
            if obj.trial.is_replay
                title_str = sprintf('Trial %i, trial group %s, replay of %i', ...
                                    obj.trial.trial_id, ...
                                    obj.trial.trial_group_label, ...
                                    obj.trial.original_trial_id);
            else
                title_str = sprintf('Trial %i, trial group %s', ...
                                    obj.trial.trial_id, ...
                                    obj.trial.trial_group_label);
            end
            
            FigureTitle(obj.h_fig, title_str);
        end
    end
    
    
    
    methods (Static = true)
        
        function val = colours(type)
            
            cols = lines(7);
            
            switch type
                case 'solenoid'
                    val = cols(3, :);
                case 'stationary'
                    val = cols(1, :);
                case 'motion'
                    val = cols(2, :);
                case 'aw'
                    val = cols(5, :);
                case 'gain'
                    val = cols(6, :);
            end
        end
    end
end