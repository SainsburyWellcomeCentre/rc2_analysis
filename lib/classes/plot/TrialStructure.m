classdef TrialStructure < handle
% TrialStructure Class for plotting the structure of RVT trials.
%
%  TrialStructure Properties:
%       mask_scale       - where to plot the max mask value on y-axis
%       solenoid_scale   - where to plot the solenoid on y-axis
%       trig_shift       - how much to shift the masks above 0 on y-axis
%       aw_scale         - where to plot the max analysis window value on y-axis
%       n_axes           - number of axes (5)
%
%  TrialStructure Methods:
%       plot_treadmill      - plots the treadmill velocity
%       plot_visual         - plots the visual motion velocity
%       plot_stage          - plots the stage velocity
%       plot_camera         - plots the camera motion energy
%       plot_general_trace  - shared function for plotting traces
%       downsample_trace    - downsamples the traces
%       x_limits            - set x-axis limits and plotboxaspectratio
%       plot_fr             - plots the firing rate trace
%       overlay_masks       - plots all the masks
%       plot_on_axes        - plots specific trace on specific axes
%       add_text            - adds time information about stationary and motion periods
%       plot                - main plotting function
%
%   The TrialStructure class takes information about a trial and plots the
%   breakdown of the trial. This includes separate axes for treadmill,
%   visual and linear stage speeds, as well as the camera motion energy and
%   firing rate of a population of clusters.
%
%   The main method is `plot`, which takes a Trial object and an array of
%   Cluster objects. The other methods are largely for internal use.
%   
%
%   See also: Trial, Cluster

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
        n_axes = 7
    end
    
    
    
    methods
        
        function obj = TrialStructure()
            
            obj.h_fig = figure();
            set(obj.h_fig, 'position', [16, 80, 1440, 890], ...
                           'paperunits', 'inches', ...
                           'papersize', [15, 10], ...
                           'paperposition', [0, 0, 15, 10], ...
                           'renderer', 'painters');
        end
        
        
        
        function plot_treadmill(obj, h_ax)
        %%plot_treadmill Plots the treadmill speed
        
            obj.plot_general_trace(h_ax, obj.trial.treadmill_speed, 'R (cm/s)');
        end
        
        
        
        function plot_visual(obj, h_ax)
        %%plot_visual Plots the speed of the visual motion
        
            obj.plot_general_trace(h_ax, obj.trial.visual_speed, 'V (cm/s)');
        end
        
        
        
        function plot_stage(obj, h_ax)
        %%plot_stage Plots the speed of the stage    
            obj.plot_general_trace(h_ax, obj.trial.stage, 'T (cm/s)');
        end
        

        function plot_pupil_diameter(obj, h_ax)
        %%TODO
            obj.plot_general_trace(h_ax, obj.mask_scale * obj.trial.pupil_diameter/max(obj.trial.pupil_diameter), 'Pupil diameter (pixels)');            

        end

        
        function plot_camera0(obj, h_ax)
        %%plot_camera Plots the motion energy from the camera
            obj.plot_general_trace(h_ax, obj.mask_scale * obj.trial.camera0/max(obj.trial.camera0), 'Camera 0 (a.u.)');            

        end
        
        
        function plot_camera1(obj, h_ax)
        %%plot_camera Plots the motion energy from the camera
            obj.plot_general_trace(h_ax, obj.mask_scale * obj.trial.camera1/max(obj.trial.camera1), 'Camera 1 (a.u.)');
            
%             lims = prctile(obj.trial.camera1, [0, 30]);
%             m = median(obj.trial.camera1(obj.trial.camera1 > lims(1) & obj.trial.camera1 < lims(2)));
%             new_cam = obj.mask_scale * (obj.trial.camera1 - m)/max(obj.trial.camera1 - m);
%             plot(h_ax, obj.trial.camera_t, new_cam, 'color', 'k');
%             ylabel(h_ax, 'Camera (a.u.)')
%             obj.x_limits(h_ax);
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
        
            n_samples = round((obj.trial.probe_t(end)-obj.trial.probe_t(1))*60);
            downsample_t = linspace(obj.trial.probe_t(1), obj.trial.probe_t(end), n_samples);
            val = interp1(obj.trial.probe_t, trace, downsample_t);
        end
        
        
        
        function x_limits(obj, h_ax)
        %%x_limits Set x-axis limits and plotboxaspectratio
        
            set(h_ax, 'plotboxaspectratio', [10, 1, 1], 'box', 'off');
            set(h_ax, 'xlim', obj.trial.probe_t([1, end]));
        end
        
        
        function plot_fr(obj, h_ax)
        %%plot_fr Plots firing rate trace for population
        
            cluster_fr = cell(1, length(obj.clusters));
            for ii = 1 : length(obj.clusters)
                cluster_fr{ii} = obj.clusters(ii).fr.get_convolution(obj.trial.probe_t);
            end
            
            population_fr = mean(cat(2, cluster_fr{:}), 2);
            obj.plot_general_trace(h_ax, population_fr, 'FR (Hz)');
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
                case 'treadmill'
                    obj.plot_treadmill(obj.h_ax{n});
                case 'visual'
                    obj.plot_visual(obj.h_ax{n});
                case 'stage'
                    obj.plot_stage(obj.h_ax{n});
                case 'camera1'
                    obj.plot_camera1(obj.h_ax{n});
                case 'camera0'
                    obj.plot_camera0(obj.h_ax{n});
                case 'pupil_diameter'
                    obj.plot_pupil_diameter(obj.h_ax{n});
                case 'fr'
                    obj.plot_fr(obj.h_ax{n});
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
            
            obj.plot_on_axes(1, 'treadmill');
            obj.plot_on_axes(2, 'visual');
            obj.plot_on_axes(3, 'stage');
            obj.plot_on_axes(4, 'camera1');
            obj.plot_on_axes(5, 'camera0');
            obj.plot_on_axes(6, 'pupil_diameter');
            obj.plot_on_axes(7, 'fr');
            
            obj.add_text(obj.h_ax{1});
            xlabel(obj.h_ax{5}, 'Time (s)');
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