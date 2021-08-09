classdef TrialStructure < handle
    
    properties
        
        h_ax
        
    end
    
    
    methods
        
        function obj = TrialStructure(trial, clusters)
            
            VariableDefault('clusters', []);
            
            if isempty(clusters)
                n_plots = 4;
            else
                n_plots = 5;
            end
            
            cols = lines(7);
            
            solenoid_col = cols(3, :);
            stat_col = cols(1, :);
            mot_col = cols(2, :);
            aw_col = cols(5, :);
            gain_col = cols(6, :);
            
            h = figure;
            set(h, 'position', [16, 80, 1440, 890], 'papersize', [15, 10], ...
                'renderer', 'painters');
            
            if ~isempty(trial.filtered_teensy)
                t = (0:length(trial.filtered_teensy)-1) * (1/trial.fs);
                f = trial.filtered_teensy;
            else
                t = (0:length(trial.filtered_teensy_2)-1) * (1/trial.fs);
                f = trial.filtered_teensy_2;
            end
            idx = find(trial.analysis_window());
            trig_shift = 3;
            
            obj.h_ax = subplot(n_plots, 1, 1);
            hold on;
            
            
            plot(t, f, 'color', 'k');
            g = nan(size(f));
            g(idx) = f(idx);
            plot(t, g, 'color', aw_col);
            plot(t, 20*trial.motion_mask+trig_shift, 'color', mot_col);
            plot(t, 20*trial.stationary_mask+trig_shift, 'color', stat_col);
            plot(t, 5*trial.solenoid+trig_shift, 'color', solenoid_col);
            plot(t, 22.5*trial.analysis_window()+trig_shift, 'color', aw_col);
            if ~isempty(trial.teensy_gain)
                plot(t, 20*trial.teensy_gain+trig_shift, 'color', gain_col);
            end
            
            
            text(t(end), max(get(gca, 'ylim')), ...
                sprintf('Stationary time: %.2f s\nMotion time: %.2f s\nAnalysis window time: %.2f s', ...
                trial.stationary_time, trial.motion_time, trial.analysis_window_time), ...
                'horizontalalignment', 'right', 'verticalalignment', 'top');
            ylabel('M (cm/s)')
            set(gca, 'plotboxaspectratio', [10, 1, 1]);
            box off;
            
            
            subplot(n_plots, 1, 2);
            hold on;
            
            if isnumeric(trial.config.enable_vis_stim)
                enable_vis_stim = trial.config.enable_vis_stim;
            else
                enable_vis_stim = str2double(trial.config.enable_vis_stim);
            end
            
            if ~isempty(trial.multiplexer_output) && enable_vis_stim
                f = trial.multiplexer_output;
                plot(t, f, 'color', 'k');
                g = nan(size(f));
                g(idx) = f(idx);
                plot(t, g, 'color', aw_col);
            else
                f = zeros(length(t), 1);
                plot(t, f, 'color', 'k');
                g = nan(size(f));
                g(idx) = f(idx);
                plot(t, g, 'color', aw_col);
            end
            
            plot(t, 20*trial.motion_mask+trig_shift, 'color', mot_col);
            plot(t, 20*trial.stationary_mask+trig_shift, 'color', stat_col);
            plot(t, 5*trial.solenoid+trig_shift, 'color', solenoid_col);
            plot(t, 22.5*trial.analysis_window()+trig_shift, 'color', aw_col);
            if ~isempty(trial.teensy_gain)
                plot(t, 20*trial.teensy_gain+trig_shift, 'color', gain_col);
            end
            ylabel('V (cm/s)')
            set(gca, 'plotboxaspectratio', [10, 1, 1]);
            box off;
            legend({'Speed', 'Analysis window', 'In motion', 'Stationary', 'Solenoid'})
            
            subplot(n_plots, 1, 3);
            hold on;
            f = trial.stage;
            plot(t, f, 'color', 'k');
            g = nan(size(f));
            g(idx) = f(idx);
            plot(t, g, 'color', aw_col);
            plot(t, 20*trial.motion_mask+trig_shift, 'color', mot_col);
            plot(t, 20*trial.stationary_mask+trig_shift, 'color', stat_col);
            plot(t, 5*trial.solenoid+trig_shift, 'color', solenoid_col);
            plot(t, 22.5*trial.analysis_window()+trig_shift, 'color', aw_col);
            if ~isempty(trial.teensy_gain)
                plot(t, 20*trial.teensy_gain+trig_shift, 'color', gain_col);
            end
            ylabel('T (cm/s)')
            set(gca, 'plotboxaspectratio', [10, 1, 1]);
            box off;
            
            
            if ~isempty(trial.camera1)
                subplot(n_plots, 1, 4);
                hold on;
                y = 20*(trial.camera1 - min(trial.camera1))/(max(trial.camera1) - min(trial.camera1));
                plot(t, y, 'color', 'k');
                g = nan(size(y));
                g(idx) = y(idx);
                plot(t, g, 'color', aw_col);
                plot(t, 20*trial.motion_mask+trig_shift, 'color', mot_col);
                plot(t, 20*trial.stationary_mask+trig_shift, 'color', stat_col);
                plot(t, 5*trial.solenoid+trig_shift, 'color', solenoid_col);
                plot(t, 22.5*trial.analysis_window()+trig_shift, 'color', aw_col);
                if ~isempty(trial.teensy_gain)
                    plot(t, 20*trial.teensy_gain+trig_shift, 'color', gain_col);
                end
                ylabel('Camera (a.u.)')
                xlabel('Time (s)');
                set(gca, 'plotboxaspectratio', [10, 1, 1]);
                box off;
            end
            
            
            
            if ~isempty(clusters)
                
                subplot(n_plots, 1, 5);
                hold on;
                
                r = {};
                for i = 1 : length(clusters)
                    fr = FiringRate(clusters(i).spike_times);
                    r{i} = fr.get_convolution(trial.probe_t);
                end
                
                y = mean(cat(2, r{:}), 2);
                plot(t, y, 'color', 'k');
                g = nan(size(y));
                g(idx) = y(idx);
                plot(t, g, 'color', aw_col);
                plot(t, 20*trial.motion_mask+trig_shift, 'color', mot_col);
                plot(t, 20*trial.stationary_mask+trig_shift, 'color', stat_col);
                plot(t, 5*trial.solenoid+trig_shift, 'color', solenoid_col);
                plot(t, 22.5*trial.analysis_window()+trig_shift, 'color', aw_col);
                if ~isempty(trial.teensy_gain)
                    plot(t, 20*trial.teensy_gain+trig_shift, 'color', gain_col);
                end
                ylabel('Firing rate (Hz)')
                xlabel('Time (s)');
                set(gca, 'plotboxaspectratio', [10, 1, 1]);
                box off;
            end
        end
    end
end