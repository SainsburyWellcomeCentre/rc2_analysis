% plot_all_trials_matched_across_replays_analysis_window.m
%
%   For each visual_flow recording, plot velocity of all trials.
%       Match trials from MVT and MV to their replays in VT and V.
%
%       Color the analysis window, determined from the MVT and MV trials
%       and matched onto the VT and V



clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'visual_flow';  % 'darkness', 'visual_flow', 'head_tilt', 'mismatch_nov20'
combination = 'protocols';

config = RC2AnalysisConfig();
figs = RC2Figures(config);

% where to save figure
figs.save_on = true;
figs.set_figure_subdir('visual_flow', 'all_trials_analysis_window');

% get details of the experiment
probe_fnames = experiment_details(experiment, combination);

trial_types = {'Coupled', 'EncoderOnly'};
fname_suffix = {'MVT', 'MV'};

for probe_i = 1 : length(probe_fnames)
    
    % all data for this probe recording
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    
    % get visual flow experiment object
    vf = VisualFlowExperiment(data, config);
    
    base_t = {};
    M_trace = {};
    V_trace = {};
    T_trace = {};
    trial_id = {};
    analysis_window = {};
    
    base_t_rep = {};
    M_trace_rep = {};
    V_trace_rep = {};
    T_trace_rep = {};
    trial_id_rep = {};
    bouts_rep = {};
    protocol_type_rep = {};
    analysis_window_rep = {};
    
    for type_i = 1 : length(trial_types)
        
        these_trials = vf.trials_of_type(trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            
            base_t{type_i}{trial_i} = these_trials(trial_i).probe_t - ...
                these_trials(trial_i).probe_t(1);
            M_trace{type_i}{trial_i} = these_trials(trial_i).filtered_teensy;
            V_trace{type_i}{trial_i} = these_trials(trial_i).multiplexer_output;
            T_trace{type_i}{trial_i} = these_trials(trial_i).stage;
            trial_id{type_i}{trial_i} = these_trials(trial_i).id;
            analysis_window{type_i}{trial_i} = these_trials(trial_i).analysis_window();
            
            % get repeats of trials
            repeats = vf.get_replay_of(these_trials(trial_i));
            
            % match
            for rep_i = 1 : length(repeats)
                
                %
                offset = vf.get_offset(these_trials(trial_i), repeats(rep_i));
                
                n_remaining = length(repeats(rep_i).probe_t) - offset;
                n_to_plot = min(n_remaining, length(these_trials(trial_i).probe_t));
                
                base_t_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).probe_t(offset+(0:n_to_plot-1)) - ...
                    repeats(rep_i).probe_t(offset);
                M_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).filtered_teensy(offset+(0:n_to_plot-1));
                V_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).multiplexer_output(offset+(0:n_to_plot-1));
                T_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).stage(offset+(0:n_to_plot-1));
                trial_id_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).id;
                protocol_type_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).protocol;
                
                analysis_window_rep{type_i}{trial_i}{rep_i} = analysis_window{type_i}{trial_i}(1:n_to_plot);
                
            end
            
        end
        
    end
    
    
    all_traces = {M_trace, V_trace, T_trace};
    all_traces_rep = {M_trace_rep, V_trace_rep, T_trace_rep};
    axis_labels = {'M', 'V', 'T'};
    
    end_pad = 10e3;
    trace_offset = 10;
    aw_col = [0.8, 0.6, 0.4];
    
    %% plot each trial
    for type_i = 1 : length(base_t)
        
        for trial_i = 1 : length(base_t{type_i})
            
            h_fig = figs.a4figure();
            
            h_ax = [];
            yL = 0;
            
            h_ax(end+1) = subplot(3, 2, 1);
            hold on;
            
            for trace_i = 1 : length(all_traces)
                
                t = base_t{type_i}{trial_i};
                y = all_traces{trace_i}{type_i}{trial_i} - (trace_i-1)*trace_offset;
                trace_plot(h_ax(end), t, y, [], []);
                
                idx = analysis_window{type_i}{trial_i};
                t = base_t{type_i}{trial_i}(idx);
                y = all_traces{trace_i}{type_i}{trial_i}(idx) - (trace_i-1)*trace_offset;
                
                plot(t, y, 'color', aw_col);
                
                set(h_ax(end), 'ytick', []);
                
                title(sprintf('Trial %i, %s', ...
                    trial_id{type_i}{trial_i}, ...
                    fname_suffix{type_i}));
                
            end
            
            line(max(get(h_ax(end), 'xlim'))*[1, 1], [5, 25], 'color', 'k');
            text(max(get(h_ax(end), 'xlim')), 15, '20cm/s', 'rotation', 90, 'horizontalalignment', 'center', 'verticalalignment', 'top');
            
            yL = max(yL, max(get(gca, 'ylim')));
            
            h_ax(end+1) = subplot(3, 2, 2);
            hold on;
            
            for trace_i = 1 : length(all_traces)
                
                idx = find(analysis_window{type_i}{trial_i}, 1, 'last');
                
                t = base_t{type_i}{trial_i}(idx + (-end_pad:end_pad));
                y = all_traces{trace_i}{type_i}{trial_i}(idx + (-end_pad:end_pad)) - (trace_i-1)*trace_offset;
                trace_plot(h_ax(end), t, y, [], []);
                
                t = base_t{type_i}{trial_i}(idx + (-end_pad:0));
                y = all_traces{trace_i}{type_i}{trial_i}(idx + (-end_pad:0)) - (trace_i-1)*trace_offset;
                
                plot(t, y, 'color', aw_col);
                
                title(sprintf('Trial %i, %s', ...
                    trial_id{type_i}{trial_i}, ...
                    fname_suffix{type_i}));
                
            end
            
            yL = max(yL, max(get(gca, 'ylim')));
            
            for rep_i = 1 : length(base_t_rep{type_i}{trial_i})
                
                if strcmp(protocol_type_rep{type_i}{trial_i}{rep_i}, 'StageOnly')
                    
                    sp_n = 3;
                    title_str = 'VT';
                    
                else
                    
                    sp_n = 5;
                    title_str = 'V';
                    
                end
                
                
                h_ax(end+1) = subplot(3, 2, sp_n);
                hold on;
                
                for trace_i = 1 : length(all_traces_rep)
                    
                    t = base_t_rep{type_i}{trial_i}{rep_i};
                    y = all_traces_rep{trace_i}{type_i}{trial_i}{rep_i} - (trace_i-1)*trace_offset;
                    trace_plot(h_ax(end), t, y, [], []);
                    
                    idx = analysis_window_rep{type_i}{trial_i}{rep_i};
                    t = base_t_rep{type_i}{trial_i}{rep_i}(idx);
                    y = all_traces_rep{trace_i}{type_i}{trial_i}{rep_i}(idx) - (trace_i-1)*trace_offset;
                    
                    plot(t, y, 'color', aw_col);
                    
                    set(h_ax(end), 'ytick', []);
                    
                    title(sprintf('Trial %i, %s (%s)', ...
                        trial_id{type_i}{trial_i}, ...
                        title_str, ...
                        fname_suffix{type_i}));
                    
                end
                
                yL = max(yL, max(get(gca, 'ylim')));
                
                h_ax(end+1) = subplot(3, 2, sp_n + 1);
                hold on;
                
                for trace_i = 1 : length(all_traces_rep)
                    
                    idx = find(analysis_window_rep{type_i}{trial_i}{rep_i}, 1, 'last');
                    
                    t = base_t_rep{type_i}{trial_i}{rep_i}(idx + (-end_pad:end_pad));
                    y = all_traces_rep{trace_i}{type_i}{trial_i}{rep_i}(idx + (-end_pad:end_pad)) - (trace_i-1)*trace_offset;
                    trace_plot(h_ax(end), t, y, [], []);
                    
                    t = base_t_rep{type_i}{trial_i}{rep_i}(idx + (-end_pad:0));
                    y = all_traces_rep{trace_i}{type_i}{trial_i}{rep_i}(idx + (-end_pad:0)) - (trace_i-1)*trace_offset;
                    
                    plot(t, y, 'color', aw_col);
                    
                    title(sprintf('Trial %i, %s (%s)', ...
                        trial_id_rep{type_i}{trial_i}{rep_i}, ...
                        title_str, ...
                        fname_suffix{type_i}));
                    
                end
                
                yL = max(yL, max(get(gca, 'ylim')));
                
            end
                
            for ax_i = 1 : length(h_ax)
                
                set(h_ax(ax_i), 'ylim', [-30, yL]);
                
            end
            
            
            FigureTitle(h_fig, sprintf('%s, %s, Trial %i', ...
                probe_fnames{probe_i}, ...
                trial_types{type_i}, ...
                trial_id{type_i}{trial_i}));
            
            figs.save_fig_to_join();
            
        end
        
        figs.join_figs( ...
            sprintf('%s_%s_all_trials_analysis_window.pdf', ...
            probe_fnames{probe_i}, ...
            fname_suffix{type_i}));
        figs.clear_figs();
        
    end
    
end






