experiment = 'darkness';

config = RC2AnalysisConfig();

figs = RC2Figures(config);
figs.save_on = true;
figs.set_figure_subdir(experiment, 'all_trials_analysis_window');

% get details of the experiment
probe_fnames = experiment_details(experiment, 'protocols');

main_trial_types    = {'Coupled', 'EncoderOnly', 'StageOnly'};

for probe_i = 1:2% : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    dark            = DarknessExperiment(data, config);
    
    base_t_pair = {};
    M_trace_pair = {};
    V_trace_pair = {};
    T_trace_pair = {};
    trial_id_pair = {};
    protocol_type_pair = {};
    analysis_window_pair = {};
    
    base_t = {};
    M_trace = {};
    V_trace = {};
    T_trace = {};
    trial_id = {};
    protocol_type = {};
    analysis_window = {};
    
    for type_i = 1 : length(main_trial_types)
        
        these_trials = dark.trials_of_type(main_trial_types{type_i});
        
        % seriously???
        base_t_pair{type_i} = {{}, {}};
        M_trace_pair{type_i} = {{}, {}};
        V_trace_pair{type_i} = {{}, {}};
        T_trace_pair{type_i} = {{}, {}};
        trial_id_pair{type_i} = {{}, {}};
        protocol_type_pair{type_i} = {{}, {}};
        analysis_window_pair{type_i} = {{}, {}};
        
        base_t{type_i} = {{}};
        M_trace{type_i} = {{}};
        V_trace{type_i} = {{}};
        T_trace{type_i} = {{}};
        trial_id{type_i} = {{}};
        protocol_type{type_i} = {{}};
        analysis_window{type_i} = {{}};
        
        for trial_i = 1 : length(these_trials)
            
            this_trial = these_trials(trial_i);
            
            replay = dark.get_replay_of(this_trial);
            replayed_trial = dark.get_replayed_trial(this_trial);
            
            if ~isempty(replay)
                store_pair = true;
            elseif ~isempty(replayed_trial)
                continue
            else
                store_pair = false;
            end
            
            if store_pair
                
                offset = dark.get_offset(this_trial, replay);
                replay = AlignedTrial(replay, this_trial, offset);
                
                base_t_pair{type_i}{1}{end+1} = this_trial.probe_t - this_trial.probe_t(1);
                M_trace_pair{type_i}{1}{end+1} = this_trial.filtered_teensy;
                if strcmp(experiment, 'darkness')
                    V_trace_pair{type_i}{1}{end+1} = zeros(length(this_trial.filtered_teensy), 1);
                else
                    V_trace_pair{type_i}{1}{end+1} = this_trial.multiplexer_output;
                end
                T_trace_pair{type_i}{1}{end+1} = this_trial.stage;
                trial_id_pair{type_i}{1}{end+1} = this_trial.id;
                protocol_type_pair{type_i}{1}{end+1} = this_trial.protocol;
                analysis_window_pair{type_i}{1}{end+1} = this_trial.analysis_window();
                
                pt = replay.probe_t;
                base_t_pair{type_i}{2}{end+1} = pt - pt(1);
                M_trace_pair{type_i}{2}{end+1} = replay.filtered_teensy;
                if strcmp(experiment, 'darkness')
                    V_trace_pair{type_i}{2}{end+1} = zeros(length(replay.filtered_teensy), 1);
                else
                    V_trace_pair{type_i}{2}{end+1} = replay.multiplexer_output;
                end
                T_trace_pair{type_i}{2}{end+1} = replay.stage;
                trial_id_pair{type_i}{2}{end+1} = replay.id;
                protocol_type_pair{type_i}{2}{end+1} = replay.protocol;
                analysis_window_pair{type_i}{2}{end+1} = replay.analysis_window();
                
            else
                
                base_t{type_i}{1}{end+1} = this_trial.probe_t - this_trial.probe_t(1);
                M_trace{type_i}{1}{end+1} = this_trial.filtered_teensy;
                if strcmp(experiment, 'darkness')
                    V_trace{type_i}{1}{end+1} = zeros(length(this_trial.filtered_teensy), 1);
                else
                    V_trace{type_i}{1}{end+1} = this_trial.multiplexer_output;
                end
                T_trace{type_i}{1}{end+1} = this_trial.stage;
                trial_id{type_i}{1}{end+1} = this_trial.id;
                protocol_type{type_i}{1}{end+1} = this_trial.protocol;
                analysis_window{type_i}{1}{end+1} = this_trial.analysis_window();
                
            end
            
        end
    end
    
    
    all_traces = {M_trace, V_trace, T_trace};
    all_traces_pair = {M_trace_pair, V_trace_pair, T_trace_pair};
    axis_labels = {'M', 'V', 'T'};
    
    end_pad = 10e3;
    trace_offset = 10;
    aw_col = [0.8, 0.6, 0.4];
    
    %% plot each trial
    for type_i = 1 : length(base_t_pair)
        
        for trial_i = 1 : length(base_t_pair{type_i}{1})
            
            h_fig = figs.a4figure();
            
            h_ax = [];
            yL = 0;
            
            for pair_i = 1 : 2
                
                h_ax(end+1) = subplot(3, 2, (pair_i-1)*2+1);
                hold on;
                
                for trace_i = 1 : 3
                    
                    t = base_t_pair{type_i}{pair_i}{trial_i};
                    y = all_traces_pair{trace_i}{type_i}{pair_i}{trial_i} - (trace_i-1)*trace_offset;
                    trace_plot(h_ax(end), t, y, [], []);
                    
                    idx = analysis_window_pair{type_i}{pair_i}{trial_i};
                    t = base_t_pair{type_i}{pair_i}{trial_i}(idx);
                    y = all_traces_pair{trace_i}{type_i}{pair_i}{trial_i}(idx) - (trace_i-1)*trace_offset;
                    
                    plot(t, y, 'color', aw_col);
                    
                    set(h_ax(end), 'ytick', []);
                    
                    title(sprintf('Trial %i, %s', ...
                        trial_id_pair{type_i}{pair_i}{trial_i}, ...
                        protocol_type_pair{type_i}{pair_i}{trial_i}));
                    
                end
                
                line(max(get(h_ax(end), 'xlim'))*[1, 1], [5, 25], 'color', 'k');
                text(max(get(h_ax(end), 'xlim')), 15, '20cm/s', 'rotation', 90, 'horizontalalignment', 'center', 'verticalalignment', 'top');
                
                yL = max(yL, max(get(gca, 'ylim')));
                
                
                h_ax(end+1) = subplot(3, 2, 2*pair_i);
                hold on;
                
                for trace_i = 1 : 3
                    
                    idx = find(analysis_window_pair{type_i}{pair_i}{trial_i}, 1, 'last');
                    
                    t = base_t_pair{type_i}{pair_i}{trial_i}(idx + (-end_pad:end_pad));
                    y = all_traces_pair{trace_i}{type_i}{pair_i}{trial_i}(idx + (-end_pad:end_pad)) - (trace_i-1)*trace_offset;
                    trace_plot(h_ax(end), t, y, [], []);
                    
                    t = base_t_pair{type_i}{pair_i}{trial_i}(idx + (-end_pad:0));
                    y = all_traces_pair{trace_i}{type_i}{pair_i}{trial_i}(idx + (-end_pad:0)) - (trace_i-1)*trace_offset;
                    
                    plot(t, y, 'color', aw_col);
                    
                    title(sprintf('Trial %i, %s', ...
                        trial_id_pair{type_i}{pair_i}{trial_i}, ...
                        protocol_type_pair{type_i}{pair_i}{trial_i}));
                    
                end
                
                yL = max(yL, max(get(gca, 'ylim')));
                
            end
            
                
            for ax_i = 1 : length(h_ax)
                
                set(h_ax(ax_i), 'ylim', [-30, yL]);
                
            end
            
            
            FigureTitle(h_fig, sprintf('%s, %i', ...
                probe_fnames{probe_i}, ...
                trial_id_pair{type_i}{1}{trial_i}));
            
            figs.save_fig_to_join();
            
        end
        
    end
    
    
    %% plot each trial
    for type_i = 1 : length(base_t)
        
        for trial_i = 1 : length(base_t{type_i}{1})
            
            if isempty(base_t{type_i}{1})
                continue
            end
            
            h_fig = figs.a4figure();
            
            h_ax = [];
            yL = 0;
            
            pair_i = 1;
            
            h_ax(end+1) = subplot(3, 2, 1);
            hold on;
            
            for trace_i = 1 : 3
                
                t = base_t{type_i}{pair_i}{trial_i};
                y = all_traces{trace_i}{type_i}{pair_i}{trial_i} - (trace_i-1)*trace_offset;
                trace_plot(h_ax(end), t, y, [], []);
                
                idx = analysis_window{type_i}{pair_i}{trial_i};
                t = base_t{type_i}{pair_i}{trial_i}(idx);
                y = all_traces{trace_i}{type_i}{pair_i}{trial_i}(idx) - (trace_i-1)*trace_offset;
                
                plot(t, y, 'color', aw_col);
                
                set(h_ax(end), 'ytick', []);
                
                title(sprintf('Trial %i, %s', ...
                    trial_id{type_i}{pair_i}{trial_i}, ...
                    protocol_type{type_i}{pair_i}{trial_i}));
                
            end
            
            line(max(get(h_ax(end), 'xlim'))*[1, 1], [5, 25], 'color', 'k');
            text(max(get(h_ax(end), 'xlim')), 15, '20cm/s', 'rotation', 90, 'horizontalalignment', 'center', 'verticalalignment', 'top');
            
            yL = max(yL, max(get(gca, 'ylim')));
            
            
            h_ax(end+1) = subplot(3, 2, 2);
            hold on;
            
            for trace_i = 1 : 3
                
                idx = find(analysis_window{type_i}{pair_i}{trial_i}, 1, 'last');
                
                t = base_t{type_i}{pair_i}{trial_i}(idx + (-end_pad:end_pad));
                y = all_traces{trace_i}{type_i}{pair_i}{trial_i}(idx + (-end_pad:end_pad)) - (trace_i-1)*trace_offset;
                trace_plot(h_ax(end), t, y, [], []);
                
                t = base_t{type_i}{pair_i}{trial_i}(idx + (-end_pad:0));
                y = all_traces{trace_i}{type_i}{pair_i}{trial_i}(idx + (-end_pad:0)) - (trace_i-1)*trace_offset;
                
                plot(t, y, 'color', aw_col);
                
                title(sprintf('Trial %i, %s', ...
                    trial_id{type_i}{pair_i}{trial_i}, ...
                    protocol_type{type_i}{pair_i}{trial_i}));
                
            end
            
            yL = max(yL, max(get(gca, 'ylim')));
            
            
            
            for ax_i = 1 : length(h_ax)
                
                set(h_ax(ax_i), 'ylim', [-30, yL]);
                
            end
            
            
            FigureTitle(h_fig, sprintf('%s, %i', ...
                probe_fnames{probe_i}, ...
                trial_id{type_i}{1}{trial_i}));
            
            figs.save_fig_to_join();
            
        end
    end
    
    
    figs.join_figs( ...
        sprintf('%s_all_trials_analysis_window.pdf', ...
        probe_fnames{probe_i}));
    figs.clear_figs();
    
end






