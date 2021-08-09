close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'visual_flow';
combination = 'protocols';

config = RC2AnalysisConfig();
figs = RC2Figures(config);

% where to save figure
figs.save_on = true;
figs.set_figure_subdir('visual_flow', 'all_trials_by_replay');

% get details of the experiment
probe_fnames = experiment_details(experiment, combination);

trial_types = {'Coupled', 'EncoderOnly'};

options = default_options();


for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    vf = VisualFlowExperiment(data, config);
    
    base_t = {};
    M_trace = {};
    V_trace = {};
    T_trace = {};
    trial_id = {};
    bouts = {};
    
    base_t_rep = {};
    M_trace_rep = {};
    V_trace_rep = {};
    T_trace_rep = {};
    trial_id_rep = {};
    bouts_rep = {};
    
    for type_i = 1 : length(trial_types)
        
        these_trials = vf.trials_of_type(trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            
            base_t{type_i}{trial_i} = these_trials(trial_i).probe_t - ...
                these_trials(trial_i).probe_t(1);
            M_trace{type_i}{trial_i} = these_trials(trial_i).filtered_teensy;
            V_trace{type_i}{trial_i} = these_trials(trial_i).multiplexer_output;
            T_trace{type_i}{trial_i} = these_trials(trial_i).stage;
            trial_id{type_i}{trial_i} = these_trials(trial_i).id;
            
            bouts_t = get_motion_bouts_by_trial(these_trials(trial_i), options.stationary);
            bouts_t = bouts_t([bouts_t(:).duration] > 2);
            
            bouts{type_i}{trial_i} = bouts_t;
            
            % find repeats of this trial
            repeats = vf.get_replay_of(these_trials(trial_i));
            
            for rep_i = 1 : length(repeats)
                
                base_t_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).probe_t - ...
                    repeats(rep_i).probe_t(1);
                M_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).filtered_teensy;
                V_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).multiplexer_output;
                T_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).stage;
                trial_id_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).id;
                
                bouts_t = get_motion_bouts_by_trial(repeats(rep_i), options.stationary);
                bouts_t = bouts_t([bouts_t(:).duration] > 2);
                bouts_rep{type_i}{trial_i}{rep_i} = bouts_t;
            end
        end
    end
    
    all_traces = {M_trace, V_trace, T_trace};
    all_traces_rep = {M_trace_rep, V_trace_rep, T_trace_rep};
    
    %% plot each trial
    for type_i = 1 : length(base_t)
        
        for trial_i = 1 : length(base_t{type_i})
            
            section_n = 1;
            sp_n = (ceil(section_n/2) - 1)*8 + mod(section_n-1, 2) + 1;
            
            if section_n == 1
                h_fig = figs.a4figure();
            end
            
            h_ax = [];
            
            yL = 0;
            for trace_i = 1 : length(all_traces) 
                h_ax(end+1) = subplot(8, 2, sp_n + (trace_i-1)*2);
                trace_plot(h_ax(end), base_t{type_i}{trial_i}, ...
                    all_traces{trace_i}{type_i}{trial_i}, ...
                    bouts{type_i}{trial_i}, ...
                    trial_id{type_i}{trial_i});
                yL = max(yL, max(get(h_ax(end), 'ylim')));
            end
            
            for rep_i = 1 : length(base_t_rep{type_i}{trial_i})
                
                section_n = section_n + 1;
                sp_n = (ceil(section_n/2) - 1)*8 + mod(section_n-1, 2) + 1;
                
                for trace_i = 1 : length(all_traces_rep)
                    
                    h_ax(end+1) = subplot(8, 2, sp_n + (trace_i - 1)*2);
                    
                    trace_plot(h_ax(end), base_t_rep{type_i}{trial_i}{rep_i}, ...
                        all_traces_rep{trace_i}{type_i}{trial_i}{rep_i}, ...
                        bouts_rep{type_i}{trial_i}{rep_i}, ...
                        trial_id_rep{type_i}{trial_i}{rep_i});
                    
                    yL = max(yL, max(get(h_ax(end), 'ylim')));
                end
            end
            
            for ax_i = 1 : length(h_ax)
                set(h_ax(ax_i), 'ylim', [-5, yL]);
            end
            
            
            FigureTitle(h_fig, sprintf('%s, %s, Trial %i', probe_fnames{probe_i}, trial_types{type_i}, trial_id{type_i}{trial_i}));
            figs.save_fig_to_join();
            
        end
        
        figs.join_figs(sprintf('%s_%s_trials_by_replay.pdf', probe_fnames{probe_i}, trial_types{type_i}));
        figs.clear_figs();
    end
end


