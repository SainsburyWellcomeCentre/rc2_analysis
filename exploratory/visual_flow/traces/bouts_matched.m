% function plot_all_visual_flow_bouts_across_conditions()

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
figs.set_figure_subdir('visual_flow', 'all_bouts_matched_to_replays');

% get details of the experiment
probe_fnames = experiment_details(experiment, combination);

trial_types = {'Coupled', 'EncoderOnly'};
fname_suffix = {'MVT', 'MV'};

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
    protocol_type_rep = {};
    
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
            
            % match
            for rep_i = 1 : length(repeats)
                
                base_t_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).probe_t - ...
                    repeats(rep_i).probe_t(1);
                M_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).filtered_teensy;
                V_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).multiplexer_output;
                T_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).stage;
                trial_id_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).id;
                protocol_type_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).protocol;
                
                offset = vf.get_offset(these_trials(trial_i), repeats(rep_i));
                bout_t = [];
                for bout_i = 1 : length(bouts{type_i}{trial_i})
                    bout_t(bout_i).start_idx = bouts{type_i}{trial_i}(bout_i).start_idx + offset;
                    bout_t(bout_i).end_idx = bouts{type_i}{trial_i}(bout_i).end_idx + offset;
                end
                
                bouts_rep{type_i}{trial_i}{rep_i} = bout_t;
                
            end
        end
    end
    
    all_traces = {M_trace, V_trace, T_trace};
    all_traces_rep = {M_trace_rep, V_trace_rep, T_trace_rep};
    axis_labels = {'M', 'V', 'T'};
    
    %% plot each trial
    for type_i = 1 : length(base_t)
        
        bout_n = 0;
        
        for trial_i = 1 : length(base_t{type_i})
            
            for bout_i = 1 : length(bouts{type_i}{trial_i})
                
                bout_n = bout_n + 1;
                
                if bout_n == 3 || bout_n == 7
                    disp('')
                end
                
                section_n = 1;
                sp_n = (ceil(section_n/2) - 1)*8 + mod(section_n-1, 2) + 1;
                
                if section_n == 1
                    
                    h_fig = figs.a4figure();
                    
                end
                
                h_ax = [];
                
                yL = 0;
                xL = [inf, -inf];
                
                for trace_i = 1 : length(all_traces)
                    
                    h_ax(end+1) = subplot(8, 2, sp_n + (trace_i-1)*2);
                    
                    bout_plot(h_ax(end), base_t{type_i}{trial_i}, ...
                        all_traces{trace_i}{type_i}{trial_i}, ...
                        bouts{type_i}{trial_i}(bout_i));
                    
                    if trace_i == 1
                        title(sprintf('Bout %i, Trial %i, %s', ...
                                      bout_n, ...
                                      trial_id{type_i}{trial_i}, ...
                                      fname_suffix{type_i}));
                    end
                    
                    yL = max(yL, max(get(h_ax(end), 'ylim')));
                    xL(1) = min(xL(1), min(get(h_ax(end), 'xlim')));
                    xL(2) = max(xL(2), max(get(h_ax(end), 'xlim')));
                    
                end
                
                for rep_i = 1 : length(base_t_rep{type_i}{trial_i})
                    
                    if strcmp(protocol_type_rep{type_i}{trial_i}{rep_i}, 'StageOnly')
                    
                        section_n = 2;
                        title_str = 'VT';
                    
                    else
                        
                        section_n = 3;
                        title_str = 'V';
                        
                    end
                    
                    sp_n = (ceil(section_n/2) - 1)*8 + mod(section_n-1, 2) + 1;
                    
                    for trace_i = 1 : length(all_traces_rep)
                        
                        h_ax(end+1) = subplot(8, 2, sp_n + (trace_i - 1)*2);
                        
                        bout_plot(h_ax(end), base_t_rep{type_i}{trial_i}{rep_i}, ...
                            all_traces_rep{trace_i}{type_i}{trial_i}{rep_i}, ...
                            bouts_rep{type_i}{trial_i}{rep_i}(bout_i));
                        
                        if trace_i == 1
                            title(sprintf('Bout %i, Trial %i, %s (%s)', ...
                                          bout_n, ...
                                          trial_id_rep{type_i}{trial_i}{rep_i}, ...
                                          title_str, ...
                                          fname_suffix{type_i}));
                        end
                        
                        yL = max(yL, max(get(h_ax(end), 'ylim')));
                        xL(1) = min(xL(1), min(get(h_ax(end), 'xlim')));
                        xL(2) = max(xL(2), max(get(h_ax(end), 'xlim')));
                        
                    end
                    
                end
                
                for ax_i = 1 : length(h_ax)
                    
                    set(gcf, 'currentaxes', h_ax(ax_i));
                    
                    set(h_ax(ax_i), 'ylim', [-5, yL], 'xlim', xL);
                    
                    text(xL(1), yL, ...
                         axis_labels{mod(ax_i-1, 3)+1}, ...
                         'fontweight', 'bold', ...
                         'fontsize', 14, ...
                         'horizontalalignment', 'left', ...
                         'verticalalignment', 'top');
                    
                end
                
                FigureTitle(h_fig, sprintf('%s, %s, Trial %i, Bout %i', ...
                    probe_fnames{probe_i}, ...
                    trial_types{type_i}, ...
                    trial_id{type_i}{trial_i}, ...
                    bout_i));
                
                figs.save_fig_to_join();
                
            end
            
        end
        
        figs.join_figs(sprintf('%s_%s_bouts_and_replays.pdf', probe_fnames{probe_i}, fname_suffix{type_i}));
        figs.clear_figs();
    end
end


