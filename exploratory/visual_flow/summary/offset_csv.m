% match trials in Coupled and EncoderOnly to the replayed trials to get
% offset... this takes a long time

% general config info
config = RC2AnalysisConfig();

% where to save figure
figs = RC2Figures(config);
figs.save_on = true;
figs.set_figure_subdir('visual_flow', 'matched_trials_overlay');

% get details of the experiment
probe_fnames = experiment_details('visual_flow', 'protocols');

trial_types = {'Coupled', 'EncoderOnly'};



for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    vf = VisualFlowExperiment(data, config);
    
    table_row = 0;
    offset_table = table();
    
    for type_i = 1 : length(trial_types)
        
        these_trials = vf.trials_of_type(trial_types{type_i});
        
        matched_trial_n = 0;
        
        for trial_i = 1 : length(these_trials)
            
            repeats = vf.get_replay_of(these_trials(trial_i));
            
            for rep_i = 1 : length(repeats)
                
                matched_trial_n = matched_trial_n + 1;
                sp_n = mod(matched_trial_n - 1, 10) + 1;
                
                if sp_n == 1
                    just_saved = 0;
                    h_fig = figs.a4figure();
                end
                
                offset = vf.match_replay_to_original(these_trials(trial_i), repeats(rep_i));
                
                table_row = table_row + 1;
                offset_table.trial_id(table_row) = these_trials(trial_i).id;
                offset_table.replay_trial_id(table_row) = repeats(rep_i).id;
                offset_table.offset(table_row) = offset;
                
                subplot(10, 1, sp_n);
                hold on;
                plot(these_trials(trial_i).multiplexer_output);
                plot(repeats(rep_i).multiplexer_output(offset:end));
                set(gca, 'plotboxaspectratio', [3, 1, 1]);
                box off;
                title(sprintf('Trial %i, replay trial %i', these_trials(trial_i).id, repeats(rep_i).id));
                
                if sp_n == 10
                    figs.save_fig_to_join();
                    just_saved = 1;
                end
            end
            
            if trial_i == length(these_trials) && ~just_saved
                figs.save_fig_to_join();
            end
        end
    end
    
    csv_fname = fullfile(config.summary_data_dir, 'match_visual_flow_trials', sprintf('%s_trial_offset_match.csv', probe_fnames{probe_i}));
    writetable(offset_table, csv_fname);
    
    figs.join_figs(sprintf('%s_trials_matched.pdf', probe_fnames{probe_i}))
    figs.clear_figs();
end
