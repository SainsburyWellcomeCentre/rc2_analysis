% match trials in Coupled and EncoderOnly to the replayed trials to get
% offset... this takes a long time


% get details of the experiment
experiment_name         = 'mismatch_jul21';
trial_types             = {'EncoderOnly'};
start_at_recording      = 2;
save_figures            = true;



%%
probe_fnames    = experiment_details(experiment_name);
config          = config_rc2_analysis();
loader          = Loader(config);

figs            = RC2Figures(config);
figs.save_on    = save_figures;
figs.set_figure_subdir(experiment_name, 'matched_trials_overlay');



for probe_i = start_at_recording : length(probe_fnames)
    
    data                = loader.formatted_data(probe_fnames{probe_i});
    exp_obj             = MismatchJul2021Session(data, config);
    
    table_row           = 0;
    offset_table        = table();
    
    for type_i = 1 : length(trial_types)
        
        these_trials    = exp_obj.trials_of_type(trial_types{type_i});
        matched_trial_n = 0;
        
        for trial_i = 1 : length(these_trials)
            
            repeats = exp_obj.get_replay_of(these_trials(trial_i));
            
            for rep_i = 1 : length(repeats)
                
                matched_trial_n = matched_trial_n + 1;
                sp_n = mod(matched_trial_n - 1, 10) + 1;
                
                if sp_n == 1
                    just_saved = 0;
                    h_fig = figs.a4figure();
                end
                
                offset = exp_obj.match_replay_to_original(these_trials(trial_i), repeats(rep_i));
                
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
    
    csv_fname = fullfile(config.summary_data_dir, 'trial_matched_offsets', sprintf('%s_trial_offset_match.csv', probe_fnames{probe_i}));
    writetable(offset_table, csv_fname);
    
    figs.join_figs(sprintf('%s_trials_matched.pdf', probe_fnames{probe_i}))
    figs.clear_figs();
end
