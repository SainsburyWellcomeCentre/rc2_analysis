config              = RC2AnalysisConfig();
probe_fnames        = experiment_details('visual_flow', 'protocols');
main_trial_types    = {'Coupled', 'EncoderOnly'};

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    vf              = VisualFlowExperiment(data, config);
    clusters        = data.selected_clusters;
    
    svm             = StationaryVsMotionTable(config, probe_fnames{probe_i});
    
    for type_i = 1 : length(main_trial_types)
        
        these_trials = vf.trials_of_type(main_trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            trial_i
            svm.add_trial(these_trials(trial_i));
            
            for cluster_i = 1 : length(clusters)
                svm.add_table_row_for_cluster(clusters(cluster_i));
            end
            
            repeats = vf.get_replay_of(these_trials(trial_i));
            
            for rep_i = 1 : length(repeats)
                
                offset = vf.get_offset(these_trials(trial_i), repeats(rep_i));
                aligned_trial = AlignedTrial(repeats(rep_i), these_trials(trial_i), offset);
                
                svm.add_trial(aligned_trial);
                
                for cluster_i = 1 : length(clusters)
                    svm.add_table_row_for_cluster(clusters(cluster_i));
                end
            end
        end
        
    end
    
    svm.save_table();
end
