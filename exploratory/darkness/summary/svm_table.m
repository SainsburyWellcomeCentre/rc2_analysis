config              = RC2AnalysisConfig();
probe_fnames        = experiment_details('darkness', 'protocols');
main_trial_types    = {'Coupled', 'EncoderOnly', 'StageOnly'};

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    dark            = DarknessExperiment(data, config);
    clusters        = data.selected_clusters;
    
    svm             = StationaryVsMotionTable(config, probe_fnames{probe_i});
    
    for type_i = 1 : length(main_trial_types)
        
        these_trials = dark.trials_of_type(main_trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            
            this_trial = these_trials(trial_i);
            
            replayed_trial = dark.get_replayed_trial(this_trial);
            
            if ~isempty(replayed_trial)
                offset = dark.get_offset(replayed_trial, this_trial);
                this_trial = AlignedTrial(this_trial, replayed_trial, offset);
            end
            
            svm.add_trial(this_trial);
            
            for cluster_i = 1 : length(clusters)
                svm.add_table_row_for_cluster(clusters(cluster_i));
            end
        end
    end
    
    svm.save_table();
end
