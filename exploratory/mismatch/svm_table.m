experiment          = 'mismatch_nov20';

config              = RC2AnalysisConfig();
probe_fnames        = experiment_details(experiment, 'protocols');

for probe_i = 1 : length(probe_fnames)
    
    data            = load_formatted_data(probe_fnames{probe_i}, config);
    exp_obj         = MismatchExperiment(data, config);
    clusters        = data.selected_clusters;
    
    svm             = StationaryVsMotionTable(config, probe_fnames{probe_i});
    
    trials          = exp_obj.trials;
    
    for trial_i = 1 : length(trials)
        
        svm.add_trial(trials(trial_i));
        
        for cluster_i = 1 : length(clusters)
            svm.add_table_row_for_cluster(clusters(cluster_i));
        end
    end
    
    svm.save_table();
end
