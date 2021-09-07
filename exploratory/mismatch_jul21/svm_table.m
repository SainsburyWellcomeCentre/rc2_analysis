experiment_name     = 'mismatch_jul21';

probe_fnames        = experiment_details(experiment_name, 'protocols');
config              = config_rc2_analysis();
loader              = Loader(config);


for probe_i = 1 : length(probe_fnames)
    
    data            = loader.formatted_data(probe_fnames{probe_i});
    exp_obj         = MismatchJul21Experiment(data, config);
    clusters        = data.selected_clusters;
    
    svm             = StationaryVsMotionTable(config, probe_fnames{probe_i});
    
    trials          = exp_obj.trials;
    
    for trial_i = 1 : length(trials)
        
        aligned_trial = exp_obj.to_aligned(trials(trial_i));
        
        svm.add_trial(aligned_trial);
        
        for cluster_i = 1 : length(clusters)
            svm.add_table_row_for_cluster(clusters(cluster_i));
        end
    end
    
    svm.save_table();
end
