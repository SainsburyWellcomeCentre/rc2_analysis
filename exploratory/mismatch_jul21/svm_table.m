experiment_group     = 'mismatch_jul21';


ctl                 = RC2Analysis();
probe_ids           = ctl.get_probe_ids(experiment_group);


for i = 1 : length(probe_ids)
    
    data            = ctl.load_formatted_data(probe_ids{i});
    exp_obj         = MismatchJul21Experiment(data, config);
    clusters        = data.selected_clusters;
    
    svm             = StationaryVsMotionTable(config, probe_ids{i});
    
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
