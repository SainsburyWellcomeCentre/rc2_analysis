experiment          = 'visual_flow';
config              = RC2AnalysisConfig();
probe_fnames        = experiment_details(experiment, 'protocols');

if strcmp(experiment, 'visual_flow')
    types               = {'Coupled', 'EncoderOnly', 'StageOnly', 'StageOnly', 'ReplayOnly', 'ReplayOnly'};
    replay_of           = {'', '', 'Coupled', 'EncoderOnly', 'Coupled', 'EncoderOnly'};
elseif strcmp(experiment, 'darkness')
    types               = {'Coupled', 'EncoderOnly', 'StageOnly'};
    replay_of           = {'', '', ''};
end

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    if strcmp(experiment, 'visual_flow')
        exp_obj = VisualFlowExperiment(data, config);
    else
        exp_obj = DarknessExperiment(data, config);
    end
    
    clusters        = data.selected_clusters;
    
    tt             = TuningTable(config, probe_fnames{probe_i});
    
    for type_i = 1 : length(types)
        
        trials = exp_obj.trials_of_type_replay_of_type(types{type_i}, replay_of{type_i});
        
        for trial_i = 1 : length(trials)
            
            atrials(trial_i) = exp_obj.to_aligned(trials(trial_i));
        end
        
        tt.add_trials(atrials, types{type_i}, replay_of{type_i});
        
        for cluster_i = 1 : length(clusters)
            
            tt.add_row_for_cluster(clusters(cluster_i));
        end
    end
    
    tt.save_table();
end
