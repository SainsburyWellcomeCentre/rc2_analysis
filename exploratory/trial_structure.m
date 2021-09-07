experiment          = 'mismatch_jul21';
start_recording     = 2;

config              = config_rc2_analysis();
loader              = Loader(config);

figs                = RC2Figures(config);
figs.save_on        = true;
figs.set_figure_subdir(experiment, 'trial_structure');

probe_fnames        = experiment_details(experiment, 'protocols');

for probe_i = start_recording : length(probe_fnames)
    
    data                = loader.formatted_data(probe_fnames{probe_i});
    clusters            = data.VISp_clusters();
    
    switch experiment
        case 'visual_flow'
            exp_obj = VisualFlowExperiment(data, config);
        case 'darkness'
            exp_obj = DarknessExperiment(data, config);
        case 'passive'
            exp_obj = PassiveExperiment(data, config);
        case 'head_tilt'
            exp_obj = HeadTiltExperiment(data, config);
        case 'mismatch_nov20'
            exp_obj = MismatchExperiment(data, config);
        case 'mismatch_jul21'
            exp_obj = MismatchJul2021Session(data, config);
    end
    
    for trial_i = 1 : length(exp_obj.trials)
        
        this_trial = exp_obj.trials(trial_i);
        
        if this_trial.is_replay && ~strcmp(this_trial.replay_of, 'Bank')
            replayed_trial = exp_obj.get_replayed_trial(this_trial);
            offset = exp_obj.get_offset(replayed_trial, this_trial);
            this_trial = AlignedTrial(this_trial, replayed_trial, offset);
        end
        
        ts = TrialStructure(this_trial, clusters);
        
        switch experiment
            case {'visual_flow', 'darkness', 'mismatch_jul21'}
                if this_trial.is_replay
                    title(ts.h_ax, sprintf('Trial %i, %s;   replay of Trial %i, %s', ...
                                            this_trial.id, ...
                                            this_trial.protocol, ...
                                            this_trial.replayed_trial_id, ...
                                            this_trial.replay_of));
                else
                    title(ts.h_ax, sprintf('Trial %i, %s',  this_trial.id, ...
                                                            this_trial.protocol));
                end
            case {'passive', 'head_tilt', 'mismatch_nov20'}
                id = exp_obj.get_trial_protocol_id(this_trial);
                label = exp_obj.label_from_id(id);
                title(ts.h_ax, sprintf('Trials %i, %s', this_trial.id, label));
        end
        
        figs.save_fig_to_join();
    end
    
    fname = sprintf('%s.pdf', probe_fnames{probe_i});
    figs.join_figs(fname);
    figs.clear_figs();
end
