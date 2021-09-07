function exp_obj = get_experiment(data, config)

if strcmp(data.experiment_group, 'visual_flow')
    exp_obj         = VisualFlowExperiment(data, config);
elseif strcmp(data.experiment_group, 'darkness')
    exp_obj         = DarknessExperiment(data, config);
elseif strcmp(data.experiment_group, 'mismatch_nov20')
    exp_obj         = MismatchExperiment(data, config);
elseif strcmp(data.experiment_group, 'passive')
    exp_obj         = PassiveExperiment(data, config);
elseif strcmp(data.experiment_group, 'head_tilt')
    exp_obj         = HeadTiltExperiment(data, config);
end
