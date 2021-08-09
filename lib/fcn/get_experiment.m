function exp_obj = get_experiment(data, config)

if strcmp(data.experiment_type, 'visual_flow')
    exp_obj         = VisualFlowExperiment(data, config);
elseif strcmp(data.experiment_type, 'darkness')
    exp_obj         = DarknessExperiment(data, config);
elseif strcmp(data.experiment_type, 'mismatch_nov20')
    exp_obj         = MismatchExperiment(data, config);
elseif strcmp(data.experiment_type, 'passive')
    exp_obj         = PassiveExperiment(data, config);
elseif strcmp(data.experiment_type, 'head_tilt')
    exp_obj         = HeadTiltExperiment(data, config);
end
