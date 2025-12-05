% Create tuning curve tables for acceleration. 
% These tables are required by any script calling `load_tuning_curves_acceleration`


experiment_groups    = {'ambient_light'};
trial_group_labels   = {'RT'};

ctl                  = RC2Analysis();

for jj = 1 : length(experiment_groups)
    probe_ids            = ctl.get_probe_ids(experiment_groups{jj});

    for ii = 1 : length(probe_ids)
        ii
        ctl.create_tuning_curves_acceleration(probe_ids{ii}, trial_group_labels);
    end
end