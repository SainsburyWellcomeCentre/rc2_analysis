% Create tuning curve tables for velocity. 
% These tables are required by any script calling `load_tuning_curves`.

experiment_groups    = {'ambient_light'};
trial_group_labels   = {'RT'};

ctl                  = RC2Analysis();

for jj = 1 : length(trial_group_labels)
    probe_ids            = ctl.get_probe_ids(experiment_groups{jj});

    for ii = 1 : length(probe_ids)
        ii
        ctl.create_tuning_curves(probe_ids{ii}, trial_group_labels);
    end
end