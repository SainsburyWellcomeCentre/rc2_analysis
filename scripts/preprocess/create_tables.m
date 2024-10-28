% Create tuning curve tables for velocity. 
% These tables are required by any script calling `load_tuning_curves`.

experiment_groups    = {'darkness', 'mismatch_darkness_oct21'};%, 'mismatch_nov20', 'mismatch_jul21', };
trial_group_labels   = {{'RT', 'R', {'T_bank', 'T_RT', 'T_R'}}, {'R', 'T', 'RT_gain_up'}};

ctl                  = RC2Analysis();

for jj = 1 : length(trial_group_labels)
    probe_ids            = ctl.get_probe_ids(experiment_groups{jj});

    for ii = 1 : length(probe_ids)
        ii
        ctl.create_tuning_curves(probe_ids{ii}, trial_group_labels);
    end
end