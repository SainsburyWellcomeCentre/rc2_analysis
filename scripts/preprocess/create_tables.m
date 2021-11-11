experiment_groups    = {'darkness'};%, 'mismatch_nov20', 'mismatch_jul21', 'mismatch_darkness_oct21'};

ctl                  = RC2Analysis();
probe_ids            = ctl.get_probe_ids(experiment_groups{:});

for ii = 1 : length(probe_ids)
    ii
    % load formatted data
%     ctl.create_replay_offsets_table(probe_ids{ii});
%     ctl.create_svm_table(probe_ids{ii});
    ctl.create_tuning_curves(probe_ids{ii}, {'RT', 'R', {'T_bank', 'T_RT', 'T_R'}});
end
