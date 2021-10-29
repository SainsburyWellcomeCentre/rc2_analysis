experiment_groups    = {'visual_flow', 'darkness', 'mismatch_nov20', 'mismatch_jul21'};

ctl                  = RC2Analysis();
probe_ids            = ctl.get_probe_ids(experiment_groups{:});


for ii = 1 : length(probe_ids)
    
    % load formatted data
    ctl.create_svm_table(probe_ids{ii});
%     ctl.create_replay_offsets_table(probe_ids{ii});
end
