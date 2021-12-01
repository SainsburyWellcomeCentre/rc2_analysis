experiment_groups    = {'visual_flow', 'darkness', 'mismatch_nov20', 'mismatch_jul21', 'mismatch_darkness_oct21'};

ctl               = RC2Preprocess();
probe_ids         = ctl.get_probe_ids(experiment_groups{:});


for ii = 1 : length(probe_ids)
    
    ctl.process_camera_data(probe_ids{ii});
end
