% check camera data is matched to the triggers
experiment_groups = {'visual_flow', 'darkness', 'mismatch_nov20', 'mismatch_jul21'};

ctl         = RC2Analysis();
probe_ids   = ctl.get_probe_ids(experiment_groups{:});

probe_id = {};
session_id = {};
cam_length = [];
trigger_length = [];

for ii = 1 : length(probe_ids)
    
    data = ctl.load_formatted_data(probe_ids{ii});
    
    for jj = 1 : length(data.sessions)
        probe_id{end+1} = probe_ids{ii};
        session_id{end+1} = data.sessions{jj}.session_id;
        cam_length(end+1) = length(data.sessions{jj}.camera0);
        trigger_length(end+1) = length(data.sessions{jj}.camera_t);
    end
end

tbl = table('probe_id', probe_id(:), ...
            'session_id', session_id(:), ...
            'cam_length', length_cam(:), ...
            'trigger_length', trigger_length(:));
        
