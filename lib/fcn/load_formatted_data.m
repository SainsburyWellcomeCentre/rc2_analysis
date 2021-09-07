function data = load_formatted_data(probe_recording, config)

formatted_data_fname = fullfile(config.formatted_data_dir, [probe_recording, '.mat']);
data = load_data(formatted_data_fname);
data.probe_recording = probe_recording;
data.experiment_group = get_experiment_group(data.sessions(1).id, config);
data = DataController(data, config);



function exp_type = get_experiment_group(session_id, config)

rec_table = readtable(fullfile(config.local_data_dir, 'summary_data', 'session_list.csv'));
idx = strcmp(rec_table.session_id, session_id);
exp_type = rec_table.short_name{idx};
