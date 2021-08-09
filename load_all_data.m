% takes about 4 minutes to load all the data waveforms
config          = config_rc2_analysis();
loader          = Loader(config);

recording_ids   = [];
recording_ids  = experiment_details('visual_flow', 'protocols');
recording_ids   = [recording_ids, experiment_details('mismatch_nov20', 'protocols')];
recording_ids  = [recording_ids, experiment_details('darkness', 'protocols')];

recording_id    = [];
cluster_id      = [];
firing_rates    = [];
durations       = [];

for rec_i = 1 : length(recording_ids)
    
    fprintf('%i/%i\n', rec_i, length(recording_ids));
        
    data(rec_i) = loader.formatted_data(recording_ids{rec_i});
end
