% takes about 4 minutes to load all the data waveforms
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids('visual_flow', 'darkness', 'mismatch_nov20', 'mismatch_jul21');

for ii = 1 : length(probe_ids)
    
    fprintf('%i/%i\n', ii, length(probe_ids));
        
    data(ii) = ctl.load_formatted_data(probe_ids{ii});
end
