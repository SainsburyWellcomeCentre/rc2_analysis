% Update selected clusters (as in the text file) by recalculating formatted data.

ctl = RC2Analysis();
preprocess = RC2Preprocess();
format = RC2Format();

probe_ids  = ctl.get_probe_ids('visual_flow', ...
                               'mismatch_nov20', ...
                               'mismatch_jul21', ...
                               'mismatch_darkness_oct21', ...
                               'darkness');


for ii = 1 : length(probe_ids)
    fprintf('Processing %s (%i/%i)\n', probe_ids{ii}, ii, length(probe_ids));
    %preprocess.create_selected_clusters_txt(probe_ids{ii});
    format.format(probe_ids{ii});
end
    