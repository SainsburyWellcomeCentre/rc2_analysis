experiment_groups    = {'visual_flow', 'darkness', 'mismatch_nov20', 'mismatch_jul21', 'mismatch_darkness_oct21'};

format               = RC2Format();
probe_ids            = format.get_probe_ids(experiment_groups{:});


% for ii = 1 : length(probe_ids)
%     format.format(probe_ids{ii});
% end


for ii = 1 : length(probe_ids)
    ii
    s.anatomy = format.format_anatomy(probe_ids{ii});
    s.clusters = format.format_clusters(probe_ids{ii});
    format.save.append_to_formatted_data(probe_ids{ii}, s);
end
