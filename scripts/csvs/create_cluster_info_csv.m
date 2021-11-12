% create a csv with each cluster

experiment_groups   = {'visual_flow', 'darkness', 'mismatch_nov20', 'mismatch_jul21', 'mismatch_darkness_oct21'};

ctl                 = RC2Analysis();
probe_ids           = ctl.get_probe_ids(experiment_groups{:});
tbl                 = table([]);

r                   = 0;

for ii = 1 : length(probe_ids)
    
    data            = ctl.load_formatted_data(probe_ids{ii});
    clusters        = data.VISp_clusters();
    
    for jj = 1 : length(clusters)
        
        r = r + 1;
        
        tbl.cluster_id(r)                   = clusters(jj).id;
        tbl.probe_id{r}                     = probe_ids{ii};
        tbl.shank_id(r)                     = clusters(jj).cluster.shank_id;
        tbl.region_str{r}                   = clusters(jj).region_str;
        tbl.depth_um(r)                     = clusters(jj).depth;
        tbl.distance_from_probe_tip(r)      = clusters(jj).distance_from_probe_tip;
        tbl.spiking_class{r}                = clusters(jj).spiking_class;
        
        [tbl.region_relative_depth(r), tbl.region_str_backup{r}] = data.get_relative_layer_depth_of_cluster(clusters(jj).id);
    end
end

tbl(:, 1) = [];
writetable(tbl, 'cluster_info.csv');