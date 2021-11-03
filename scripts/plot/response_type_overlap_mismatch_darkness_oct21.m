% 
experiment_group     = 'mismatch_darkness_oct21';

ctl                 = RC2Analysis();
probe_ids           = ctl.get_probe_ids(experiment_group);

probe_ids_store     = {};
cluster_ids_store   = [];
R_response          = [];
T_response          = [];
RT_MM_response      = [];

for i = 1 : length(probe_ids)
    
    % load formatted data
    data = ctl.load_formatted_data(probe_ids{i});
    
    % VISp clusters
    cluster_ids = data.selected_cluster_ids();
    
    for j = 1 : length(cluster_ids)
        
        probe_ids_store{end+1} = i;
        cluster_ids_store(end+1) = cluster_ids(j);
        [~, ~, R_response(end+1)] = data.is_stationary_vs_motion_significant(cluster_ids(j), 'R');
        [~, ~, T_response(end+1)] = data.is_stationary_vs_motion_significant(cluster_ids(j), 'T');
        [~, ~, RT_MM_response(end+1)] = data.is_mismatch_significant(cluster_ids(j), 'RT_gain_up');
    end
end



%% save csv
csvs                    = CSVManager(ctl);
csvs.save_on            = true;
csvs.set_csv_fulldir('C:\Users\lee\Desktop');
csvs.create_table(  'probe_id',          probe_ids_store(:), ...
                    'cluster_id',        cluster_ids_store(:), ...
                    'R significant',     R_response(:), ...
                    'T_signficant',      T_response(:), ...
                    'RT_MM_signficant',  RT_MM_response(:));
csvs.save('R_T_MM_overlap');
