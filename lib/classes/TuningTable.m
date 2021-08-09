classdef TuningTable < handle
    
    properties
        
        config
        
        probe_fname
        
        vb
        vtc
        
        n_rows
        tt_table = table([]);
        
        current_protocol
        current_replay
    end
    
    
    
    methods
        
        function obj = TuningTable(config, probe_fname)
            
            obj.config = config;
            obj.probe_fname = probe_fname;
        end
        
        
        
        function val = get.n_rows(obj)
        
            val = size(obj.tt_table, 1);            
        end
        
        
        
        function add_trials(obj, trials, prot_name, replay_of)
            
            obj.vb = VelocityBins(trials);
            obj.vtc = VelocityTuningCurve(trials, obj.vb);
            
            obj.current_protocol = prot_name;
            obj.current_replay = replay_of;
        end
        
        
        
        function add_row_for_cluster(obj, cluster)
            
            [tuning, timing, stat_rate, stat_time] = obj.vtc.fr_curve(cluster);
            
            table_row = obj.n_rows + 1;
            
            % fill the table
            obj.tt_table.probe_name{table_row} = obj.probe_fname;
            obj.tt_table.cluster_id(table_row) = cluster.id;
            obj.tt_table.cluster_region{table_row} = cluster.region_str;
            obj.tt_table.cluster_depth(table_row) = cluster.depth;
            obj.tt_table.cluster_from_tip(table_row) = cluster.distance_from_probe_tip;
            
            obj.tt_table.protocol{table_row} = obj.current_protocol;
            obj.tt_table.replay_of{table_row} = obj.current_replay;
            obj.tt_table.trial_ids{table_row} = [obj.vb.trials(:).id];
            obj.tt_table.tuning{table_row} = tuning;
            obj.tt_table.timing{table_row} = timing;
            obj.tt_table.bin_edges{table_row} = obj.vb.bin_edges;
            obj.tt_table.stationary_fr{table_row} = stat_rate;
            obj.tt_table.stationary_time{table_row} = stat_time;
        end
        
        
        
        function save_table(obj)
            
            mat_fname = fullfile(obj.config.summary_data_dir, 'tuning_table', sprintf('%s.mat', obj.probe_fname));
            tuning_table = obj.tt_table; %#ok<NASGU>
            save(mat_fname, 'tuning_table'); 
        end
    end
end
