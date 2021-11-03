classdef TuningTable < handle
    
    properties
        
        probe_id
        
        trials
        velocity_bins
        vtc
        
        n_rows
        tbl = table([]);
    end
    
    
    
    methods
        
        function obj = TuningTable(probe_id)
            
            obj.probe_id = probe_id;
        end
        
        
        
        function val = get.n_rows(obj)
        
            val = size(obj.tbl, 1);            
        end
        
        
        
        function add_trials(obj, trials)
            
            obj.trials = trials;
            obj.velocity_bins = VelocityBins(trials);
            obj.vtc = VelocityTuningCurve(trials, obj.velocity_bins);
        end
        
        
        
        function add_row_for_cluster(obj, cluster)
            
            [tuning, timing, stat_rate, stat_time] = obj.vtc.fr_curve(cluster);
            
            table_row = obj.n_rows + 1;
            
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            
            % fill the table
            obj.tbl.probe_name{table_row} = obj.probe_id;
            obj.tbl.cluster_id(table_row) = cluster.id;
            obj.tbl.cluster_region{table_row} = cluster.region_str;
            obj.tbl.cluster_depth(table_row) = cluster.depth;
            obj.tbl.cluster_from_tip(table_row) = cluster.distance_from_probe_tip;
            
            obj.tbl.trial_ids{table_row} = cellfun(@(x)(x.trial_id), obj.trials);
            obj.tbl.tuning{table_row} = tuning;
            obj.tbl.timing{table_row} = timing;
            obj.tbl.bin_edges{table_row} = obj.velocity_bins.bin_edges;
            obj.tbl.stationary_fr{table_row} = stat_rate;
            obj.tbl.stationary_time{table_row} = stat_time;
            
            warning('on', 'MATLAB:table:RowsAddedExistingVars');
        end
    end
end
