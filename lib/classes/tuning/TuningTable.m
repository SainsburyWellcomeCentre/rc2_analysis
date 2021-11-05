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
            
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            
            % fill the table
            for ii = 1 : length(obj.trials)
                
                table_row = obj.n_rows + 1;
                
                obj.tbl.probe_name{table_row} = obj.probe_id;
                obj.tbl.cluster_id(table_row) = cluster.id;
                obj.tbl.cluster_region{table_row} = cluster.region_str;
                obj.tbl.cluster_depth(table_row) = cluster.depth;
                obj.tbl.cluster_from_tip(table_row) = cluster.distance_from_probe_tip;
                
                obj.tbl.binned_trial_ids{table_row} = cellfun(@(x)(x.trial_id), obj.trials);
                obj.tbl.trial_id(table_row) = obj.trials{ii}.trial_id;
                obj.tbl.trial_group_label{table_row} = obj.trials{ii}.trial_group_label;
                obj.tbl.tuning{table_row} = tuning(:, ii);
                obj.tbl.timing{table_row} = timing(:, ii);
                obj.tbl.bin_edges{table_row} = obj.velocity_bins.bin_edges;
                obj.tbl.stationary_fr{table_row} = stat_rate(ii);
                obj.tbl.stationary_time{table_row} = stat_time(ii);
            end
            
            warning('on', 'MATLAB:table:RowsAddedExistingVars');
        end
    end
end
