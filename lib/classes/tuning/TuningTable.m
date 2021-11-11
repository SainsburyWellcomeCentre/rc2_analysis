classdef TuningTable < handle
    
    properties
        
        probe_id
        
        trials
        trial_group_labels
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
        
        
        
        function add_trials(obj, trials, trial_group_labels)
            
            obj.trial_group_labels = trial_group_labels;
            obj.trials = trials;
            obj.velocity_bins = VelocityBins(trials);
            obj.vtc = VelocityTuningCurve(trials, obj.velocity_bins);
        end
        
        
        
        function curve_info = tuning_curve(obj, cluster)
            
            [tuning, timing, stat_rate, stat_time] = obj.vtc.fr_curve(cluster);
            
            shuff = ShuffleTuning(tuning, obj.velocity_bins.bin_centers);
            
            curve_info.probe_id = obj.probe_id;
            curve_info.cluster_id = cluster.id;
            curve_info.trial_group_labels = obj.trial_group_labels;
            curve_info.trial_ids = cellfun(@(x)(x.trial_id), obj.trials);
            curve_info.tuning = tuning;
            curve_info.timing = timing;
            curve_info.stationary_fr = stat_rate;
            curve_info.stationary_time = stat_time;
            curve_info.bin_edges = obj.velocity_bins.bin_edges;
            curve_info.bin_centers = obj.velocity_bins.bin_centers;
            curve_info.prc_per_bin = obj.velocity_bins.prc_per_bin;
            curve_info.shuffled = shuff.get_summary();
            
            
%             warning('off', 'MATLAB:table:RowsAddedExistingVars');
%             
%             % fill the table
%             for ii = 1 : length(obj.trials)
%                 
%                 table_row = obj.n_rows + 1;
%                 
%                 obj.tbl.probe_name{table_row} = obj.probe_id;
%                 obj.tbl.cluster_id(table_row) = cluster.id;
%                 obj.tbl.cluster_region{table_row} = cluster.region_str;
%                 obj.tbl.cluster_depth(table_row) = cluster.depth;
%                 obj.tbl.cluster_from_tip(table_row) = cluster.distance_from_probe_tip;
%                 
%                 obj.tbl.binned_trial_ids{table_row} = cellfun(@(x)(x.trial_id), obj.trials);
%                 obj.tbl.trial_id(table_row) = obj.trials{ii}.trial_id;
%                 obj.tbl.trial_group_label{table_row} = obj.trials{ii}.trial_group_label;
%                 obj.tbl.tuning{table_row} = tuning(:, ii);
%                 obj.tbl.timing{table_row} = timing(:, ii);
%                 obj.tbl.bin_edges{table_row} = obj.velocity_bins.bin_edges;
%                 obj.tbl.stationary_fr{table_row} = stat_rate(ii);
%                 obj.tbl.stationary_time{table_row} = stat_time(ii);
%             end
%             
%             warning('on', 'MATLAB:table:RowsAddedExistingVars');
        end
    end
end
