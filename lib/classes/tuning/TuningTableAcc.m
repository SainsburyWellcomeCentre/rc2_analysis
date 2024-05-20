classdef TuningTableAcc < handle
% TODO
    properties
        
        probe_id
        
        trials
        trial_group_labels
        acceleration_bins
        atc
        
        n_rows
        tbl = table([]);
        
        mode = "all"
    end
    
    
    
    methods
        
        function obj = TuningTableAcc(probe_id)
        % TuningTable
        %
        %   TuningTable(PROBE_ID) creates the object. PROBE_ID is a string
        %   with the probe recording ID.
        
            obj.probe_id = probe_id;
        end
        
        
        
        function val = get.n_rows(obj)
        %%number of rows of the table
            val = size(obj.tbl, 1);            
        end
        
        
        
        function add_trials(obj, trials, trial_group_labels)
        %TODO
        
            obj.trial_group_labels = trial_group_labels;
            obj.trials = trials;
            obj.acceleration_bins = AccelerationBins(trials, obj.mode);
            obj.atc = AccelerationTuningCurve(trials, obj.acceleration_bins);
        end
        
        
        
        function curve_info = tuning_curve(obj, cluster)
        %%TODO
        
            [tuning, timing, stat_rate, stat_time] = obj.atc.fr_curve(cluster);
            
            
            % Remove the last bin - hacky way
            tuning = tuning(2:end-1, :);
            timing = timing(2:end-1, :);
            if length(obj.acceleration_bins.bin_edges) == 21
                obj.acceleration_bins.bin_edges = obj.acceleration_bins.bin_edges(2:end-1);
            end
            % end of hack
            
            shuff = ShuffleTuning(tuning, obj.acceleration_bins.bin_centers);
            
            curve_info.probe_id = obj.probe_id;
            curve_info.cluster_id = cluster.id;
            curve_info.trial_group_labels = obj.trial_group_labels;
            curve_info.trial_ids = cellfun(@(x)(x.trial_id), obj.trials);
            curve_info.tuning = tuning;
            curve_info.timing = timing;
            curve_info.stationary_fr = stat_rate;
            curve_info.stationary_time = stat_time;
            curve_info.bin_edges = obj.acceleration_bins.bin_edges;
            curve_info.bin_centers = obj.acceleration_bins.bin_centers;
            curve_info.prc_per_bin = obj.acceleration_bins.prc_per_bin;
            curve_info.shuffled = shuff.get_summary();
            
            
        end
    end
end
