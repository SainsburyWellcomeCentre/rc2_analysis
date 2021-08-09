classdef StationaryVsMotionTable < handle
    
    properties
        
        config
        
        probe_fname
        
        n_rows
        svm_table = table([]);
        
        current_trial
        current_stationary_mask
        current_motion_mask
        current_stationary_time
        current_motion_time
    end
    
    
    methods
        
        function obj = StationaryVsMotionTable(config, probe_fname)
            
            obj.config = config;
            obj.probe_fname = probe_fname;
            
        end
        
        
        function val = get.n_rows(obj)
        
            val = size(obj.svm_table, 1);
            
        end
        
        
        function add_trial(obj, trial)
            
            obj.current_trial = trial;
            obj.current_stationary_mask = trial.stationary_mask;
            obj.current_motion_mask = trial.motion_mask;
            obj.current_stationary_time = trial.stationary_time;
            obj.current_motion_time = trial.motion_time;
            
        end
        
        
        function add_table_row_for_cluster(obj, cluster)
            
            % get convolved firing rate
            fr = FiringRate(cluster.spike_times);
            fr_conv = fr.get_convolution(obj.current_trial.probe_t);
            
            stationary_rate = mean(fr_conv(obj.current_stationary_mask));
            motion_rate = mean(fr_conv(obj.current_motion_mask));
            
            table_row = obj.n_rows + 1;
            
            % fill the table
            obj.svm_table.probe_name{table_row} = obj.probe_fname;
            obj.svm_table.cluster_id(table_row) = cluster.id;
            obj.svm_table.cluster_region{table_row} = cluster.region_str;
            obj.svm_table.cluster_depth(table_row) = cluster.depth;
            obj.svm_table.cluster_from_tip(table_row) = cluster.distance_from_probe_tip;
            obj.svm_table.protocol{table_row} = obj.current_trial.protocol;
            obj.svm_table.replay_of{table_row} = obj.current_trial.replay_of;
            obj.svm_table.trial_id(table_row) = obj.current_trial.id;
            obj.svm_table.stationary_firing_rate(table_row) = stationary_rate;
            obj.svm_table.time_stationary(table_row) = obj.current_stationary_time;
            obj.svm_table.motion_firing_rate(table_row) = motion_rate;
            obj.svm_table.time_motion(table_row) = obj.current_motion_time;
            
            if isfield(obj.current_trial.config, 'enable_vis_stim')
                if isnumeric(obj.current_trial.config.enable_vis_stim)
                    obj.svm_table.vis_stim(table_row) = obj.current_trial.config.enable_vis_stim;
                else
                    obj.svm_table.vis_stim(table_row) = str2double(obj.current_trial.config.enable_vis_stim);
                end
            else
                obj.svm_table.vis_stim(table_row) = nan;
            end
            
        end
        
        
        function save_table(obj)
            
            csv_fname = fullfile(obj.config.summary_data_dir, 'stationary_vs_motion_fr', sprintf('%s.csv', obj.probe_fname));
            writetable(obj.svm_table, csv_fname); 
        end
        
        
        
    end
end