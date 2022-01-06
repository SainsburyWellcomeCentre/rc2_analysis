classdef StationaryVsMotionTable < handle
% StationaryVsMotionTable Class for handling the table of stationary and
% motion firing rates for trials/clusters
%
%   StationaryVsMotionTable Properties:
%       n_rows                  - number of rows currently in the table
%       svm_table               - the table
%       current_trial           - an instance of class Trial
%       current_stationary_mask - a boolean mask associated with the Trial
%                                 indicating the stationary periods
%       current_motion_mask     - a boolean mask associated with the Trial
%                                 indicating the motion periods
%       current_stationary_time - total time in stationary for the Trial
%       current_motion_time     - total time in motion for the Trial
%
%   StationaryVsMotionTable Methods:
%       add_trial               - adds a trial and fills in the masks
%       add_table_row_for_cluster - adds information about firing rates for
%                                   a cluster to the table
%
%   After object creation, a trial is added with the `add_trial` method.
%   Next several clusters can be added iterative using the
%   `add_table_row_for_cluster` which adds a row to the table for a
%   cluster.
%   
%   This approach removes the need to compute the stationary and motion
%   times for each cluster.
%
%
%   See also: FormattedData.create_svm_table

    properties
        
        n_rows = 0
        svm_table = table([]);
        
        current_trial
        current_stationary_mask
        current_motion_mask
        current_stationary_time
        current_motion_time
    end
    
    
    methods
        
        function obj = StationaryVsMotionTable()
        % StationaryVsMotionTable
        %
        %   StationaryVsMotionTable() creates the object, but doesn't do
        %   anything.
        
        end
        
        
        
        function val = get.n_rows(obj)
        %%return number of rows currently in the table
            val = size(obj.svm_table, 1);
        end
        
        
        
        function add_trial(obj, trial)
        %%add_trial Adds a trial and fills in the masks
        %
        %   add_trial(TRIAL) takes a trial, TRIAL (object of class Trial)
        %   and computes and stores masks for the stationary and motion
        %   periods for that trial.
        
            obj.current_trial = trial;
            obj.current_stationary_mask = trial.stationary_mask;
            obj.current_motion_mask = trial.motion_mask;
            obj.current_stationary_time = trial.stationary_time;
            obj.current_motion_time = trial.motion_time;
        end
        
        
        
        function add_table_row_for_cluster(obj, cluster)
        %%add_table_row_for_cluster Adds information about firing rates for
        %%a cluster to the table 
        %
        %   add_table_row_for_cluster(CLUSTER) takes a cluster, CLUSTER
        %   (object of class Cluster) and computes the firing rates during
        %   the stationary and motion periods of the trial in property
        %   `current_trial`. It then stores this information as a row in
        %   the main `svm_table`.
        %
        %   See also: add_trial
        
            % get convolved firing rate
            fr = FiringRate(cluster.spike_times);
            fr_conv = fr.get_convolution(obj.current_trial.probe_t);
            
            stationary_rate = mean(fr_conv(obj.current_stationary_mask));
            motion_rate = mean(fr_conv(obj.current_motion_mask));
            
            table_row = obj.n_rows + 1;
            
            % fill the table
            obj.svm_table.probe_id{table_row} = obj.probe_id;
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
    end
end