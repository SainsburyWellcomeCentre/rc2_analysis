classdef RestrictClusters < handle
% RestrictClusters Class for
%
%  RestrictClusters Properties:
%       isi_viol        - ISI violation threshold (default = 0.15)
%       isolation_dist  - isolation distance threshold (default = 20)
%       amp_cutoff      - amplitude cutoff threshold (default = 0.1)
%       max_drift       - max. drift allowed (default = 80)
%       ctl             - instance of RC2Preprocess
%       probe_id        - string with probe recording ID
%
%  RestrictClusters Methods:
%       restrict_metrics_table -  loads the metrics.csv and restricts it to
%       clusters which satisfy certain quality criteria

    properties
        
        isi_viol        = 0.15
        isolation_dist  = 20
        amp_cutoff      = 0.1
        max_drift       = 80
    end
    
    properties (SetAccess = private)
        
        ctl
        probe_id
    end
    
    
    methods
        
        function obj = RestrictClusters(ctl, probe_id)
        %%RestrictClusters
        %
        %   RestrictClusters(CTL, PROBE_ID). CTL is instance of
        %   RC2Preprocess and PROBE_ID is a string with the probe recording
        %   ID.
        
            obj.ctl = ctl;
            obj.probe_id = probe_id;
        end
        
        
        
        function new_tbl = restrict_metrics_table(obj)
        %%restrict_metrics_table Loads the metrics.csv and restricts it to
        % clusters which satisfy certain quality criteria
        %
        %   TABLE = restrict_metrics_table() restricts the metrics.csv to
        %   clusters which satisfy quality of criteria. 
        %       < `isi_viol`
        %       > `isolation_distance`
        %       < `amp_cutoff`
        %       < `max_drift`
        
            metrics = obj.ctl.load.metrics_csv(obj.probe_id);
            cluster_groups = obj.ctl.load.cluster_groups(obj.probe_id);
            
            idx = metrics.isi_viol < obj.isi_viol & ...
                  metrics.isolation_distance > obj.isolation_dist & ...
                  metrics.amplitude_cutoff < obj.amp_cutoff & ...
                  metrics.max_drift < obj.max_drift;
            
            good_idx = false(size(metrics, 1), 1);
              
            for ii = 1 : size(metrics, 1)
                
                clust_idx = find(metrics.cluster_id(ii) == cluster_groups.cluster_id);
                
                if isempty(clust_idx)
                    continue
                end
                
                good_idx(ii) = strcmp(cluster_groups.group(clust_idx), 'good');
            end
            
            idx = idx & good_idx;
            
            new_tbl = metrics(idx, :);
        end
    end
end
