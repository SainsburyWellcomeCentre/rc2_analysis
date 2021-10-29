classdef RestrictClusters < handle
    
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
            
            obj.ctl = ctl;
            obj.probe_id = probe_id;
        end
        
        
        
        function new_tbl = restrict_metrics_table(obj)
            
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
