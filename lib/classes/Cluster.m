classdef Cluster < handle
    
    properties (Hidden = true, SetAccess = private)
        
        cluster
    end
    
    
    properties (SetAccess = private)
        
        fr
    end
    
    
    properties (Dependent = true)
        
        id
        spike_times
        region_str
        depth
        spiking_class
        distance_from_probe_tip;
        duration
        overall_firing_rate
    end
    
    
    
    methods
        
        function obj = Cluster(cluster)
        %%class for interacting a little with a cluster
            obj.cluster = cluster;
            obj.fr = FiringRate(cluster.spike_times);
        end
        
        
        
        function val = get.id(obj)
        %%cluster id
            val = obj.cluster.id;
        end
        
        
        
        function val = get.spike_times(obj)
            val = obj.cluster.spike_times;
        end
        
        
        
        function val = get.region_str(obj)
        %%which region is the cluster in
            val = obj.cluster.region_str;
        end
        
        
        
        function val = get.depth(obj)
            val = obj.cluster.depth;
        end
        
        
        
        function val = get.distance_from_probe_tip(obj)
        %%distance of cluster from probe tip
            val = obj.cluster.distance_from_probe_tip;
        end
           
        
        
        function val = get.duration(obj)
            val = obj.cluster.duration;
        end
        
        
        
        function val = get.overall_firing_rate(obj)
            val = obj.cluster.firing_rate;
        end
        
        
        
        function val = get.spiking_class(obj)
        %%return the spiking class of the cluster
            if obj.cluster.duration < constants('spiking_class_threshold_ms')
                val = 'narrow';
            else
                val = 'wide';
            end
        end
        
        
        
        function val = is_VISp(obj)
        %%return whether cluster is in VISp
            val = ~isempty(regexp(obj.region_str, 'VISp\d', 'once'));
        end
    end
end
