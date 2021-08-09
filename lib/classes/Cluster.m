classdef Cluster < handle
    
    properties (SetAccess = private)
        
        cluster
    end
    
    properties
        
        spike_class_threshold = 0.45
    end
    
    properties (Dependent = true)
        
        spiking_class
    end
    
    
    
    methods
        
        function obj = Cluster(cluster)
            
            obj.cluster = cluster;
        end
        
        
        
        function val = get.spiking_class(obj)
            
            if obj.cluster.duration < obj.spike_class_threshold
                val = 'FS';
            else
                val = 'RS';
            end
        end
    end
end
