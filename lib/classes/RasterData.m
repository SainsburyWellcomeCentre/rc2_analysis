classdef RasterData < handle
    
    properties (SetAccess = private)
        
        cluster
    end
    
    properties
        
        limits = [-1, 1]
        trigger_times
    end
    
    properties (Dependent = true)
        
        common_t
        n_clusters
    end
    
    properties (Hidden = true)
        
        fs = 60
    end
    
    
    
    methods
        
        function obj = RasterData(cluster)
        %%collection of methods for handling raster data
            obj.cluster = cluster;
        end
        
        
        
        function val = get.common_t(obj)
        %%common timebase for the convolution
            val = linspace(obj.limits(1), obj.limits(2), round(diff(obj.limits)*obj.fs)+1);
        end
        
        
        
        
        function val = spike_array(obj)
        %%cell array of spike times around each trigger
            val = cell(length(obj.trigger_times), 1);
            for ii = 1 : length(obj.trigger_times)
                these_lims = obj.limits + obj.trigger_times(ii);
                val{ii} = obj.cluster.fr.restrict_times(these_lims) - obj.trigger_times(ii);
            end
        end
        
        
        
        function val = spike_convolution_avg(obj)
        %%gets average of spike convolutions around triggers
            val = obj.spike_convolutions();
            val = mean(val, 2);
        end
        
        
        
        function val = spike_convolutions(obj)
        %%gets spike convolution triggered on obj.trigger_times
            val = nan(length(obj.common_t), length(obj.trigger_times));
            for ii = 1 : length(obj.trigger_times)
                this_t = obj.common_t + obj.trigger_times(ii);
                val(:, ii) = obj.cluster.fr.get_convolution(this_t);
            end
        end
    end
end
