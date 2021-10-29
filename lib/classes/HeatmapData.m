classdef HeatmapData < handle
    
    properties (SetAccess = private)
        
        raster_data
        heatmap
    end
    
    properties
        
        baseline = [-0.4, 0]
        response = [0, 0.4]
        
        baseline_display = [-1, 0]
        response_display = [0, 1]
    end
    
    properties (Dependent = true)
        
        baseline_idx
        response_idx
        baseline_display_idx
        common_t
    end
    

    
    methods
        
        function obj = HeatmapData(raster_data)
        %%collection of methods for handling heatmaps
            obj.raster_data = raster_data;
            obj.compute_heatmap();
        end
        
        
        
        function val = get.baseline_display_idx(obj)
        %%mask on common_t picking out the baseline period - used for
        %%display
            val = obj.common_t >= obj.baseline_display(1) & obj.common_t < obj.baseline_display(2);
        end
        
        
        
        function val = get.baseline_idx(obj)
        %%mask on common_t picking out the baseline period - used for
        %%repsonse calculation
            val = obj.common_t >= obj.baseline(1) & obj.common_t < obj.baseline(2);
        end
        
        
        
        function val = get.response_idx(obj)
        %%mask on common_t picking out the response period - used for
        %%repsonse calculation
            val = obj.common_t >= obj.response(1) & obj.common_t < obj.response(2);
        end
        
        
        
        function val = get.common_t(obj)
        %%common timebase for all traces
            val = obj.raster_data{1}.common_t;
        end
        
        
        
        function compute_heatmap(obj)
        %%calculate the heatmap index
            n_clusters = length(obj.raster_data);
            
            obj.heatmap = nan(n_clusters, length(obj.common_t));
            for ii = 1 : n_clusters
                obj.heatmap(ii, :) = obj.raster_data{ii}.spike_convolution_avg();
            end
            
            % subtract baseline
            obj.heatmap = bsxfun(@minus, obj.heatmap, mean(obj.heatmap(:, obj.baseline_display_idx), 2));
        end
        
        
        
        function idx = heatmap_response_order(obj)
        %%use the baseline and response periods to compute the delta FR for
        %%each cluster around the trigger, then return the order of that
            
            baseline_fr = mean(obj.heatmap(:, obj.baseline_idx), 2);
            response_fr = mean(obj.heatmap(:, obj.response_idx), 2);
            
            delta_fr = response_fr - baseline_fr;
            
            [~, idx] = sort(delta_fr, 'ascend');
        end
        
        
        
        function val = sort_heatmap(obj, idx)
        %%sort heatmap according to index
            val = obj.heatmap(idx, :);
        end
    end
end