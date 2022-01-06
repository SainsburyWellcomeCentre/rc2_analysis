classdef HeatmapData < handle
% HeatmapData Class for helping with creation of the heatmaps
%
%   HeatmapData Properties:
%       raster_data         - #clusters x 1 cell array of instances of class RasterData (one for each cluster)
%       heatmap             - the #clusters x #time points firing rate heatmap
%       baseline            - window specifying the time to take for the baseline (to subtract from each trace)
%       response            - window specifying the time to take for the resposne period
%       baseline_display    - window specifying the amount of time to show pre-event
%       response_display    - window specifying the amount of time to show post-event
%
%     Dependent
%       baseline_idx        - boolean vector masking the `baseline` period
%       response_idx        - boolean vector masking the `response` period
%       baseline_display_idx - boolean vector masking the `baseline_display` period
%       common_t            - timebase for the x-axis of the heatmap
%
%   HeatmapData Methods:
%       compute_heatmap         - computes the heatmap given the `raster_data`
%       heatmap_response_order  - use `baseline` and `response` periods to sort the heatmap
%       sort_heatmap            - perform the sorting of the heatmap
%       
%   See also: RasterData

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
        % HeatmapData
        %
        %   HeatmapData(RASTER_DATA_ARRAY) collection of methods for
        %   handling heatmaps. Takes a cell array with each entry an
        %   instance of RasterData class.
       
            obj.raster_data = raster_data;
            obj.compute_heatmap();
        end
        
        
        
        function val = get.baseline_display_idx(obj)
        %%mask on common_t picking out the baseline period - used for display
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
        %%compute_heatmap Computes the heatmap given the `raster_data`
        %
        %   compute_heatmap() 
        %
        %   See also: RasterData.spike_convolution_avg()
        
            n_clusters = length(obj.raster_data);
            
            obj.heatmap = nan(n_clusters, length(obj.common_t));
            for ii = 1 : n_clusters
                obj.heatmap(ii, :) = obj.raster_data{ii}.spike_convolution_avg();
            end
            
            % subtract baseline
            obj.heatmap = bsxfun(@minus, obj.heatmap, mean(obj.heatmap(:, obj.baseline_display_idx), 2));
        end
        
        
        
        function idx = heatmap_response_order(obj)
        %%heatmap_response_order Get the index after sorting the rows of
        %%the heatmap according to response magnitude
        %
        %   INDEX = heatmap_response_order()
        %   uses the baseline and response periods to compute the delta FR for
        %   each cluster around the event, then return the INDEX which
        %   orders the data from highest magnitude to lowest.
        %
        %   See also: sort_heatmap
            
            baseline_fr = mean(obj.heatmap(:, obj.baseline_idx), 2);
            response_fr = mean(obj.heatmap(:, obj.response_idx), 2);
            
            delta_fr = response_fr - baseline_fr;
            
            [~, idx] = sort(delta_fr, 'ascend');
        end
        
        
        
        function val = sort_heatmap(obj, idx)
        %%sort_heatmap Sorts the heatmap according to an index
        %
        %   SORTED_HEATMAP = sort_heatmap(INDEX)
        %   sorts the heatmap in `heatmap` property, by index in INDEX.
        %   INDEX should be a #clusters x 1 vector specifying the new
        %   ordering. This function allows an arbitrary ordering of the
        %   heatmap according to different criteria.
        %
        %   See also: heatmap_response_order
        
            val = obj.heatmap(idx, :);
        end
    end
end
