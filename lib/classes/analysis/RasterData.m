classdef RasterData < handle
% RasterData Class for handling firing rate data around a set of events
%
%   RasterData Properties:
%       cluster             - instance of class Cluster
%       limits              - 1x2 vector giving time limits, in seconds, around each event to look 
%                             for spikes and take the spike convolution
%                             (default, [-1, 1])
%       trigger_times       - #events x 1, vector of times (in probe time) around which 
%                             to get raster data
%       common_t            - #samples x 1 vector giving a common timebase
%                             for firing rate traces
%       fs                  - sample rate, Hz, to sample the firing rate
%                             convolution (default, 60Hz)
%
%   RasterData Methods:
%       spike_array         - return a cell array with spike times around the events
%       spike_convolutions  - return a matrix with spike FR convolutions around the events
%
%   
%   After creation of the object, the event times should be set by setting
%   the `trigger_times` property
%
%   See also: Cluster, FiringRate


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
        % RasterData
        %
        %   RasterData(CLUSTER)
        %   collection of methods for handling raster data. Takes as
        %   argument CLUSTER an instance of class Cluster.
        
            obj.cluster = cluster;
        end
        
        
        
        function val = get.common_t(obj)
        %%common timebase for the firing rate convolution
            val = linspace(obj.limits(1), obj.limits(2), round(diff(obj.limits)*obj.fs)+1);
        end
        
        
        
        
        function val = spike_array(obj)
        %%spike_array Return a cell array with spike times around the
        %%events
        %
        %   SPIKE_TIME_ARRAY = spike_array()
        %   for each event looks for spikes for the `cluster` around the
        %   event. Returns these in relative time (relative to the time of
        %   the event) as a cell array of vectors, SPIKE_TIME_ARRAY, of
        %   size #events x 1 with each entry containing the times of the
        %   spikes around the event. 
        
            val = cell(length(obj.trigger_times), 1);
            for ii = 1 : length(obj.trigger_times)
                these_lims = obj.limits + obj.trigger_times(ii);
                val{ii} = obj.cluster.fr.restrict_times(these_lims) - obj.trigger_times(ii);
            end
        end
        
        
        
        function val = spike_convolution_avg(obj)
        %%spike_convolution_avg Return average of firing rate convolutions
        %%around the events
        %
        %   AVG_FIRING_RATE = spike_convolution_avg()
        %
        %   Averages the firing rate convolutions returned by
        %   `spike_convolutions` across events, to give a event triggered
        %   average firing rate.
        %
        %   See also: spike_convolutions
        
            val = obj.spike_convolutions();
            val = mean(val, 2);
        end
        
        
        
        function val = spike_convolutions(obj)
        %%spike_convolutions Return a matrix with spike FR convolutions
        %%around the events
        %
        %   FIRING_RATE_MTX = spike_convolutions()
        %   for each event convolves the spikes for the `cluster` around the
        %   event. Returns these firing rate traces as a matrix of size
        %   #samples x #events, where #samples is determined by the
        %   `limits` and `fs` properties (see also `common_t`)
        %
        %   See also: FiringRate
        
            val = nan(length(obj.common_t), length(obj.trigger_times));
            for ii = 1 : length(obj.trigger_times)
                this_t = obj.common_t + obj.trigger_times(ii);
                val(:, ii) = obj.cluster.fr.get_convolution(this_t);
            end
        end
    end
end
