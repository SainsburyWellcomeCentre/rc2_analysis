classdef Cluster < handle
% Cluster Class for handling cluster data such as spike times
%
%   Cluster Properties:
%       cluster         - a `cluster` structure from the formatted data file
%       id              - an integer specifying the ID of the cluster
%       spike_times     - #spikes x 1 vector of spike times in probe time
%       region_str      - string, acronym of brain region to which cluster assigned
%       depth           - depth of the cluster below the pia in um
%       spiking_class   - string, whether the cluster is 'wide' or 'narrow' spiking
%       distance_from_probe_tip - distance of the cluster from the probe tip in um
%       duration        - peak-to-trough duration of the average waveform in ms
%       overall_firing_rate - average firing rate of cluster across the probe recording
%       fr              - instance of class FiringRate handling spike rates of cluster across recording
%       shank_id        - shank ID on which cluster appears
%
%   Cluster Methods:
%       is_VISp         - boolean, is cluster in VISp
%
%   See also: FiringRate

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
        distance_from_probe_tip
        duration
        overall_firing_rate
        shank_id
    end
    
    
    
    methods
        
        function obj = Cluster(cluster)
        % Cluster
        %
        %   Cluster(CLUSTER) class for interacting a little with a cluster
        %   structure in CLUSTER, from the formatted data.
        
            obj.cluster = cluster;
            obj.fr = FiringRate(cluster.spike_times);
        end
        
        
        
        function val = get.id(obj)
        %%cluster id
            val = obj.cluster.id;
        end
        
        
        
        function val = get.spike_times(obj)
        %%return cluster spike times (probe timebase)
            val = obj.cluster.spike_times;
        end
        
        
        
        function val = get.region_str(obj)
        %%which region is the cluster in
            val = obj.cluster.region_str;
        end
        
        
        
        function val = get.depth(obj)
        %%depth of the cluster from pia in um
            val = obj.cluster.depth;
        end
        
        
        
        function val = get.distance_from_probe_tip(obj)
        %%distance of cluster from probe tip in um
            val = obj.cluster.distance_from_probe_tip;
        end
           
        
        
        function val = get.duration(obj)
        %%peak-to-trough duration of the average waveform in ms
            val = obj.cluster.duration;
        end
        
        
        
        function val = get.overall_firing_rate(obj)
        %%average firing rate of cluster across the probe recording
            val = obj.cluster.firing_rate;
        end
        
        
        
        function val = get.spiking_class(obj)
        %%string, whether the cluster is 'wide' or 'narrow' spiking
            if obj.cluster.duration < constants('spiking_class_threshold_ms')
                val = 'narrow';
            else
                val = 'wide';
            end
        end
        
        
        
        function val = get.shank_id(obj)
        %%shank ID on which cluster appears
            val = obj.cluster.shank_id;
        end
        
        
        
        function val = is_VISp(obj)
        %%is_VISp Return whether cluster is in VISp
        %
        %   IN_VISP = is_VISp()
        
            val = ~isempty(regexp(obj.region_str, 'VISp\d', 'once'));
        end
    end
end
