classdef MUACluster < handle
% MUACluster Class for handling data from multiple clusters and treat them
% as a single "MUA cluster"
%
%   MUACluster Properties:
%       clusters        - an array of `Cluster`  objects
%       spike_times     - #spikes x 1 vector of spike times for all the clusters
%       fr              - instance of class FiringRate handling spike rates of all clusters

%   MUACluster Methods:
%
%   See also: Cluster, FiringRate

    properties (SetAccess = private)
        
        clusters
        spike_times
        fr
    end
    
    
    
    methods
        
        function obj = MUACluster(clusters)
            
            obj.clusters = clusters;
            
            % concatenate all spike times from all clusters
            obj.spike_times = sort(cat(1, clusters(:).spike_times));
            
            obj.fr = FiringRate(obj.spike_times);
        end
    end
end