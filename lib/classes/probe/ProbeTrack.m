classdef ProbeTrack < handle
    
    properties
        
        offset = 0  % offset compared to electrophysiology
    end
    
    properties (Hidden = true, SetAccess = private)
        
        tbl
        probe
    end
    
    properties (Dependent = true)
        
        from_probe_tip
        region_id
        region_str
    end
    
    
    
    methods
        
        function obj = ProbeTrack(tbl, probe_type)
        %%for handling information contained in the probe track file
        %   tbl is the MATLAB table acquired by loading the track file
        %   (track_0.csv, track_1.csv etc) with 'readtable'
        %
        %   Can apply an offset on the boundaries
        %       +ve means boundaries shift up relative to probe
        
        
            VariableDefault('probe_type', '3A');
            
            obj.tbl = tbl;
            obj.probe = Probe(probe_type);
        end
        
        
        
        function val = get.from_probe_tip(obj)
        %%for each point in the probe track table get the distance from the
        %   probe tip, not accounting for the electrophysiology offset
            val  = obj.tbl.Position(end) - obj.tbl.Position;
        end
        
        
        
        function val = get.region_id(obj)
        %%region id at each mm of the probe track
        %   converts any NaN values into a numeric value -10
            val = obj.tbl.RegionID;
            val(isnan(val)) = -10;
        end
        
        
        
        function val = get.region_str(obj)
        %%region str at each mm of the probe track
            val = obj.tbl.RegionAcronym;
        end
        
        
        
        function [boundaries, region_id, region_str] = region_boundaries(obj)
        %%return where on the probe the region boundaries occur
        % as well as the region strings for those boundaries
        % without any correction for the electrophysiology
        %   Outputs:
        %       boundaries - 1x(N+1) vector containing the distance of a region
        %                   bounday from the probe tip in um, where N is
        %                   the number of regions in the probe track file
        %       region_id - 1xN vector containing the (Allen) region ID
        %                   jth entry is region ID between boundary j and
        %                   j+1
        %       region_str - 1xN cell array containing the (Allen) region acronym
        %                   jth entry is region ID between boundary j and
        %                   j+1
            
            % find where the region IDs change
            idx = find(abs(diff(obj.region_id)) > 0);
            
            boundaries = obj.from_probe_tip([1, idx(:)'+1, end]);
            region_id = obj.region_id([1, idx(:)'+1]);
            region_str = obj.region_str([1, idx(:)'+1]);
        end
        
        
        
        function [boundaries, region_id, region_str] = region_boundaries_adjusted(obj)
        %%return region boundaries relative to probe after correction for
        %   the difference between the position of electrophysiology and
        %   anatomy mid-L5
        
            [boundaries, region_id, region_str] = obj.region_boundaries();
            boundaries = boundaries + obj.offset;
        end
        
        
        
        function [region_id, region_str, from_pia] = get_region_of_point_from_tip(obj, from_tip)
        %%for a point 'from_tip' microns from the probe tip, get the region ID, region str and
        %   distance from pia of that point
        
            [region_id, region_str, from_pia] = obj.region_from_tip_shared(from_tip, 0);
        end
        
        
        
        function [region_id, region_str, from_pia] = get_region_of_point_from_tip_adjusted(obj, from_tip)
        %%for a point 'from_tip' microns from the probe tip, get the region ID, region str and
        %   distance from pia of that point
            
            [region_id, region_str, from_pia] = obj.region_from_tip_shared(from_tip, obj.offset);
        end
        
        
        
        function [region_id, region_str, from_pia] = region_from_tip_shared(obj, from_tip, offset)
        %%shared function
            
            pia_from_tip = obj.get_pia_from_tip + offset;
            from_pia = pia_from_tip - from_tip;
            
            if from_pia < 0
                region_id = -10;
                region_str = 'Not found in brain';
            elseif from_tip < offset
                region_id = -11;
                region_str = 'Beyond tracked region';
            else
                % get closest point
                [~, idx] = min(abs(from_tip - (obj.from_probe_tip + offset)));
                % return region ID and region acronym
                region_id = obj.region_id(idx);
                region_str = obj.region_str{idx};
            end
        end
        
        
        
        function val = mid_l5_visp(obj)
        %%get distance of mid layer 5 (VISp5) from the probe tip before any
        %   correction for electrophysiology L5 peak
            
            [boundaries, ~, region_str] = obj.region_boundaries(); %#ok<*PROP>
            
            idx = find(strcmp(region_str, 'VISp5'));
            
            if length(idx) == 1 && ~ismember(idx, [1, length(region_str)])
                val = mean(boundaries([idx, idx+1]));
            else
                % one mouse in which boundaries run across VISpm and VISp
                val = nan;
            end
        end
        
        
        
        function val = get_pia_from_tip(obj)
        %%return the distance of the pial surface from the probe tip
        %   ASSUMES THAT THE TOP MOST REGION IN THE TRACK FILE IS OF THE
        %   FORM VISp*1
            [boundaries, ~, regions] = obj.region_boundaries();
            visp1_cmp =  regexp(regions, 'VISp*\w1');
            visp1_idx = cellfun(@(x)(~isempty(x)), visp1_cmp);
            val = boundaries(visp1_idx);
        end
        
        
        
        function val = get_pia_from_tip_adjusted(obj)
        %%return the distance of the pial surface from the probe tip
        %   ASSUMES THAT THE TOP MOST REGION IN THE TRACK FILE IS OF THE
        %   FORM VISp*1
            val = obj.get_pia_from_tip();
            val = val + obj.offset;
        end
        
        
        
        function val = get_region_relative_depth(obj, from_tip)
            
            region_id = obj.get_region_of_point_from_tip(from_tip); %#ok<*PROPLC>
            [boundaries, region_ids] = obj.region_boundaries();
            
            val = obj.relative_depth_shared(from_tip, boundaries, region_ids, region_id);
        end
        
        
        
        function val = get_region_relative_depth_adjusted(obj, from_tip)
          
            region_id = obj.get_region_of_point_from_tip_adjusted(from_tip); %#ok<*PROPLC>
            [boundaries, region_ids] = obj.region_boundaries_adjusted();
            
            val = obj.relative_depth_shared(from_tip, boundaries, region_ids, region_id);
        end
        
        
        
        function val = relative_depth_shared(obj, from_tip, boundaries, region_ids, region_id)
            
            idx = find(region_ids == region_id);
            if length(idx) ~= 1
                val = nan; return
            end
            
            upper = boundaries(idx);
            lower = boundaries(idx+1);
            
            val = (upper - from_tip) / (upper - lower);
        end
    end
end
