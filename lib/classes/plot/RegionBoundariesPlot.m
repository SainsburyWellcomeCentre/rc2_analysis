classdef RegionBoundariesPlot < handle
% RegionBoundariesPlot Class for plotting anatomy region boundaries from ProbeTrack information
%
%  RegionBoundariesPlot Properties:
%       probe_track     - object of type ProbeTrack
%       boundaries      - (#regions+1) x 1 vector of boundaries
%       region_str      - #regions x 1 cell array of strings with region names
%       middle_l5       - location of the middle of L5 from anatomy
%       
%  RegionBoundariesPlot Methods:
%       plot_boundary_lines     - plot the boundary lines on a figure
%       apply_offset            - apply an offset to the boundaries
%       print_region_labels     - print the labels to the figure
%       plot_middle_of_VISp5    - use a line to indicate middle of VISp5
%       highlight_fp            - highlight the region 'fp' (white matter)
%
%  See also: ProbeTrack, MIDepthPlot

    properties (SetAccess = private)
        
        probe_track
        boundaries
        region_str
        middle_l5
        
        h_boundary_lines = {}
        h_middle_l5
        h_region_str_txts = {}
        h_fp_fill
        
        current_offset = 0
    end
    
    
    
    methods
        
        function obj = RegionBoundariesPlot(probe_track)
        %%RegionBoundariesPlot
        %
        %   RegionBoundariesPlot(PROBE_TRACK) takes an object of type
        %   ProbeTrack and uses the information in there to plot region
        %   boundaries on a given axis.
            
            if isempty(probe_track); return; end
            
            obj.probe_track = probe_track;
            [obj.boundaries, ~, obj.region_str] = probe_track.region_boundaries();
            obj.middle_l5 = obj.probe_track.mid_l5_visp();
        end
        
        
        
        function plot_boundary_lines(obj, h_ax)
        %%plot_boundary_lines Plots the boundary lines on an axis
        %
        %   plot_boundary_lines(AXIS_HANDLE), plots the boundaries on the
        %   axis given by AXIS_HANDLE, which must exist.
        
            obj.h_boundary_lines = {};
            
            yl = get(h_ax, 'ylim');
            
            for ii = 1 : length(obj.boundaries)-1
                
                % boundary
                obj.h_boundary_lines{ii} = line(h_ax, obj.boundaries(ii)*[1, 1]+obj.current_offset, yl, ...
                    'color', 'k', ...
                    'linestyle', '--', ...
                    'linewidth', 0.25);
            end
        end
        
        
        
        function apply_offset(obj, offset)
        %%apply_offset Apply an offset to the boundaries in `h_boundary_lines`
        %
        %   apply_offset(OFFSET) shifts the boundary lines which have been
        %   created using `plot_boundary_lines`, so that they have an
        %   offset of OFFSET. Note: offset applied is absolute and not relative 
        %   to the current offset.
        
            obj.current_offset = offset;
            
            if ~isempty(obj.h_boundary_lines)
                for ii = 1 : length(obj.boundaries)-1
                    set(obj.h_boundary_lines{ii}, 'xdata', obj.boundaries(ii)*[1, 1] + obj.current_offset);
                end
            end
            
            if ~isempty(obj.h_region_str_txts)
                for ii = 1 : length(obj.boundaries)-1
                    current_pos = get(obj.h_region_str_txts{ii}, 'position');
                    set(obj.h_region_str_txts{ii}, ...
                        'position', [mean(obj.boundaries([ii, ii+1]))+ obj.current_offset, current_pos(2), 0]);
                end
            end
            
            if ~isempty(obj.h_middle_l5)
                set(obj.h_middle_l5, 'xdata', obj.middle_l5*[1, 1] + obj.current_offset);
            end
            
            if ~isempty(obj.h_fp_fill)
                set(obj.h_fp_fill, 'xdata', obj.boundaries([ii, ii+1, ii+1, ii]) + obj.current_offset);
            end
        end
        
        
        
        function print_region_labels(obj, h_ax)
        %%print_region_labels Prints the region labels on the axis
        %
        %   print_region_labels(AXIS_HANDLE), plots the region labels on the
        %   axis given by AXIS_HANDLE, which must exist.
        
            obj.h_region_str_txts = {};
            
            for ii = 1 : length(obj.boundaries)-1
                obj.h_region_str_txts{ii} = ...
                    text(h_ax, sum(obj.boundaries([ii, ii+1]))/2 + obj.current_offset, max(get(h_ax, 'ylim')), obj.region_str{ii}, ...
                    'verticalalignment', 'bottom', ...
                    'horizontalalignment', 'center');
            end
        end
        
        
        
        function plot_middle_of_VISp5(obj, h_ax)
        %%plot_middle_of_VISp5 Plot line to indicate middle of VISp5
        %
        %   plot_middle_of_VISp5(AXIS_HANDLE), plots the region labels on the
        %   axis given by AXIS_HANDLE, which must exist.
        
            obj.h_middle_l5 = line(h_ax, [1, 1]*obj.middle_l5 + obj.current_offset, get(h_ax, 'ylim'), 'color', 'r');
        end
        
        
        
        function highlight_fp(obj, h_ax)
        %%highlight_fp Highlight the region 'fp' (white matter)
        %
        %   plot_middle_of_VISp5(AXIS_HANDLE), plots the region labels on the
        %   axis given by AXIS_HANDLE, which must exist.
            
            idx = find(strcmp(obj.region_str, 'fp'));
            yl = get(h_ax, 'ylim');
            obj.h_fp_fill = patch(h_ax, obj.boundaries([idx, idx+1, idx+1, idx]) + obj.current_offset, ...
                yl([1, 1, 2, 2]), [207, 185, 151]/255, 'edgecolor', 'none', 'facealpha', 0.6);
        end
    end
end
