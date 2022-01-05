classdef UnityPlot < RC2Axis
% UnityPlot Class for helping with unity plots.
%
%   UnityPlot Properties:
%       x       - N x 1, x-axis value of data points
%       y       - N x 1, y-axis value of data points
%       histogram_on - true or false (default), whether to plot a diagonal histogram
%       h_patch - handle to MATLAB patch objects for the histogram
%       h_line  - handle to line along the diagonal
%
%   UnityPlot Methods:
%       min     - get minimum of both axes (x- and y-)
%       max     - get maximum of both axes (x- and y-)
%       xlim    - get or set x- and y-axis limits (calls xylim)
%       ylim    - get or set x- and y-axis limits (calls xylim)
%       xylim   - get or set x- and y-axis limits
%       add_histogram - include a diagonal histogram on the plot
%       remove_histogram - remove the diagonal istogram on the plot
%
%   UnityPlot is used for plotting paired values against each other when
%   one is interested in seeing if those values are similar. Therefore, the
%   limits of the x- and y-axis should be the same, and we have a line
%   along 'unity'.
%
%   UnityPlot does not have a `plot` function but rather expects to be
%   subclassed by a class with a `plot` function.
%
%   See also: UnityPlotPopulation, UnityPlotSingleCluster
%
%   TODO:   1. have `plot` function here which get's overwritten, otherwise
%              h_line, h_txt may not exist


    properties
        
        h_line
        h_txt
        
        x
        y
        
        min
        max
        
        histogram_on = false;
        h_patch = matlab.graphics.primitive.Patch.empty()
        bin_width
    end
    
    
    
    methods
        
        function obj = UnityPlot(x, y, h_ax)
        %%UnityPlot
        %
        %   UnityPlot(X, Y, AXIS_HANDLE) prepares an object of UnityPlot
        %   class containing X- and Y- data to be plotted on the same scale
        %   on each axis. AXIS_HANDLE is optional, if supplied it should be a
        %   handle to an axis object. Otherwise, an axis will be created.
        %
        %   X and Y are the x-axis and y-axis data respectively and should
        %   be vectors of the same length.
        
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
            
            obj.x = x;
            obj.y = y;
        end
        
        
        
        function val = get.min(obj)
        %%get minimum of both axes
            val = min([get(obj.h_ax, 'xlim'), get(obj.h_ax, 'ylim')]); %#ok<*CPROP>
        end
        
        
        
        function val = get.max(obj)
        %%get maximum of both axes
            val = max([get(obj.h_ax, 'xlim'), get(obj.h_ax, 'ylim')]);
        end
        
        
        
        function val = ylim(obj, val)
        %%ylim Get or set x and y limits to have the same range
        %
        %   VALUEOUT = ylim(VALUEIN) is a direct call to xylim.
        %   See also: xylim
        
            VariableDefault('val', []);
            val = obj.xylim(val);
        end
        
        
        
        function val = xlim(obj, val)
        %%xlim Get or set x and y limits to have the same range
        %
        %   VALUEOUT = xlim(VALUEIN) is a direct call to xylim.
        %   See also: xylim
        
            VariableDefault('val', []);
            val = obj.xylim(val);    
        end
        
        
        
        function val = xylim(obj, val)
        %%xylim Get or set x and y limits to have the same range
        %
        %   VALUEOUT = xylim(VALUEIN) returns VALUEOUT which is the x- and y-axis
        %   limits of the axes. If VALUEIN is not supplied or empty, the
        %   current x- and y-axis limits will be returned. If VALUEIN is supplied, it
        %   should be a 1x2 array [min, max] where max > min, and the
        %   x- and y-axis limits are updated to this value. In this case, VALUEOUT
        %   will be the same as VALUEIN.
        
            val = obj.get_set_limits(val, 'ylim');
            val = obj.get_set_limits(val, 'xlim');
            
            set(obj.h_line, 'xdata', val, 'ydata', val);
            set(obj.h_txt, 'position', [val(2), val(1), 0]);
            
            xtick = get(obj.h_ax, 'xtick');
            ytick = get(obj.h_ax, 'ytick');
            
            if length(xtick) > length(ytick)
                set(obj.h_ax, 'ytick', xtick);
            elseif length(xtick) < length(ytick)
                set(obj.h_ax, 'xtick', ytick);
            end
            
            if obj.histogram_on
                obj.remove_histogram();
                obj.add_histogram();
            end
        end
        
        
        
        function add_histogram(obj, bin_width)
        %%add_histogram Adds a diagonal histogram to the axis
        %
        %   add_histogram(BIN_WIDTH) adds a histogram running diagonally
        %   across the unity plot, to look at the distance of each point
        %   from the unity line. BIN_WIDTH determines the bin width of
        %   the histogram. If BIN_WIDTH is not supplied, the value in the
        %   `bin_width` property is used. If neither BIN_WIDTH or
        %   `bin_width` have a value, the MATLAB histcount function is
        %   called without an argument and bins are computed
        %   automatically.
        
            VariableDefault('bin_width', []);
            
            obj.remove_histogram();
            
            % distance of all points perpendicular to unity line
            points = [obj.x(:), obj.y(:)];
            
            dist_to_unity = nan(size(points, 1), 1);
            
            for p_idx = 1 : size(points, 1)
                
                % if both x and y are zero distance is zero
                %   otherwise calculation would return nan
                if points(p_idx, 1) == 0 && points(p_idx, 2) == 0
                    dist_to_unity(p_idx) = 0;
                    continue
                end
                
                r = norm(points(p_idx, :));
                theta = acos(sum(points(p_idx, :))/(r*sqrt(2)));
                
                if points(p_idx, 1) > points(p_idx, 2)
                    dist_to_unity(p_idx) = r * sin(theta);
                else
                    dist_to_unity(p_idx) = - r * sin(theta);
                end
                
            end
            
            if isempty(bin_width) && isempty(obj.bin_width)
                [bin_count, bin_edges] = histcounts(dist_to_unity);
            elseif ~isempty(bin_width)
                obj.bin_width = bin_width;
                [bin_count, bin_edges] = histcounts(dist_to_unity, 'binwidth', bin_width);
            else
                [bin_count, bin_edges] = histcounts(dist_to_unity, 'binwidth', obj.bin_width);
            end
            
            M = obj.max();
            m = obj.min();
            
            % h_limit is the middle of the axis to a corner of the axis
            h_limit = (M - m) / sqrt(2);
            
            % scale the bin counts by the h_limit
            bin_height = bin_count * (h_limit / max(bin_count)) / 2; %#ok<CPROPLC>
            
            % for each bin of the histogram
            for bin_i = 1 : length(bin_height)
                
                % get edges of the bar
                p1 = [(M+m)/2, (M+m)/2] + [1, -1]*bin_edges(bin_i)/sqrt(2);
                p2 = [(M+m)/2, (M+m)/2] + [1, -1]*bin_edges(bin_i+1)/sqrt(2);
                p3 = p2 + bin_height([bin_i, bin_i])/sqrt(2);
                p4 = p1 + bin_height([bin_i, bin_i])/sqrt(2);
                
                % plot patch
                obj.h_patch(bin_i) = patch(obj.h_ax, [p1(1), p2(1), p3(1), p4(1)], ...
                    [p1(2), p2(2), p3(2), p4(2)], 'k');
                
                set(obj.h_patch(bin_i), 'facecolor', 'none', 'edgecolor', 'k');
                
            end
            
            obj.histogram_on = true;
        end
        
        
        
        function remove_histogram(obj)
        %%remove_histogram Removes the diagonal histogram from the axis
        %
        %   remove_histogram() removes the histogram set with
        %   `add_histogram` from the axis.
        
            if ~isempty(obj.h_patch)
                
                for p = 1 : length(obj.h_patch)
                    delete(obj.h_patch(p));
                end
                obj.h_patch = matlab.graphics.primitive.Patch.empty();
                
            end
            
            obj.histogram_on = false;
        end
    end
end
