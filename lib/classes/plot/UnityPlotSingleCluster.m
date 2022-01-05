classdef UnityPlotSingleCluster < UnityPlot
% UnityPlotSingleCluster Class for plotting unity plots containing trial-by-trial data for a single cluster
% (c.f. UnityPlotPopulation)
%
%   UnityPlotSingleCluster Properties:
%       p            - # clusters x 1, vector of p-values containing the result of the
%                      statistical test on response of each cluster
%       change_direction  - # clusters x 1, vector with either 1, 0 or -1
%                      indicating whether response was +ve, no change or -ve respectively
%       marker_style - valid MATLAB marker style (default = 'o')
%       h_dot        - array of handles to the dots in the unity plot
%       dot_size     - size of median dot (default = 10 points)
%
%   UnityPlotSingleCluster Methods:
%       plot         - plot the unity plot
%
%     Internal:
%       sync_axes       - make sure axes are synced
%       add_unity_line  - adds the unity line
%       add_median_dot  - adds the median dot
%       get_p_colour    - gets the colour of the median dot
%       print_stats     - print statistics on the axes
%
%   UnityPlotSingleCluster is a subclass of UnityPlot containing information
%   about the trial to trial responses of a cluster. Each individual point
%   in `x` and `y` is coloured grey, and the median gets a colour depending
%   on the `p` and `change_direction` properties.
%
%   See also: UnityPlot, UnityPlotPopulation

    properties
        
        p
        change_direction
        marker_style = 'o'
        h_dot    
        dot_size = 10
    end
    
    
    
    methods
        
        function obj = UnityPlotSingleCluster(x, y, p, change_direction, h_ax)
        %%UnityPlotSingleCluster
        %
        %   UnityPlotSingleCluster(X, Y, P_VALUE, DIRECTION, AXIS_HANDLE)
        %   prepares an object of UnityPlotSingleCluster
        %   class containing X- and Y- data to be plotted on the same scale
        %   on each axis. AXIS_HANDLE is optional, if supplied it should be a
        %   handle to an axis object. Otherwise, an axis will be created.
        %
        %   X and Y are the x-axis and y-axis data respectively. They are,
        %   for instance, the firing rate of a single cluster to two
        %   separate conditions (which can be paired in some way) for each
        %   trial. The median value of the X and Y data is also plotted.
        %
        %   P_VALUE is the result of a statistical test for a each cluster,
        %   e.g. Wilcoxon sign rank test of the x-values vs. the y-values.
        % 
        %   DIRECTION indicates whether the response increased (1), decrease (-1)
        %   or stayed the same (note we have this extra information because in a subset of
        %   cases, median values can be 0 for X and Y, but the result of a Wilcoxon
        %   sign rank test is significant).
        %
        %   X, Y should be vectors of the same length. Unlike UnityPlotPopulation P_VALUE and
        %   DIRECTION are single values.
        
            VariableDefault('h_ax', []);
            
            len = [length(x), length(y)];
            assert(length(unique(len)) == 1, 'inputs do not have equal lengths');
            
            obj = obj@UnityPlot(x, y, h_ax);
            
            obj.p = p;
            obj.change_direction = change_direction;
        end
        
        
        function plot(obj)
        %%plot Creates the unity plot
        %
        %   plot() plots the data in `y` vs. data in `x`.  All points are
        %   grey. The median is also plotted. If `p` is >= 0.05 the median dot is black. 
        %   If `p` is < 0.05 the median dot is red if `direction` is 1 and blue if `direction` is -1.
        
        
            n_trials = length(obj.x);
            
            if sum(isnan(obj.x)) > 0 || sum(isnan(obj.y)) > 0
                warning('There are %i x and %i y NaN FR values', sum(isnan(obj.x)), sum(isnan(obj.x)));
            end
            
            obj.h_dot = cell(n_trials, 1);
            
            scatter(obj.h_ax, obj.x, obj.y, obj.dot_size, [0.5, 0.5, 0.5], 'fill');
            
            obj.add_median_dot();
            obj.sync_axes();
            obj.add_unity_line();
            obj.print_stats();
        end
        
        
        function sync_axes(obj)
        %%sync_axes Ensure x and y axis have same range
        
            m = obj.min;
            M = obj.max;
            obj.xlim([m, M]);
            obj.ylim([m, M]);
        end
        
        
        
        function add_unity_line(obj)
        %%add_unity_line Add unity line to plot
        
            axis_limits = obj.xlim();
            if ~isempty(obj.h_line)
                delete(obj.h_line);
            end
            obj.h_line = line(obj.h_ax, axis_limits, axis_limits, 'linestyle', '--', 'color', 'k');
        end
        
        
        
        function add_median_dot(obj)
        %%add_median_dot Adds the median dot to the figure
        
            x_med = nanmedian(obj.x);
            y_med = nanmedian(obj.y);
            x_iqr = prctile(obj.x, [25, 75]);
            y_iqr = prctile(obj.y, [25, 75]);
            col = obj.get_p_colour();
            
            scatter(obj.h_ax, x_med, y_med, obj.dot_size, col, 'fill')
            line(obj.h_ax, x_iqr, y_med([1, 1]), 'color', col);
            line(obj.h_ax, x_med([1, 1]), y_iqr, 'color', col);
        end
        
        
        
        function col = get_p_colour(obj)
        %%get_p_colour Get the colour of the median dot from p-value
        %
        %   COLOUR = get_p_colour() uses the `p` value and `direction` to
        %   compute the colour of the median dot.
        
            if obj.p < 0.05 && obj.change_direction == -1
                col = [30,144,255]/255;
            elseif obj.p < 0.05 && obj.change_direction == 1
                col = 'r';
            elseif obj.p < 0.05 && obj.change_direction == 0
                col = 'm';
            else
                col = 'k';
            end
        end
        
        
        function print_stats(obj)
        %%print_stats Print the statistics
       
            m = obj.min;
            M = obj.max;
            
            if obj.p < 0.01
                format_str = 'p = %.2e \nx median = %.2f (n=%i) \ny median = %.2f (n=%i)\n# nan (x=%i,y=%i)';
            else
                format_str = 'p = %.2f \nx median = %.2f (n=%i) \ny median = %.2f (n=%i)\n# nan (x=%i,y=%i)';
            end
            
            txt_str = sprintf(format_str, ...
                    obj.p, ...
                    nanmedian(obj.x), ...
                    sum(~isnan(obj.x)), ...
                    nanmedian(obj.y), ...
                    sum(~isnan(obj.y)), ...
                    sum(isnan(obj.x)), ...
                    sum(isnan(obj.y)));
            
            obj.h_txt = text(obj.h_ax, M, m, txt_str, 'verticalalignment', 'bottom', 'horizontalalignment', 'right', ...
                'fontsize', 4);
        end
    end
end
