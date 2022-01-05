classdef RC2Axis < handle
% RC2Axis Class for general handling of axes
%
%   Many of the plotting classes are subclasses of this class.
%
%  RC2Axis Properties:
%       h_fig           - handle to the figure containing the axes
%       h_ax            - handle to the axes
%       h_title         - handle to the title text
%       blue            - code for plotting blue objects
%       red             - code for plotting red objects
%       black           - code for plotting black objects
%   
%  RC2Axis Methods:
%       xlabel          - set the x-axis label
%       ylabel          - set the y-axis label
%       xlim            - get or set the x-axis limits
%       ylim            - get or set the y-axis limits
%       ymin            - return min on the y axis
%       xmin            - return min on the x axis
%       ymax            - return max on the y axis
%       xmax            - return max on the x axis
%       get_set_limits  - get or set limits of the axes
%       title           - set title of axis
%
%   TODO:   1. The xlabel and ylabel methods may not work

    properties (SetAccess = protected)
        
        h_fig
        h_ax
        h_title
    end
    
    properties
        
        blue = [30,144,255]/255
        red = 'r'
        black = 'k'
    end
    
    
    
    methods
        
        function obj = RC2Axis(h_ax)
        %%RC2Axis
        %
        %   RC2Axis(AXIS_HANDLE) associates the object with the axis
        %   referenced by AXIS_HANDLE, if supplied. Otherwise, if empty, an
        %   axis is generated.
        
            if isempty(h_ax)
                obj.h_fig = figure();
                obj.h_ax = axes;
            else
                obj.h_fig = h_ax.Parent;
                obj.h_ax = h_ax;
            end
            
            % hold the axes
            hold(obj.h_ax, 'on');
            set(obj.h_ax, 'fontsize', 8);
        end
        
        
        function xlabel(obj, str)
        %%xlabel Sets x-axis label.
        %
        %   xlabel(STRING) sets x-axis label to string in STRING.
        
%             txt = text('string', str, 'fontsize', 8);
%             set(obj.h_ax, 'xlabel', txt, 'interpreter', 'none');
            xlabel(str, 'interpreter', 'none', 'fontsize', 8);
        end
        
        
        function ylabel(obj, str)
        %%ylabel Sets y-axis label.
        %
        %   ylabel(STRING) sets y-axis label to string in STRING.
        
%             txt = text('string', str, 'fontsize', 8);
%             set(obj.h_ax, 'ylabel', txt, 'interpreter', 'none');
            ylabel(str, 'interpreter', 'none', 'fontsize', 8);
        end
        
        
        function val = ylim(obj, val)
        %%ylim Get or set the y-axis limits
        %
        %   VALUEOUT = ylim(VALUEIN) returns VALUEOUT which is the y-axis
        %   limits of the axis. If VALUEIN is not supplied or empty, the
        %   y-axis limits will be returned. If VALUEIN is supplied, it
        %   should be a 1x2 array [ymin, ymax] where ymax > ymin, and the
        %   y-axis limit is updated to this value. In this case, VALUEOUT
        %   will be the same as VALUEIN.
        %
        %   See also: ylim
        
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'ylim');
        end
        
        
        function val = xlim(obj, val)
        %%xlim Get or set the x-axis limits
        %
        %   VALUEOUT = xlim(VALUEIN) returns VALUEOUT which is the x-axis
        %   limits of the axis. If VALUEIN is not supplied or empty, the
        %   x-axis limits will be returned. If VALUEIN is supplied, it
        %   should be a 1x2 array [xmin, xmax] where xmax > xmin, and the
        %   x-axis limit is updated to this value. In this case, VALUEOUT
        %   will be the same as VALUEIN.
        %
        %   See also: xlim
        
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'xlim');
        end
        
        
        function val = ymin(obj)
        %return min on the y axis
            val = min(get(obj.h_ax, 'ylim'));
        end
        
        
        
        function val = xmin(obj)
        %return min on the x axis
            val = min(get(obj.h_ax, 'xlim'));
        end
        
        
        
        function val = ymax(obj)
        %return max on the y axis
            val = max(get(obj.h_ax, 'ylim'));
        end
        
        
        
        function val = xmax(obj)
        %return max on the x axis
            val = max(get(obj.h_ax, 'xlim'));
        end
        
        
        
        function val = get_set_limits(obj, val, str)
        %%get_set_limits Get or set limits of the axes
        %
        %   get_set_limits(VALUE, STRING)
        %   Args:
        %       VALUE - 1 x 2 vector of axis limits (if empty then the current
        %               value is returned)
        %       STRING - 'xlim' or 'ylim' determining axis to get or set
        
            if isempty(val)
                val = get(obj.h_ax, str);
                return
            end
            
            % replace nan values with existing limits
            I = isnan(val);
            if any(I)
                l = get(obj.h_ax, str);
                val(I) = l(I);
            end
                
            set(obj.h_ax, str, val);
        end
        
        
        
        function title(obj, str)
        %%title Set title of axis
        %
        %   title(STRING) sets title to STRING.
        
            obj.h_title = title(obj.h_ax, str, 'interpreter', 'none', 'fontsize', 10);
        end
        
    end
end