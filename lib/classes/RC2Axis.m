classdef RC2Axis < handle
    
    properties
        
        h_fig
        h_ax
        h_title
    
    end
    
    properties (Constant = true)
        
        blue = [30,144,255]/255
        red = 'r'
        black = 'k'
        
    end
    
    methods
        
        function obj = RC2Axis(h_ax)
            
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
            txt = text('string', str, 'fontsize', 8);
            set(obj.h_ax, 'xlabel', txt);
        end
        
        
        function ylabel(obj, str)
            txt = text('string', str, 'fontsize', 8);
            set(obj.h_ax, 'ylabel', txt);
        end
        
        
        function val = ylim(obj, val)
            
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'ylim');
        end
        
        
        function val = xlim(obj, val)
            
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'xlim');
        end
        
        
        
        function val = get_set_limits(obj, val, str)
            
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
            obj.h_title = title(obj.h_ax, str, 'interpreter', 'none', 'fontsize', 10);
        end
        
    end
end