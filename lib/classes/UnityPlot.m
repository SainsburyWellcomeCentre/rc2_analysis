classdef UnityPlot < RC2Axis
    
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
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
            
            obj.x = x;
            obj.y = y;
            
        end
        
        
        
        function val = get.min(obj)
            
            val = min(get(obj.h_ax, 'xlim')); %#ok<*CPROP>
        end
        
        
        function val = get.max(obj)
            
            val = max(get(obj.h_ax, 'xlim'));
        end
        
        
        function val = ylim(obj, val)
            
            VariableDefault('val', []);
            obj.xylim(val);
            
        end
        
        
        function val = xlim(obj, val)
            
            VariableDefault('val', []);
            obj.xylim(val);
            
        end
        
        
        function val = xylim(obj, val)
            
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