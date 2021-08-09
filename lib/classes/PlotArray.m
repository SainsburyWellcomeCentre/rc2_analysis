classdef PlotArray < handle
    
    properties
        
        nx
        ny
        paper_type
        
        ax_size_cm = 2.3
    end
    
    properties
        
        nx_total = 6
        ny_total = 6
    end
    
    methods
        
        function obj = PlotArray(ny, nx, paper_type)
            
            VariableDefault('paper_type', 'A4');
            
            obj.nx = nx;
            obj.ny = ny;
            obj.paper_type = paper_type;
            
        end
        
        
        function pos = get_position(obj, i)
            
            assert(i <= obj.nx*obj.ny, 'Index is too high');
            
            pos = obj.grid();
            row = ceil(i / obj.nx);
            col = mod(i-1, obj.nx) + 1;
            
            idx = sub2ind([obj.nx_total, obj.ny_total], col, row);
            
            pos = pos{idx};
        end
        
        
        function pos = grid(obj)
            
            if strcmp(obj.paper_type, 'A3')
                x_size = 29.7;
                y_size = 42;
            elseif strcmp(obj.paper_type, 'A4')
                x_size = 21;
                y_size = 29.7;
            end
            
            x_gap = (x_size - (obj.nx_total * obj.ax_size_cm)) / (obj.nx_total + 1);
            y_gap = 2*x_gap;
            y_top_gap = (y_size - x_size)/2;
            
            pos = {};
            for y = 1 : obj.ny_total
                for x = 1 : obj.nx_total
                    pos{end+1} = [x*x_gap + (x-1)*obj.ax_size_cm, ...
                                  y_size - y_top_gap - y*obj.ax_size_cm - (y-1)*y_gap, ...
                                  obj.ax_size_cm, obj.ax_size_cm];
                end
            end
        end
    end
end