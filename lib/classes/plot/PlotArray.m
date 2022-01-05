classdef PlotArray < handle
% PlotArray Class for handling the arrangement of axes in an A4 portrait figure
%
%   PlotArray Properties:
%       nx              - number of axes along a row (width), must be < `nx_total`
%       ny              - number of axes along a column (height), must be < `ny_total`
%       paper_type      - 'A4' 
%       ax_size_cm      - size of each axis in the grid
%       nx_total        - number of axes in the grid
%       ny_total        -
%   
%   PlotArray Methods:
%       get_position    - get the position in the array of an axis
%       grid            - underlying grid of axes positions
%
%   This class gets the positions of an axis as if it were part of a grid
%   of axes.  `nx_total` and `ny_total` determine the number of axes in the
%   underlying grid. `nx` and `ny` is the number of axes in the grid I want to work with.
%
%   It's done like this so that if I want a 6 x 6 grid of axes for one case, and a 3 x 2
%   grid for another, yet I want the axes of the 3 x 2 grid to appear in
%   the same location as the top left hand axes of the 6 x 6. 
%   However, it could probably have been done in a clearer way.
%
%   Currently this only works for A4 portrait figure.
%
%   TODO:   1. checks on the value of nx, ny, ax_size_cm
%           2. change name of class

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
        %%PlotArray
        %
        %   PlotArray(N_AXES_HEIGHT, N_AXES_WIDTH) prepares the class to
        %   work with a grid of axes with N_AXES_HEIGHT axes in height and
        %   N_AXES_WIDTH axes in width, on top of an underlying axis of
        %   size `ny_total` height and `nx_total` width.
        
            VariableDefault('paper_type', 'A4');
            
            obj.nx = nx;
            obj.ny = ny;
            obj.paper_type = paper_type;
            
        end
        
        
        function pos = get_position(obj, i)
        %%get_position Get the position in the array of an axis
        %
        %   get_position(INDEX) gets the position in the axis array of the
        %   axis with index INDEX. INDEX is a linear axis starting at the
        %   top-left axis, moving along the top row **until it reaches `nx`** 
        %   then wraps around to the left-most axis of the next row etc. until it gets to 
        %   the axis at (`nx`, `ny`) (do not confuse with index of
        %   the underlying grid).
        
            assert(i <= obj.nx*obj.ny, 'Index is too high');
            
            pos = obj.grid();
            row = ceil(i / obj.nx);
            col = mod(i-1, obj.nx) + 1;
            
            idx = sub2ind([obj.nx_total, obj.ny_total], col, row);
            
            pos = pos{idx};
        end
        
        
        function pos = grid(obj)
        %%grid Underlying grid of axes positions
        %
        %   POSITIONS = grid() computes the positions of a grid of axes on
        %   an A4 portrait figure each with size `ax_size_cm`. The grid
        %   will contain `nx_total` axes across the page width and `ny_total`
        %   axes up the page.
        %   POSITIONS is a cell array with the position of each axis. Each
        %   entry is of the form [from left, from right, width height] with
        %   the position of the axis in 'cm' on the page.
        %
        %   POSITIONS is a linear cell array, and the axis index goes from
        %   top-left, along the top row, then wraps to the next row etc.
        %   until it gets to bottom-right.
        
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