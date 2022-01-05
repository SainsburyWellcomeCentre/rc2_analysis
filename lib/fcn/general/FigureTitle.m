classdef FigureTitle < handle
% FigureTitle Class for giving a figure a title not associated with any
% axis.
%
%   FigureTitle Properties:
%       string          - the title string
%       y_position      - the normalized position of the string
%       h_fig           - the parent figure handle
%       h_ax            - the hidden axis the text is associated with
%       h_text          - the handle to the title string
%
%   FigureTitle Methods:
%       set_string      - set title string
%       set_position    - set position of title in height
%
%   
%   TODO:   `set_string` doesn't set `string` property, neither does
%           `set_position` set `y_position`...

    properties
        string
        y_position
    end
    
    properties (Hidden = true)
        h_fig
        h_ax
        h_text
    end
    
    
    
    methods
        
        function obj = FigureTitle(h_fig, str)
        %%FigureTitle
        %
        %   FigureTitle(FIGURE_HANDLE, STRING) associates a text object with the figure with handle
        %   FIGURE_HANDLE, by creating a new axis covering the entire
        %   figure and putting a text object at the top. If STRING is
        %   supplied, it will be put at the top of the axis. If not, an
        %   empty string will be put there which can be filled later.
        
            VariableDefault('str', '');
            
            obj.h_fig = h_fig;
            figure(obj.h_fig);
            obj.h_ax = axes('position', [0, 0, 1, 1]);
            
            axis off;
            set(obj.h_ax, 'xlim', [0, 1], 'ylim', [0, 1]);
            obj.h_text = text(0.5, 0.97, str, 'horizontalalignment', 'center', 'interpreter', 'none', 'fontweight', 'bold', 'fontsize', 12);
            
            obj.string = str;
            obj.y_position = 0.97;
            
            % move to back
            all_axes = obj.h_fig.Children;
            set(obj.h_fig, 'Children', [all_axes(2:end); all_axes(1)])
            axis normal;
        end
        
        
        
        function set.string(obj, val)
        %%set title string to the string 'val'
            obj.set_string(val);
            obj.string = val;
        end
        
        
        
        function set.y_position(obj, val)
        %set height of the title string to 'val' a value between 0 and 1
            obj.set_position(val);
            obj.y_position = val;
        end
        
        
        
        function set_string(obj, val)
        % set title string to 'val'
            obj.h_text.String = val;
        end
        
        
        
        function set_position(obj, val)
        % set height of the title string to 'val' a value between 0 and 1
            obj.h_text.Position = [0.5, val];
        end
    end
end