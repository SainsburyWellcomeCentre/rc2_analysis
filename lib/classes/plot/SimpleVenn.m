classdef SimpleVenn < RC2Axis
% SimpleVenn Class for plotting a simple Venn diagram with two groups
%
%   The aim is to make sure the areas of circles and area of overlap are
%   proportional to the number of observations in each of the categories.
%
%  SimpleVenn Properties:
%       A       - number of observations in group A
%       B       - number of observations in group B
%       A_and_B - number of observations in group A and B
%       N       - total number of observations
%       name_A  - the name to give to group A
%       name_B  - the name to give to group B
%       all_col - the colour to give to all observations
%       A_col   - the colour to give to group A observations
%       B_col   - the colour to give to group B observations
%       
%
%  SimpleVenn Methods:
%       plot            - plot the Venn diagram
%       intersect_area  - the area of intersection
%
%   Note: currently requires the Symbolic Math Toolbox in MATLAB

    properties
        
        A
        B
        A_and_B
        N
        
        name_A = 'A'
        name_B = 'B'
        
        all_col = 'k'
        A_col = 'b'
        B_col = 'r'
    end
    
    properties (SetAccess = private)
        
        radius_A
        radius_B
        target_intersect
        
        d_A
        d_B
        
        consistent
        
        h_all_patch
        h_A_patch
        h_B_patch
    end
    
    
    methods
        
        function obj = SimpleVenn(h_ax)
        %%SimpleVenn
        %
        %   SimpleVenn(AXIS_HANDLE) pair object with the axis referrenced
        %   by AXIS_HANDLE. If it is not supplied or empty a new axis will
        %   be created.
        
            VariableDefault('h_ax', []);
            obj = obj@RC2Axis(h_ax);
        end
        
        
        
        function val = get.radius_A(obj)
        %%radius of the circle for group A
            val = sqrt(obj.A/(obj.N*pi));
        end
        
        
        
        function val = get.radius_B(obj)
        %%radius of the circle for group B
            val = sqrt(obj.B/(obj.N*pi));
        end
        
        
        function val = get.target_intersect(obj)    
        %%target
            val = obj.A_and_B/obj.N;
        end
        
        
        function val = get.d_A(obj)
        %%solve for the centre of group A circle
            syms x
            d = real(vpasolve(obj.intersect_area(x) == 0, x));
            val = (obj.radius_A^2 - obj.radius_B^2 + d.^2)./(2*d);
        end
        
        
        function val = get.d_B(obj)
        %%solve for the centre of group B circle    
            syms x
            d = real(vpasolve(obj.intersect_area(x) == 0, x));
            val = (obj.radius_B^2 - obj.radius_A^2 + d.^2)./(2*d);
        end
        
        
        function val = get.consistent(obj)
        %%make sure numbers supplied are consistent   
            val = obj.A_and_B <= obj.A & ...
                obj.A_and_B <= obj.B & ...
                obj.A <= obj.N & ...
                obj.B <= obj.N & ...
                obj.A + obj.B - obj.A_and_B <= obj.N;
        end
        
        function plot(obj, suppress_text)
        %%plot Plot the Venn diagram
        %
        %   plot(SUPPRESS_TEXT) plots the Venn diagram. SUPPRESS_TEXT
        %   indicates whether to prevent plotting of the text (default =
        %   false, i.e. do not suppress text).
        
            VariableDefault('suppress_text', false);
            
            if ~obj.consistent
                error('numbers are not consistent')
            end
            
            t = 0:0.01:2*pi+0.01;
            
            x_A = double(-obj.d_A + obj.radius_A*cos(t));
            y_A = double(obj.radius_A*sin(t));
            
            x_B = double(obj.d_B + obj.radius_B*cos(t));
            y_B = double(obj.radius_B*sin(t));
            
%             figure
%             hold on
            
            
%             xl = get(gca, 'xlim');
%             yl = get(gca, 'ylim');
%             
%             X = range(xl);
%             Y = 1/X;
%             
%             if Y < 2*max(obj.radius_A, obj.radius_B)
%                 xl = double([(-obj.d_A - obj.radius_A), (obj.d_B + obj.radius_B)]);
%                 Y = 1/range(xl);
%             end
%             
%             xlim(xl)
%             ylim([-Y/2, Y/2])
            
            set(obj.h_ax, 'dataaspectratio', [1, 1, 1], 'clipping', 'off')
%             line(xl, -Y([1, 1])/2, 'color', 'k');
%             line(xl, Y([1, 1])/2, 'color', 'k');
%             line(xl([1, 1]), [-Y/2, Y/2], 'color', 'k');
%             line(xl([2, 2]), [-Y/2, Y/2], 'color', 'k');
%             axis off
            
            pAB = ((obj.d_B + obj.radius_B) + (-obj.d_A - obj.radius_A))/2;
            big_circle_x = double(pAB + sqrt(1/pi)*cos(t));
            big_circle_y = double(sqrt(1/pi)*sin(t));
            
            obj.h_all_patch = patch(obj.h_ax, 'xdata', big_circle_x, 'ydata', big_circle_y, 'facecolor', obj.all_col);
            obj.h_A_patch = patch(obj.h_ax, 'xdata', x_A, 'ydata', y_A, 'facecolor', obj.A_col);
            obj.h_B_patch = patch(obj.h_ax, 'xdata', x_B, 'ydata', y_B, 'facecolor', obj.B_col);
            
            pA = ((obj.d_B - obj.radius_B) + (-obj.d_A - obj.radius_A))/2;
            pB = ((-obj.d_A + obj.radius_A) + (obj.d_B + obj.radius_B))/2;
            pAB = ((obj.d_B - obj.radius_B) + (-obj.d_A + obj.radius_A))/2;
            
            if ~suppress_text
                text(obj.h_ax, pA, 0, sprintf('%s only = %i', obj.name_A, obj.A - obj.A_and_B), 'verticalalignment', 'middle', 'horizontalalignment', 'center');
                text(obj.h_ax, pB, 0, sprintf('%s only = %i', obj.name_B, obj.B - obj.A_and_B), 'verticalalignment', 'middle', 'horizontalalignment', 'center');
                text(obj.h_ax, pAB, 0, sprintf('%i', obj.A_and_B), 'verticalalignment', 'middle', 'horizontalalignment', 'center');
    %             text(xl(2), Y/2, sprintf('N = %i', obj.N), 'verticalalignment', 'top', 'horizontalalignment', 'right');
                text(obj.h_ax, pAB, sqrt(1/pi), sprintf('N = %i', obj.N), 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
            end
        end
        
        
        function A = intersect_area(obj, d)
        %%intersect_area Intersection area of the two circles depending on
        %%distance of centres from each other
        %   solve this to position the two circles so that we get the
        %   correct overlap area
        %
        %   intersect_area(DISTANCE) computes overlap area of two circles
        %   DISTANCE apart.
        
            rA = obj.radius_A;
            rB = obj.radius_B;
            
            dA = (rA^2 - rB^2 + d.^2)./(2*d);
            thA = acos(dA/rA);
            AA = thA*rA^2 - dA.*sqrt(rA^2 - dA.^2);
            
            dB = d - dA;
            thB = acos(dB/rB);
            AB = thB*rB^2 - dB.*sqrt(rB^2 - dB.^2);
            
            A = AA + AB - obj.target_intersect;
        end
    end
end