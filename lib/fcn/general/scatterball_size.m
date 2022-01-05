function dot_area_points = scatterball_size(dot_diameter_mm)
% SCATTERBALL_SIZE Convert a circle with diamter in mm to area in points
%
%   AREA_IN_POINTS = scatterball_size(DIAMETER_IN_MM) converts the diameter
%   of a circle DIAMETER_IN_MM to its area in points AREA_IN_POINTS. This
%   is used with the MATLAB `scatter` function which specifies size as area
%   in points, but more often I find it easier to specify the dot size as
%   having a diamter in mm.

% mm_in_points = 0.352778;
points_per_mm = 72/25.4;
dot_diameter_points = dot_diameter_mm * points_per_mm;
dot_area_points = pi * (dot_diameter_points/2)^2;

% dot_area_mm = pi * (dot_diameter_mm / 2)^2;
% dot_area_in_points = 
% dot_diameter_pts = (dot_diameter_mm/mm_per_point)^2;
% 