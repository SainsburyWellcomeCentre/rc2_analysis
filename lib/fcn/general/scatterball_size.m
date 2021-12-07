function dot_area_points = scatterball_size(dot_diameter_mm)
%%dot_diameter_pts = scatterball_size(dot_diameter_mm)
%   Get the dot diameter for a scatter plot in points (convert from mm)

% mm_in_points = 0.352778;
points_per_mm = 72/25.4;
dot_diameter_points = dot_diameter_mm * points_per_mm;
dot_area_points = pi * (dot_diameter_points/2)^2;

% dot_area_mm = pi * (dot_diameter_mm / 2)^2;
% dot_area_in_points = 
% dot_diameter_pts = (dot_diameter_mm/mm_per_point)^2;
% 