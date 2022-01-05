function h_ax = a4axis(h_fig, axis_position)
% A4AXIS Create an axis on an A4 figure
%
%   AXIS_HANDLE = a4axis(FIGURE_HANDLE, POSITION) creates an axis on the
%   A4 figure with handle FIGURE_HANDLE in position specified by POSITION.
%   POSITION is a 1x4 array with a position on the A4 page specified in
%   millimeters [mm_from_left, mm_from_right, mm_width, mm_height].
%
%   If FIGURE_HANDLE is not supplied or empty, a new a4figure is created.
%   If POSITION is not supplied or empty, an axis with postion [65, 123.5,
%   80, 50] is created in the middle of the A4 figure.
%
%   The handle to the axis is returned in AXIS_HANDLE.
%
%   See also: a4figure

VariableDefault('h_fig', []);
VariableDefault('axis_position', [65, 123.5, 80, 50]);

if isempty(h_fig)
    h_fig = a4figure();
end

h_ax = axes(h_fig);

set(h_ax, ...
    'Color',                    'w', ...
    'XColor',                   'k', ...
    'YColor',                   'k', ...
    'Units',                    'normalized', ...
    'PositionConstraint',       'innerposition', ...
    'InnerPosition',            mmpos2normpos(axis_position), ...
    'PlotBoxAspectRatioMode',   'auto');

hold on;
