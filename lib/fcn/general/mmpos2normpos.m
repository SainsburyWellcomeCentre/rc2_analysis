function norm_pos = mmpos2normpos(pos)
% MMPOS2NORMPOS Convert a position in mm to a normalized position for A4 figure
%
%   NORMALIZED_POSITION = mmpos2normpos(MM_POSITION) takes a box position
%   MM_POSITION of the form [mm_from_left, mm_from_right, mm_width,
%   mm_height] on an A4 figure, and converts it to a normalized position
%   NORMALIZED_POSITION where the X and Y dimensions of the figure run from
%   0 to 1 (left to right, or bottom to top).
%
%   See also: normpos2mmpos
%
%   TODO: this can be greatly simplified (see normpos2mmpos).

a4_size = [210, 297];

from_left_mm = pos(1);
from_bottom_mm = pos(2);
width_mm = pos(3);
height_mm = pos(4);

from_left_norm = from_left_mm / a4_size(1);
from_bottom_norm = from_bottom_mm / a4_size(2);
width_norm = width_mm / a4_size(1);
height_norm = height_mm / a4_size(2);

norm_pos = [from_left_norm, from_bottom_norm, width_norm, height_norm];
