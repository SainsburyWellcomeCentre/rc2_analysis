function mm_pos = normpos2mmpos(norm_pos)
% NORMPOS2MMPOS Convert a normalized position on an A4 figure to a position
% in mm
%
%   MM_POSITION = mmpos2normpos(NORMALIZED_POSITION) takes a box position
%   in normalized coordinates on an A4 figure (NORMALIZED_POSITION) of the
%   form [from_left, from_right, width, height] and where the coordinates
%   run from 0 to 1 (left to right, bottom to top). Returns the position of
%   the box in millimeters.
%
%   See also: mmpos2normpos

a4_size = [210, 297];
mm_pos = norm_pos .* a4_size([1, 2, 1, 2]);
