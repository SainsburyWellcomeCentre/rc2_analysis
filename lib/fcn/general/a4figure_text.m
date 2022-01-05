function h_text = a4figure_text(text_str, h_fig, text_position)
% A4FIGURE_TEXT Create text at a particular position on an A4 figure
%
%   TEXT_HANDLE = a4figure_text(STRING, FIGURE_HANDLE, POSITION) creates
%   text in STRING at a position POSITION on a A4 figure with handle FIGURE_HANDLE.
%   It does this by creating an a4axis across the entire A4 figure, and
%   then placing the text in this axis. The handle to the text is returned
%   in TEXT_HANDLE.
%
%   Text is left aligned horizontally, and bottom aligned vertically.
%
%   See also: a4axis, a4figure
%
%   TODO:  order of arguments is strange.

h_ax = a4axis(h_fig, [0, 0, 210, 297]);
axis off;
norm_pos = mmpos2normpos(text_position);
h_text = text(norm_pos(1), norm_pos(2), text_str, 'horizontalalignment', 'left', 'verticalalignment', 'bottom');