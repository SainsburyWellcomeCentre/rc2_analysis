function h_text = a4figure_text(text_str, h_fig, text_position)

h_ax = a4axis(h_fig, [0, 0, 210, 297]);
axis off;
norm_pos = mmpos2normpos(text_position);
h_text = text(norm_pos(1), norm_pos(2), text_str, 'horizontalalignment', 'left', 'verticalalignment', 'bottom');