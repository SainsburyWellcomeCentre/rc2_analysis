function mm_pos = normpos2mmpos(norm_pos)

a4_size = [210, 297];
mm_pos = norm_pos .* a4_size([1, 2, 1, 2]);
