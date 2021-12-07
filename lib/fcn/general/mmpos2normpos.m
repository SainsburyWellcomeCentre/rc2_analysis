function norm_pos = mmpos2normpos(pos)

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
