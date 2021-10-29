function [boundaries, id, str] = find_layer_boundaries(probe_track)

pos = probe_track.pos;
region_id = probe_track.region_id;
region_str = probe_track.region_str;

region_id(isnan(region_id)) = -1;

idx = find(abs(diff(region_id)) > 0);
boundaries = [pos(1), pos(idx+1), pos(end)];

id = [region_id(1), region_id(idx+1)];
str = [region_str(1), region_str(idx+1)];
