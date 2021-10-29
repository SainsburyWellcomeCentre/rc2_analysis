function anatomy = format_anatomy(ks_dir)

% read probe file
probe_track       = read_track_csv(ks_dir);

% correct probe track
if ~isempty(probe_track)
    track_offset      = read_track_offset(ks_dir);
    probe_track.pos   = probe_track.pos + track_offset;
else
    anatomy.region_boundaries = [];
    anatomy.region_id = [];
    anatomy.region_str = [];
    return
end

% find boundaries (um from probe tip)
[boundaries, id, str] = find_layer_boundaries(probe_track);

% store
anatomy.region_boundaries = boundaries;
anatomy.region_id = id;
anatomy.region_str = str;
