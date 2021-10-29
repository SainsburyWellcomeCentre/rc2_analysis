function probe_track = read_track_csv(ks_dir, track_fname)

VariableDefault('track_fname', 'track_0');

track_fname = fullfile(ks_dir, 'tracks', [track_fname, '.csv']);

% if probe track file doesn't exist return empty array
if ~exist(track_fname, 'file')
    probe_track = [];
    warning('No track file named %s', track_fname);
    return
end


fid = fopen(track_fname, 'r');
a = fread(fid);
a = char(a');
fclose(fid);

b = strsplit(a, '\n');

for i = 2 : length(b)-1
    c = strsplit(b{i}, ',');
    probe_track.pos(i-1) = str2double(c{1});
    probe_track.region_id(i-1) = str2double(c{2});
    probe_track.region_str{i-1} = c{3};
end

% swap for ease
probe_track.pos = probe_track.pos(end)-probe_track.pos;
