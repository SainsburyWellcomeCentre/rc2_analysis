function offset = read_track_offset(ks_dir)

correction_fname = fullfile(ks_dir, 'tracks', 'offset.txt');

% if probe track file doesn't exist return empty array
if ~exist(correction_fname, 'file')
    offset = nan;
    warning('No offset file called %s', correction_fname);
    return
end

fid = fopen(correction_fname, 'r');
offset = textscan(fid, '%f%[^\n\r]');
offset = offset{1};
fclose(fid);
