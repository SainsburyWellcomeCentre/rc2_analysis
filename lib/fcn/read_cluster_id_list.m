function ids = read_cluster_id_list(fname)

format_spec = '%f%[^\n\r]';

fid = fopen(fname, 'r');

data = textscan(fid, format_spec, 'texttype', 'string', ...
    'headerlines', 0, 'returnonerror', false);

fclose(fid);

ids = data{1};
