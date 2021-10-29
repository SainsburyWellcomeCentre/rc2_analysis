function data = read_cluster_groups(fname)

delimiter = '\t';
start_row = 2;

%% Format for each line of text:
%   column1: double (%f)
%	column2: string (%s)
% For more information, see the TEXTSCAN documentation.
format_spec = '%f%s%[^\n\r]';

fid = fopen(fname,'r');

data = textscan(fid, format_spec, 'delimiter', delimiter, 'texttype', 'string', ...
    'headerlines', start_row-1, 'returnonerror', false, 'endofline', '\r\n');

% remove last entry
data = cat(2, data{1:2});

%% Close the text file.
fclose(fid);
