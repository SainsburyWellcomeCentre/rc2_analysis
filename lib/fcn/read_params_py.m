function data = read_params_py(fname)

% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
format_spec = '%s%s%s%[^\n\r]';

%% Open the text file.
fid = fopen(fname,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
text_str = textscan(fid, format_spec, inf, ...
    'delimiter', ' ', ...
    'multipledelimsasone', true, ...
    'texttype', 'string', ...
    'headerlines', 0, ...
    'returnonerror', false, ...
    'endofline', '\r\n');

for i = 1 : length(text_str{1})
    data.(text_str{1}{i}) = text_str{3}{i};
end

data.n_channels_dat = str2double(data.n_channels_dat);
data.offset = str2double(data.offset);
data.sample_rate = str2double(data.sample_rate);

% close the text file.
fclose(fid);
