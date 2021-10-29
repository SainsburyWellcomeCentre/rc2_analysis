function data = read_bin(fname, n_chan)
%data = READ_BIN(fname, n_chan)
% Reads binary data saved by rc2 (int16)
%   Channels along the columns, samples along the rows.

fid = fopen(fname, 'r');
data = fread(fid, 'int16');
data = reshape(data, n_chan, [])';
fclose(fid);