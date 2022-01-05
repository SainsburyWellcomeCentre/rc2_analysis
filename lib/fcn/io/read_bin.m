function data = read_bin(fname, n_chan)
% READ_BIN Reads a .bin file with int16 data stored
%
%   DATA = read_bin(FILENAME, NUMBER_OF_CHANNELS) reads int16
%   data in a .bin file. FILENAME is the full path to the .bin file and
%   NUMBER_OF_CHANNELS specifies the number of channels which were saved.
%   The int16 data is returned in DATA, with channels along the columns, sample points along the rows.

fid = fopen(fname, 'r');
data = fread(fid, 'int16');
data = reshape(data, n_chan, [])';
fclose(fid);
