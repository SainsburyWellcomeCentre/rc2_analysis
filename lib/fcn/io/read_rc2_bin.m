function [data, dt, chan_names, config] = read_rc2_bin(fname_bin)
% READ_RC2_BIN Read data in a .bin file and associated .cfg saved by
% rc2_matlab.
%
%   [DATA, DT, CHANNEL_NAMES, CONFIG] = read_rc2_bin(FILENAME)
%   reads and transform data from a .bin & .cfg pair.
%       Inputs:
%           FILENAME        - string, filename of .bin
%       Outputs:
%           DATA            - # sample points x # channels voltage value
%           DT              - the time difference between sample points
%           CHANNEL_NAMES   - the string IDs given to the channels
%           CONFIG          - the full configuration information contained
%                             in the .cfg file
%
%   Assumes that the .bin is paired with a .cfg in the same directory.
%   Uses the .cfg to obtain # channels and the offsets and
%   scales of each channel.
%   20200304 - currently assumes that the data was saved by transforming
%   voltages between -10 and 10V into 16-bit signed integers.
%
%   See also: read_rc2_config, read_bin, data_transform

% Assume .cfg exists in same directory
fname_cfg       = regexprep(fname_bin, '.bin\>', '.cfg');
% Parse the .cfg
config          = read_rc2_config(fname_cfg);

% Number of channels
chan_names      = config.nidaq.ai.channel_names;
n_channels      = length(chan_names);

% Read the binary data
data            = read_bin(fname_bin, n_channels);

% Use config file to determine offset and scale of each channel
offsets         = config.nidaq.ai.offset;
scales          = config.nidaq.ai.scale;

% Transform the data.
%   Convert 16-bit integer to volts between -10 and 10.
data            = -10 + 20*(data + 2^15)/2^16;
%   Use offset and scale to transform to correct units (cm/s etc.)
data            = data_transform(data, offsets, scales);

% Also return sampling period
dt              =  1/config.nidaq.ai.rate;
