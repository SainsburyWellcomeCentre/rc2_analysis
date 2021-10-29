function [lfp_power, p, channel_list, params] = get_lfp_power(bin_fname)
%%lfp_power = GET_LFP_POWER(bin_fname)
%   
%   Gets the LFP power between 500 and 5000 (or max from sample rate), on
%   each channel of an .imec.lf.bin recording.
%   
%   Assumes that the trigger channel has been saved.
%       
%   Reads the .imec.lf.meta file for channel info and sample rate, so this
%   file must exist.
%
%   Assumes a set of 'bad_channels' which are interpolated over at the
%   end... this needs to become an input argument or automatically
%   determined.

% number of channels recorded (incl. trigger)
meta_fname          = regexprep(bin_fname, '.bin\>', '.meta');
config              = read_spikeglx_config(meta_fname);
n_channels          = str2double(config.nSavedChans);
fs                  = str2double(config.imSampRate);

%TODO: needs to become parameter or automatic detection
% bad_chans = [1, 27, 37, 70, 72, 76, 98, 105, 113, 124, 152, 155, 156, 157, 168, 189];
bad_chans = [1, 2, 3, 37, 76, 113, 152, 189];

% create memory map to lfp file
finfo = dir(bin_fname);
n_samples = finfo.bytes/(2*n_channels);
mmf = memmapfile(bin_fname, 'format', {'int16', [n_channels, n_samples], 'x'});

% where to assess the power
offset_s = 60:5*60:(n_samples/fs);
duration_s = 10;

% get the power on each channel between 500 and
n_channels_to_proc = n_channels - 1;
offset_samps = offset_s * fs;
duration_samps = duration_s * fs;

% use 'bad_channels' to get list of good channels
good_chans = 1:n_channels_to_proc;
good_chans(bad_chans) = [];


p = nan(n_channels_to_proc, length(offset_samps));

for o_idx = 1 : length(offset_samps)
    
    o_idx
    
    v = offset_samps(o_idx) + (0:duration_samps-1)';
    if v(end) > n_samples
        p(:, end) = [];
        continue
    end
    for i = 1 : n_channels_to_proc
        p(i, o_idx) = bandpower(single(mmf.Data.x(i, v)), fs, [500, min(fs/2, 5000)]);
    end
    
    % interpolate the power over bad channels
    p(:, o_idx) = smooth(interp1(good_chans, p(good_chans, o_idx), 1:length(p)), 10);
end
lfp_power = mean(p, 2);


channel_list = 1:n_channels_to_proc;

params.offset_s = offset_s;
params.duration_s = duration_s;
params.bad_channels = bad_chans;
params.bin_fname = bin_fname;
