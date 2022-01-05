function trace = filter_trace(trace, fs, cut_off, order)
% FILTER_TRACE Low-pass filter the velocity traces
%
%   FILTERED_TRACE = filter_trace(TRACE, FS, CUT_OFF, ORDER)
%   forward and backward filters TRACE (with sample rate FS Hz) with lowpass
%   Butterworth filter with cut off CUT_OFF and order ORDER.
%   Default FS = 10000 Hz, cut_off = 50 Hz, order = 3
%
%   See also: butter

% Filter
VariableDefault('fs', 10000);
VariableDefault('cut_off', 50);
VariableDefault('order', 3);

Wn = cut_off/(fs/2);

[b, a] = butter(order, Wn, 'low');

trace = filter(b, a, trace);
trace = filter(b, a, trace(end:-1:1));
trace = trace(end:-1:1);
