function trace = filter_trace(trace, fs, cut_off, order)
% forward and backward filters trace with sample rate of 'fs' with lowpass
% filter with cut off 'cut_off' and order 'order'.
%   Default fs = 10000, cut_off = 50, order = 3

%% Filter
VariableDefault('fs', 10000);
VariableDefault('cut_off', 50);
VariableDefault('order', 3);

Wn = cut_off/(fs/2);

[b, a] = butter(order, Wn, 'low');

trace = filter(b, a, trace);
trace = filter(b, a, trace(end:-1:1));
trace = trace(end:-1:1);
