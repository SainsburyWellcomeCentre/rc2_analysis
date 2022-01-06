function trace_decreases = does_mismatch_running_decrease(t, traces, window_1, window_2, n_sds)
% DOES_MISMATCH_RUNNING_DECREASE Compare the running traces in two windows
% and determine if there is a significant decrease in running.
%
%   IS_DECREASE = does_mismatch_running_decrease(TIMEBASE, TRACES, WINDOW_1, WINDOW_2, NUMBER_SDS)
%   takes a set of running traces in TRACES (#samples x #traces) with
%   timebase TIMEBASE (#samples x 1), and looks at the average and S.D. of the speed in
%   WINDOW_1 (a 1x2 array defining the time in which to look at the speed).
%   It then calcuates the average speed in WINDOW_2 (also a 1x2 array
%   defining the time in which to look at the speed) for each trace. It
%   then compares the average speed in WINDOW_2 to the (average in WINDOW_1) -
%   NUMBER_SDS*(standard deviation of speed in WINDOW_1) and if it is lower
%   consider this to be a decrease in speed.
%
%   IS_DECREASE is a 1 x #traces boolean array with true if the trace
%   decreased.

baseline_idx = t >= window_1(1) & t < window_1(2);
response_idx = t >= window_2(1) & t < window_2(2);

baseline_avg = mean(traces(baseline_idx, :), 1);
baseline_std = std(traces(baseline_idx, :), 1);
response_avg = mean(traces(response_idx, :), 1);

trace_decreases = response_avg < (baseline_avg - n_sds * baseline_std);
