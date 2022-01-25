function out = zscore(timebase, trace, baseline_window)
%%ZSCORE Computes the z-score of a trace
%
%   ZSCORE_TRACE = zscore(TIMEBASE, TRACE, BASELINE_WINDOW)
%   takes a trace, TRACE, a Nx1 vector of values with associated timebase
%   TIMEBASE (also a Nx1 vector specifying the time in seconds of each
%   sample point in TRACE) and computes a z-score by first calcuating the
%   mean and standard deviation of the trace in the time window given by
%   BASELINE_WINDOW (a 1x2 vector specifying the start and end points of
%   the window to take). The trace is then computed by
%
%       ZSCORE_TRACE = (TRACE - mean_baseline) / std_baseline

idx = timebase >= baseline_window(1) & ...
      timebase < baseline_window(2);
 
baseline_avg = mean(trace(idx));
baseline_std = std(trace(idx));

out = (trace - baseline_avg) / baseline_std;
