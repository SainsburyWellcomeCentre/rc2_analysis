function trace_decreases = does_mismatch_running_decrease(t, traces, window_1, window_2, n_sds)

baseline_idx = t >= window_1(1) & t < window_1(2);
response_idx = t >= window_2(1) & t < window_2(2);

baseline_avg = mean(traces(baseline_idx, :), 1);
baseline_std = std(traces(baseline_idx, :), 1);
response_avg = mean(traces(response_idx, :), 1);

trace_decreases = response_avg < (baseline_avg - n_sds * baseline_std);
