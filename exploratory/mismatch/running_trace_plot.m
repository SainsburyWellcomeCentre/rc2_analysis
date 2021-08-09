function running_trace_plot(h_ax, t, traces)

yM                      = 60;
n_sample_points         = length(t);
downsample_by           = 10;
idx_to_show             = 1:downsample_by:n_sample_points;

hold on
plot(h_ax, t(idx_to_show), traces(idx_to_show, :), 'color', [0.6, 0.6, 0.6]);
avg = mean(traces, 2);
plot(h_ax, t(idx_to_show), avg(idx_to_show), 'k');
ylim(h_ax, [0, yM]);
xlim(h_ax, [-1, 1]);
line(h_ax, [0, 0], [0, yM], 'color', 'k')
text(h_ax, 0, yM, sprintf('n = %i', size(traces, 2)), ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
box off
