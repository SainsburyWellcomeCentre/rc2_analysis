function bout_plot(h_ax, base_t, trace, bout)

pad = 5000;
hold on;
idx = (bout.start_idx-pad):(bout.end_idx+pad);
plot(base_t(idx)-base_t(bout.start_idx), trace(idx), 'color', 'k');
set(h_ax, 'plotboxaspectratio', [3, 1, 1]);
box off;
