function spiking_trace_plot(h_ax, t, traces)

hold on;
yL                      = [-4, 6];
n_sample_points         = length(t);
downsample_by           = 10;
idx_to_show             = 1:downsample_by:n_sample_points;
bsl                     = t > -1 & t < 0;

m_rm                    = bsxfun(@minus, traces, mean(traces(bsl, :), 1));
m                       = nanmean(m_rm, 2)';
s                       = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));

plot(h_ax, t(idx_to_show), m(idx_to_show), 'r');
plot(h_ax, t(idx_to_show), m(idx_to_show)-s(idx_to_show), 'color', [1, 0.2, 0.2])
plot(h_ax, t(idx_to_show), m(idx_to_show)+s(idx_to_show), 'color', [1, 0.2, 0.2])

ylim(h_ax, yL);
line(h_ax, [0, 0], yL, 'color', 'k');
box off;
