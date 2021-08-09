function trace_plot(h_ax, base_t, trace, bouts, id)

hold on;
plot(base_t, trace, 'k');
if ~isempty(bouts)
    for bout_i = 1 : length(bouts)
        idx = bouts(bout_i).start_idx:bouts(bout_i).end_idx;
        plot(h_ax, base_t(idx), trace(idx));
        x_txt = mean(base_t(idx));
        y_txt = max(trace(idx));
        text(x_txt, y_txt, num2str(bout_i));
    end
end
set(gca, 'plotboxaspectratio', [3, 1, 1]);
box off;
if ~isempty(id)
    title(sprintf('Trial %i', id));
end

end