clear all
close all

window_t                = 0.2;
n_sds                   = 2;
window_1                = -0.15 + [0, window_t];
window_2                = 0.1 + [0, window_t];
display_window          = [-1, 1];
protocol                = 2;


%%
config                  = RC2AnalysisConfig();

protocols               = MismatchExperiment.protocol_ids;
protocol_labels         = MismatchExperiment.protocol_label;

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('mismatch', datestr(now, 'yyyymmdd'), 'running_change_at_mm', sprintf('%ims_%isd_offset_protocol_%i', 1e3*window_t, n_sds, protocol));

probe_fnames            = experiment_details('mismatch_nov20', 'protocol');

running                 = cell(length(probe_fnames));

store_running_up          = [];
store_running_down        = [];
store_running_none        = [];
store_running_all         = [];
store_spikes_up           = [];
store_spikes_none         = [];
store_spikes_down         = [];
store_spikes_all          = [];
cluster_c                 = 0;

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.VISp_clusters;
    
    exp_obj             = MismatchExperiment(data, config);
    
    for prot_i = protocol
        
        %%
        h_fig                   = figs.a4figure();
        
        [running, t] = exp_obj.running_around_mismatch_by_protocol(prot_i, display_window);
        
        spike_rate = {};
        for cluster_i = 1 : length(clusters)
            spike_rate{cluster_i} = exp_obj.firing_around_mismatch_by_protocol(clusters(cluster_i), prot_i, display_window);
        end
        
        trials = exp_obj.trials_of_type(prot_i);
        trial_class = nan(length(trials), 1);
        
        for trial_i = 1 : size(running, 2)
            
            h_ax = subplot(5, 5, trial_i);
            hold on;
            
            baseline_idx = t > window_1(1) & t < window_1(2);
            response_idx = t >= window_2(1) & t < window_2(2);
            
            m_before = mean(running(baseline_idx, trial_i));
            sd = std(running(baseline_idx, trial_i));
            
            m_after = mean(running(response_idx, trial_i));
            
            changed_up_idx = any(m_after > m_before + n_sds*sd);
            changed_down_idx = any(m_after < m_before - n_sds*sd);
            
            line(h_ax, window_1, m_before([1, 1]), 'color', 'k')
            line(h_ax, [window_1(1), window_2(2)], m_before([1, 1]) + n_sds*sd, 'color', [0.5, 0.5, 0.5])
            line(h_ax, [window_1(1), window_2(2)], m_before([1, 1]) - n_sds*sd, 'color', [0.5, 0.5, 0.5])
            
            line(h_ax, window_2, m_after([1, 1]), 'color', 'k')
            
            
            if changed_up_idx | changed_down_idx
                plot(h_ax, t, running(:, trial_i), 'r');
            else
                plot(h_ax, t, running(:, trial_i), 'k');
            end
            
            if changed_up_idx
                store_running_up = [store_running_up, running(:, trial_i)];
                trial_class(trial_i) = 1;
            elseif changed_down_idx
                store_running_down = [store_running_down, running(:, trial_i)];
                trial_class(trial_i) = 2;
            else
                store_running_none = [store_running_none, running(:, trial_i)];
                trial_class(trial_i) = 3;
            end
            
            store_running_all = [store_running_all, running(:, trial_i)];
            
            line([0, 0], get(gca, 'ylim'), 'color', 'k')
            
            if trial_i == 1
                xlabel('Time from MM onset (s)')
                ylabel('Running (cm/s)')
            end
            
            box off
            
            title(sprintf('Trial ID %i', trials(trial_i).id));
        end
        
        
        for cluster_i = 1 : length(clusters)
%             store_spikes_up = [store_spikes_up, mean(spike_rate{cluster_i}(:, trial_class == 1), 2)];
            store_spikes_down = [store_spikes_down, mean(spike_rate{cluster_i}(:, trial_class == 2), 2)];
            store_spikes_none = [store_spikes_none, mean(spike_rate{cluster_i}(:, ismember(trial_class, [1, 3])), 2)];
            store_spikes_all = [store_spikes_all, mean(spike_rate{cluster_i}, 2)];
        end
        
        FigureTitle(gcf, sprintf('%s, %s', probe_fnames{probe_i}, protocol_labels{prot_i}));
        
%         figs.save_fig_to_join();
    end
    
%     figs.join_figs(sprintf('%s.pdf', probe_fnames{probe_i}));
%     figs.clear_figs();
end

% combine positive
store_running_none = [store_running_none, store_running_up];




%% PLOT AVERAGES
yM              = 60;
grey_col        = 0.6;
h_fig           = figs.a4figure();


subplot(2, 4, 1)
hold on
plot(t, store_running_down, 'color', grey_col([1, 1, 1]));
plot(t, mean(store_running_down, 2), 'k', 'linewidth', 2);
ylim([0, yM]);
line([0, 0], [0, yM], 'color', 'k')
text(0, yM, sprintf('n = %i', size(store_running_down, 2)), 'verticalalignment', 'top', 'horizontalalignment', 'left');
title('Negative change trials');
box off

subplot(2, 4, 2)
hold on
plot(t, store_running_none, 'color', grey_col([1, 1, 1]));
plot(t, mean(store_running_none, 2), 'k', 'linewidth', 2);
ylim([0, yM]);
line([0, 0], [0, yM], 'color', 'k');
text(0, yM, sprintf('n = %i', size(store_running_none, 2)), 'verticalalignment', 'top', 'horizontalalignment', 'left');
title('No change trials');
box off

subplot(2, 4, 3)
hold on
plot(t, store_running_all, 'color', grey_col([1, 1, 1]));
plot(t, mean(store_running_all, 2), 'k', 'linewidth', 2);
ylim([0, yM]);
line([0, 0], [0, yM], 'color', 'k');
text(0, yM, sprintf('n = %i', size(store_running_all, 2)), 'verticalalignment', 'top', 'horizontalalignment', 'left');
title('All trials');
box off

FigureTitle(gcf, 'Average running speed around MM onset');
figs.save_fig('averages.pdf');



yL = [-4, 6];
bsl = t > -1 & t < 0;
m_rm = bsxfun(@minus, store_spikes_all, mean(store_spikes_all(bsl, :), 1));
m_all = nanmean(m_rm, 2)';
s_all = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));

subplot(2, 4, 5)
hold on
m_rm = bsxfun(@minus, store_spikes_down, mean(store_spikes_down(bsl, :), 1));
m = nanmean(m_rm, 2)';
s = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));
h = fill([t, t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], 'r');
set(h, 'facealpha', 0.6);
plot(t, m, 'r');
ylim(yL);
line([0, 0], yL, 'color', 'k');
title('Negative change trials');
box off

subplot(2, 4, 6)
hold on
m_rm = bsxfun(@minus, store_spikes_none, mean(store_spikes_none(bsl, :), 1));
m = nanmean(m_rm, 2)';
s = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));
h = fill([t, t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], 'r');
set(h, 'facealpha', 0.6);
plot(t, m, 'r');
ylim(yL);
line([0, 0], yL, 'color', 'k');
title('No change trials');
box off

subplot(2, 4, 7)
hold on
fill([t, t(end:-1:1)], [m_all-s_all, m_all(end:-1:1)+s_all(end:-1:1)], [0.7, 0.7, 0.7]);
plot(t, m_all, 'k');
ylim(yL);
line([0, 0], yL, 'color', 'k');
title('All trials');
box off;

FigureTitle(gcf, 'Average running speed around MM onset');
figs.save_fig('averages_with_spikes.pdf');
