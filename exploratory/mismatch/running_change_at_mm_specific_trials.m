% 

clear all
close all

probe_id = [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5];
trial_id = [8, 33, 15, 24, 28, 38, 46, 1, 32, 32, 51];

window_t                  = 0.2;
n_sds                     = 2;
window_1                  = -0.15 + [0, window_t];
window_2                  = 0.1 + [0, window_t];
display_window            = [-1, 1];
protocol                  = 2;

%%
config                    = RC2AnalysisConfig();

protocols                 = MismatchExperiment.protocol_ids;
protocol_labels           = MismatchExperiment.protocol_label;

figs                      = RC2Figures(config);
figs.save_on              = true;
figs.set_figure_subdir('mismatch', 'running_change_at_mm', 'specific_trials');

probe_fnames              = experiment_details('mismatch_nov20', 'protocol');

running                   = cell(length(probe_fnames));

store_running_down        = [];
store_running_remaining   = [];
store_running_all         = [];
store_spikes_down         = [];
store_spikes_remaining    = [];
store_spikes_all          = [];
cluster_c                 = 0;

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.VISp_clusters;
    
    exp_obj             = MismatchExperiment(data, config);
    
    for prot_i = protocol%1 : length(protocols)
        
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
            
            if any((probe_id == probe_i) & (trial_id == trials(trial_i).id))
                is_unchanged = true;
            else
                is_unchanged = false;
            end
            
%             changed_up_idx = any(m_after > m_before + n_sds*sd);
%             changed_down_idx = any(m_after < m_before - n_sds*sd);
            
%             line(h_ax, window_1, m_before([1, 1]), 'color', 'k')
%             line(h_ax, [window_1(1), window_2(2)], m_before([1, 1]) + n_sds*sd, 'color', [0.5, 0.5, 0.5])
%             line(h_ax, [window_1(1), window_2(2)], m_before([1, 1]) - n_sds*sd, 'color', [0.5, 0.5, 0.5])
            
%             line(h_ax, window_2, m_after([1, 1]), 'color', 'k')
            
            if ~is_unchanged
                plot(h_ax, t, running(:, trial_i), 'r');
            else
                plot(h_ax, t, running(:, trial_i), 'k');
            end
            
            if ~is_unchanged
                store_running_down = [store_running_down, running(:, trial_i)];
                trial_class(trial_i) = 1;
            else
                store_running_remaining = [store_running_remaining, running(:, trial_i)];
                trial_class(trial_i) = 2;
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
            store_spikes_down = [store_spikes_down, mean(spike_rate{cluster_i}(:, trial_class == 1), 2)];
            store_spikes_remaining = [store_spikes_remaining, mean(spike_rate{cluster_i}(:, trial_class == 2), 2)];
            store_spikes_all = [store_spikes_all, mean(spike_rate{cluster_i}, 2)];
        end
        
        FigureTitle(gcf, sprintf('%s, %s', probe_fnames{probe_i}, protocol_labels{prot_i}));
        
        figs.save_fig_to_join();
    end
    
    figs.join_figs(sprintf('%s.pdf', probe_fnames{probe_i}));
    figs.clear_figs();
end



%% PLOT AVERAGES
n_down = size(store_running_down, 2);
yM = 40;

h_fig                   = figs.a4figure();

subplot(2, 3, 1)
hold on
if n_down > 0
    m = mean(store_running_down, 2)';
    s = std(store_running_down, [], 2)';
    fill([t, t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], [0.7, 0.7, 0.7]);
    plot(t, m, 'k');
end
ylim([0, yM]);
xlim([-1, 1]);
line([0, 0], [0, yM], 'color', 'k')
text(0, yM, sprintf('n = %i', n_down), 'verticalalignment', 'top', 'horizontalalignment', 'left');
xlabel('Time from MM onset (s)')
ylabel('Running (cm/s)')
title('Negative change trials');
box off

subplot(2, 3, 2)
hold on
m = mean(store_running_remaining, 2)';
s = std(store_running_remaining, [], 2)';
fill([t, t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], [0.7, 0.7, 0.7]);
plot(t, m, 'k');
ylim([0, yM]);
line([0, 0], [0, yM], 'color', 'k')
text(0, yM, sprintf('n = %i', size(store_running_remaining, 2)), 'verticalalignment', 'top', 'horizontalalignment', 'left');
title('Remaining trials');
box off

subplot(2, 3, 3)
hold on
m = mean(store_running_all, 2)';
s = std(store_running_all, [], 2)';
fill([t, t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], [0.7, 0.7, 0.7]);
plot(t, m, 'k');
ylim([0, yM]);
line([0, 0], [0, yM], 'color', 'k');
text(0, yM, sprintf('n = %i', size(store_running_all, 2)), 'verticalalignment', 'top', 'horizontalalignment', 'left');
title('All trials');
box off

FigureTitle(gcf, 'Average running speed around MM onset');
figs.save_fig('averages.pdf');


n_down = size(store_spikes_down, 2);
yL = [-4, 6];
bsl = t > -1 & t < 0;
m_rm = bsxfun(@minus, store_spikes_all, mean(store_spikes_all(bsl, :), 1));
m_all = nanmean(m_rm, 2)';
s_all = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));


subplot(2, 3, 4)
hold on
if n_down > 0
    
    m_rm = bsxfun(@minus, store_spikes_down, mean(store_spikes_down(bsl, :), 1));
    
    m = nanmean(m_rm, 2)';
    s = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));
    
    h = fill([t, t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], 'r');
    set(h, 'facealpha', 0.6);
    plot(t, m, 'r');
end
ylim(yL);
line([0, 0], yL, 'color', 'k');
ylabel('\Delta Hz');
title('Negative change trials');
box off


subplot(2, 3, 5)
hold on
m_rm = bsxfun(@minus, store_spikes_remaining, mean(store_spikes_remaining(bsl, :), 1));
m = nanmean(m_rm, 2)';
s = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));
h = fill([t, t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], 'r');
set(h, 'facealpha', 0.6);
plot(t, m, 'r');
ylim(yL);
line([0, 0], yL, 'color', 'k');
title('Remaining trials');
box off


subplot(2, 3, 6)
hold on
fill([t, t(end:-1:1)], [m_all-s_all, m_all(end:-1:1)+s_all(end:-1:1)], [0.7, 0.7, 0.7]);
plot(t, m_all, 'k');
ylim(yL);
line([0, 0], yL, 'color', 'k');
title('All trials');
box off;


FigureTitle(gcf, 'Average running speed around MM onset');
figs.save_fig('averages_with_spikes.pdf');





