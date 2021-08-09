% For each mismatch trial, takes the running speed around the mismatch and
% looks for periods most similar to that running speed in the period before
% the mismatch.

clear all
close all

window_t                = 0.2;
n_sds                   = 2;
window_1                = -0.15 + [0, window_t];
window_2                = 0.1 + [0, window_t];
window_t                = [0, 0.15];
display_window          = [-1, 1];
protocol                = 2;
use_average_trace       = true;  % use an average trace for template matching or each mismatch trial?



%%
config                  = RC2AnalysisConfig();

protocols               = MismatchExperiment.protocol_ids;
protocol_labels         = MismatchExperiment.protocol_label;

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('mismatch', 'running_change_at_mm', 'trial_matched');

probe_fnames            = experiment_details('mismatch_nov20', 'protocol');

running                 = cell(length(probe_fnames));

store_running_mm        = [];
store_running_matched   = [];
store_spikes_mm         = [];
store_spikes_matched    = [];
store_details           = [];

fname = '1s';
if use_average_trace
    
    load('avg_to_compare_1s.mat', 'avg_to_compare');
    
    % template velocity
    template_velocity = avg_to_compare(:);
    
    % the number of samples to compare on each window
    n_samples_to_compare = length(template_velocity);
    
    % minimum amplitude we will attempt to find
    template_amplitude = range(template_velocity);
end


for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.VISp_clusters;
    
    exp_obj             = MismatchExperiment(data, config);
    
    for prot_i = protocol
        
        % get trials for this protocol
        trials = exp_obj.trials_of_type(prot_i);
        n_trials = length(trials);
        
        % get running traces for each trial to display
        [display_running, display_t] = exp_obj.running_around_mismatch_by_protocol(prot_i, display_window);
        
        % make sure it's consistent
        assert(n_trials == size(display_running, 2), ...
            'Size of running traces is not equal to number of trials');
        
        % preallocate cell array for spike rates of each cluster
        mm_spike_rate = cell(1, length(clusters));
        matched_spike_rate = cell(1, length(clusters));
        
        for cluster_i = 1 : length(clusters)
            
            % for this cluster, get the spike rate around the mismatch
            mm_spike_rate{cluster_i} = exp_obj.firing_around_mismatch_by_protocol(clusters(cluster_i), prot_i, display_window);
            
            % make sure it is the same number of trials
            assert(n_trials == size(display_running, 2), ...
                    'Size of spike matrix traces is not equal to number of trials');
            
            % get FiringRate object for later use in getting the firing
            % rate traces
            cluster_fr(cluster_i) = FiringRate(clusters(cluster_i).spike_times);
        end
        
        include_trial = false(1, n_trials);
        
        
        for trial_i = 1 : n_trials
            
            % determine whether the trial has changed velocity after
            % mismatch
            baseline_idx = display_t > window_1(1) & display_t < window_1(2);
            response_idx = display_t >= window_2(1) & display_t < window_2(2);
            
            m_before = mean(display_running(baseline_idx, trial_i));
            sd = std(display_running(baseline_idx, trial_i));
            
            m_after = mean(display_running(response_idx, trial_i));
            
            changed_down_idx = any(m_after < m_before - n_sds*sd);
            
            if ~changed_down_idx
                continue
            end
            
            % if we are not using the average trace, get the template from
            % the mismatch running profile
            if ~use_average_trace
                template_idx = display_t > 0 & display_t < 1;
                template_velocity = display_running(template_idx, trial_i);
                n_samples_to_compare = length(template_velocity);
            end
            
            % since, when we find a match we are looking back -1s
            % we need to start the search this many samples into
            % trial(x).velocity
            n_samples_before_match_to_show = find(display_t > 0, 1);
            n_samples_after_match_to_show = sum(display_t > 0);
            
            % onset of mismatch (trial sample index)
            mm_onset_t      = trials(trial_i).mismatch_onset_t();
            mm_onset_idx    = find(trials(trial_i).probe_t > mm_onset_t, 1, 'first');
            
            % we search in this space
            %   n_samples_before_match_to_show:n_samples_after_match
            search_idx = n_samples_before_match_to_show : mm_onset_idx - n_samples_after_match_to_show; %%%%%%
            
            % velocity to search
            search_velocity = trials(trial_i).velocity(search_idx);
            
            % number of samples we have to search
            n_samples_to_search = length(search_velocity) - n_samples_to_compare;
            
            % preallocate error and 
            err = nan(n_samples_to_search, 1);
            subsearch_amplitude = nan(n_samples_to_search, 1);
            
            for sample_i = 1 : n_samples_to_search
                
                subsearch_velocity = search_velocity(sample_i + (0:n_samples_to_compare-1));
                err(sample_i) = sum((subsearch_velocity - template_velocity).^2);
                subsearch_amplitude(sample_i) = range(subsearch_velocity);
            end
            
            [~, I] = min(err);
            
            % sort the errors
            [~, sort_idx] = sort(err, 'ascend');
            % sort the amplitudes of the 
            subsearch_amplitude_sorted = subsearch_amplitude(sort_idx);
            
            % find the first index at which amplitude is larger than
            % template, with minimum error
            I = find(subsearch_amplitude_sorted > template_amplitude, 1);
            
            % if there is no such index
            if isempty(I)
                % just take the minimum index of err
                [~, I] = min(err);
            else
                % take the minimum index of err
                I = sort_idx(I);
            end
            
%             matched_idx = I + (-first_sample_to_compare+1:n_samples_after_match-1);
%             matched_idx(matched_idx < 1) = 1;
            
            % I is the index in err with minimum distance between template
            % and search velocities
            %   the first point of 'err' is at 'n_samples_before_match_to_show'
            %   of trial(x).velocities
            %   so if I = 1, then we want to take from
            %       1 of the original velocity trace up through size of
            %       display_running
            matched_idx = I + (0:size(display_running, 1)-1);
            matched_idx(matched_idx < 1) = 1;
            
            store_running_mm = [store_running_mm, display_running(:, trial_i)];
            store_running_matched = [store_running_matched, trials(trial_i).velocity(matched_idx)];
            store_details = [store_details; probe_i, trial_i];
            
            
            
            n_samples = range(display_window) * trials(1).fs;
            common_t = display_window(1) + (0:n_samples-1)*(1/trials(trial_i).fs);
            matched_t = trials(trial_i).probe_t(I + n_samples_before_match_to_show);
            time_base = common_t + matched_t;
            
            for cluster_i = 1 : length(clusters)
                matched_spike_rate{cluster_i}(:, end+1) = cluster_fr(cluster_i).get_convolution(time_base);
            end
            
            include_trial(trial_i) = true;
        end
        
        
        for cluster_i = 1 : length(clusters)
            store_spikes_mm = [store_spikes_mm, mean(mm_spike_rate{cluster_i}(:, include_trial), 2)];
            store_spikes_matched = [store_spikes_matched, mean(matched_spike_rate{cluster_i}, 2)];
        end
    end
end




%% PLOT AVERAGES
yM = 60;

h_fig                   = figs.a4figure();

subplot(2, 2, 1)
hold on

plot(display_t, store_running_mm, 'color', [0.6, 0.6, 0.6]);
plot(display_t, mean(store_running_mm, 2), 'k');
ylim([0, yM]);
xlim([-1, 1]);
line([0, 0], [0, yM], 'color', 'k')
text(0, yM, sprintf('n = %i', size(store_running_mm, 2)), 'verticalalignment', 'top', 'horizontalalignment', 'left');
xlabel('Time from MM onset (s)')
ylabel('Running (cm/s)')
title('Mismatch trials');
box off


subplot(2, 2, 2)
hold on
plot(display_t, store_running_matched, 'color', [0.6, 0.6, 0.6]);
plot(display_t, mean(store_running_matched, 2), 'k');
ylim([0, yM]);
line([0, 0], [0, yM], 'color', 'k')
text(0, yM, sprintf('n = %i', size(store_running_matched, 2)), 'verticalalignment', 'top', 'horizontalalignment', 'left');
title('Matched trials');
box off


yL = [-4, 6];
bsl = display_t > -1 & display_t < 0;

subplot(2, 2, 3)
hold on

m_rm = bsxfun(@minus, store_spikes_mm, mean(store_spikes_mm(bsl, :), 1));
m = nanmean(m_rm, 2)';
s = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));
h = fill([display_t, display_t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], 'r');
set(h, 'facealpha', 0.6);
plot(display_t, m, 'r');
ylim(yL);
line([0, 0], yL, 'color', 'k');
ylabel('\Delta Hz');
title('Mismatch');
box off


subplot(2, 2, 4)
hold on
m_rm = bsxfun(@minus, store_spikes_matched, mean(store_spikes_matched(bsl, :), 1));
m = nanmean(m_rm, 2)';
s = nanstd(m_rm, [], 2)'/sqrt(sum(~isnan(m_rm(1, :))));
h = fill([display_t, display_t(end:-1:1)], [m-s, m(end:-1:1)+s(end:-1:1)], 'r');
set(h, 'facealpha', 0.6);
plot(display_t, m, 'r');
ylim(yL);
line([0, 0], yL, 'color', 'k');
title('Matched');
box off

FigureTitle(gcf, 'Average running speed around MM onset');
figs.save_fig(sprintf('averages_with_spikes_%s.pdf', fname));





