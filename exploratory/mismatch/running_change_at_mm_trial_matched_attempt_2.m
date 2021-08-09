% For each mismatch trial, takes the running speed around the mismatch and
% looks for periods most similar to that running speed in the period before
%  and after the mismatch.

% Second attempt:
%   Allow the search for each mismatch period to be over ALL trials instead 
%   of the same trial (then take the # trials to compare best matches)

clearvars -except data

n_sds                        = 2;
window_1                     = [-0.15, 0.05];
window_2                     = [0.1, 0.3];
display_window               = [-1, 1];
protocol                     = 2;



%%
config                       = RC2AnalysisConfig();

figs                         = RC2Figures(config);
figs.save_on                 = true;
figs.set_figure_subdir('mismatch', 'running_change_at_mm', 'trial_matched');

recording_ids                = experiment_details('mismatch_nov20', 'protocol');

store_running_mm_change      = [];
store_running_mm_no_change   = [];
store_running_mm_all         = [];
store_running_matched        = [];

store_spikes_mm_change       = [];
store_spikes_mm_no_change    = [];
store_spikes_mm_all          = [];
store_spikes_matched         = [];

fname                        = '1s';

% average trace is a 1s long trace from -1 to 1s around mismatch onset
load( 'avg_to_compare_1s_incl_baseline.mat', 'avg_to_compare', 't_avg');

% template velocity
avg_to_take_start_t     = -0.5;
avg_to_take_end_t       = 1;
avg_to_take             = t_avg > avg_to_take_start_t & t_avg < avg_to_take_end_t;
template_velocity       = avg_to_compare(avg_to_take);

% the number of samples to compare on each window
n_samples_to_compare = length(template_velocity);
n_baseline_samples_in_avg_template = sum(t_avg > avg_to_take_start_t & t_avg < 0);


for rec_i = 1 : length(recording_ids)
    
    this_data           = get_data_for_recording_id(data, recording_ids{rec_i});
    clusters            = this_data.VISp_clusters([], 'any');
    
    if isempty(clusters)
        continue
    end
    
    exp_obj             = MismatchExperiment(this_data, config);
    
    for prot_i = protocol
        
        % get trials for this protocol
        trials = exp_obj.trials_of_type(prot_i);
        
        n_trials = length(trials);
        
        % get running traces for each trial to display
        [display_running, display_t] = exp_obj.running_around_mismatch_by_protocol(prot_i, display_window);
        
        n_samples = size(display_running, 1);
        common_t = display_window(1) + (0:n_samples-1)*(1/10e3);
        
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
        
        n_samples_before_mismatch_to_show = find(display_t > 0, 1); % -1
        n_samples_after_mismatch_to_show = sum(display_t > 0);
        
        % search these trials (with translation)
        trials_to_search = [exp_obj.trials_of_type(1), exp_obj.trials_of_type(2)];
        
        % store the data each loop
        search_window_velocities = {};
        search_window_probe_timebase = {};
        search_window_trial_sample_points = {};
        search_window_trial_number = [];
        
        % gather all running data for this mouse
        for trial_i = 1 : length(trials_to_search)
            
            this_trial = trials_to_search(trial_i);
            
            mismatch_onset_time     = this_trial.mismatch_onset_t();
            mismatch_onset_sample   = find(this_trial.probe_t > mismatch_onset_time, 1, 'first');
            
            analysis_window_mask    = this_trial.analysis_window();
            
            % we search between 
            % 1. 'n_samples_before_mismatch_to_show' (as we will have to take this many samples before any match for display)
            % 2. 'mismatch_onset_sample - n_samples_after_mismatch_to_show'
            %   (as we have to show this many samples after a match before the true mismatch onset)
            search_window_1_samples = (n_samples_before_mismatch_to_show + 1) : ...
                                      mismatch_onset_sample - n_samples_after_mismatch_to_show;
            
            % we also search between 
            % 1. 3s after mismatch onset
            % 2. last sample point of the analysis window - n_samples_after_mismatch_to_show
            search_window_2_samples = (mismatch_onset_sample + 3 * this_trial.fs) : ...
                                      find(analysis_window_mask, 1, 'last') - n_samples_after_mismatch_to_show;
            
            % velocities
            search_window_velocities{end+1} = this_trial.velocity(search_window_1_samples);
            search_window_velocities{end+1} = this_trial.velocity(search_window_2_samples);
            
            % corresponding probe times
            search_window_probe_timebase{end+1} = this_trial.probe_t(search_window_1_samples);
            search_window_probe_timebase{end+1} = this_trial.probe_t(search_window_2_samples);
            
            % corresponding sample points
            search_window_trial_sample_points{end+1} = search_window_1_samples;
            search_window_trial_sample_points{end+1} = search_window_2_samples;
            
            % trial # from which window comes
            search_window_trial_number(end+1) = trial_i;
            search_window_trial_number(end+1) = trial_i;
        end
        
        % only take search_windows with samples more than the number of
        % samples to compare
        invalid_search_window_mask = cellfun(@(x)(length(x) < n_samples_to_compare), search_window_velocities);
        search_window_velocities(invalid_search_window_mask) = [];
        search_window_probe_timebase(invalid_search_window_mask) = [];
        search_window_trial_sample_points(invalid_search_window_mask) = [];
        search_window_trial_number(invalid_search_window_mask) = [];
        
        % number of search windows
        n_search_windows = length(search_window_velocities);
        
        error_at_sample_point_of_search_window = cell(n_search_windows, 1);
        
        % do the search
        for search_i = 1 : n_search_windows
            
            %
            n_samples_in_search_window = length(search_window_velocities{search_i});
            
            % number of samples we have to search
            n_samples_to_search = n_samples_in_search_window - n_samples_to_compare + 1;
            
            % preallocate error
            error_at_sample_point_of_search_window{search_i} = nan(n_samples_to_search, 1);
            
            fprintf('Doing the search... ');
            for sample_i = 1 : n_samples_to_search
                
                if mod(sample_i, 2000) == 1
                    str = sprintf('%i/%i', sample_i, n_samples_to_search);
                    fprintf('%s\n', str);
                end
                
                velocity_to_compare = search_window_velocities{search_i}(sample_i + (0:n_samples_to_compare-1));
                
                error_at_sample_point_of_search_window{search_i}(sample_i) = sum((velocity_to_compare - template_velocity).^2);
            end
            fprintf('done\n');
        end
        
        
        trial_changed_velocity = false(1, n_trials);
        
        for trial_i = 1 : n_trials
            
            % was there a change in velocity for this trial?
            trial_changed_velocity(trial_i) = exp_obj.trial_changed_velocity(trials(trial_i).id, window_1, window_2, n_sds);
            
            % skip trial if there was no change in velocity
            if trial_changed_velocity(trial_i)
                store_running_mm_change = [store_running_mm_change, display_running(:, trial_i)];
            else
                store_running_mm_no_change = [store_running_mm_no_change, display_running(:, trial_i)];
            end
            store_running_mm_all = [store_running_mm_all, display_running(:, trial_i)];
        end
        
        n_trials_to_match = sum(trial_changed_velocity);
        
        min_error_in_search_window          = nan(1, n_search_windows);
        search_window_location_of_min_error = nan(1, n_search_windows);
        
        for search_i = 1 : n_search_windows
            
            [min_error_in_search_window(search_i), search_window_location_of_min_error(search_i)] = ...
                min(error_at_sample_point_of_search_window{search_i});
        end
        
        % sort these errors
        [~, search_window_idx_sorted_by_error] = sort(min_error_in_search_window, 'ascend');
        
        % which search windows to take
        search_window_idx_sorted_by_error = search_window_idx_sorted_by_error(1:n_trials_to_match);
        
        
        for search_i = 1 : n_trials_to_match
            
            fprintf(' Search trial %i/%i\n', search_i, n_trials_to_match);
            
            % 
            this_search_window_idx                  = search_window_idx_sorted_by_error(search_i);
            
            
            this_search_window_velocity             = search_window_velocities{this_search_window_idx};
            this_search_window_probe_timebase       = search_window_probe_timebase{this_search_window_idx};
            this_search_window_trial_sample_points  = search_window_trial_sample_points{this_search_window_idx};
            this_trial                              = trials_to_search(search_window_trial_number(this_search_window_idx));
            
            sample_to_look_in_this_search_window    = search_window_location_of_min_error(this_search_window_idx);
            
            
            
            trial_sample_point_starting_at_match    = ...
                this_search_window_trial_sample_points(sample_to_look_in_this_search_window);
            
            trial_sample_points_to_take_for_velocity = ...
                trial_sample_point_starting_at_match - n_samples_before_mismatch_to_show + (0 : length(display_t)-1) + ...
                n_baseline_samples_in_avg_template;
            
            vel = this_trial.velocity(trial_sample_points_to_take_for_velocity);
            
            store_running_matched = [store_running_matched, vel(:)];
            
            time_base = common_t + this_search_window_probe_timebase(sample_to_look_in_this_search_window + n_baseline_samples_in_avg_template);
            
            for cluster_i = 1 : length(clusters)
                
                fprintf('  Cluster %i/%i\n', cluster_i, length(clusters));
                matched_spike_rate{cluster_i}(:, end+1) = cluster_fr(cluster_i).get_convolution(time_base);
            end
        end
       
        
        for cluster_i = 1 : length(clusters)
            
            store_spikes_mm_change      = [store_spikes_mm_change, mean(mm_spike_rate{cluster_i}(:, trial_changed_velocity), 2)];
            store_spikes_mm_no_change   = [store_spikes_mm_no_change, mean(mm_spike_rate{cluster_i}(:, ~trial_changed_velocity), 2)];
            store_spikes_mm_all         = [store_spikes_mm_all, mean(mm_spike_rate{cluster_i}, 2)];
            store_spikes_matched        = [store_spikes_matched, mean(matched_spike_rate{cluster_i}, 2)];
        end
    end
end




%% PLOT AVERAGES

h_fig                   = figs.a4figure();

h_ax                    = subplot(2, 4, 1);
running_trace_plot(h_ax, display_t, store_running_mm_all)
xlabel(h_ax, 'Time from MM onset (s)')
ylabel(h_ax, 'Running (cm/s)')
title(h_ax, {'Mismatch trials', '(all)'});

h_ax                    = subplot(2, 4, 2);
running_trace_plot(h_ax, display_t, store_running_mm_change)
title(h_ax, {'Mismatch trials', '(change)'});

h_ax                    = subplot(2, 4, 3);
running_trace_plot(h_ax, display_t, store_running_mm_no_change)
title(h_ax, {'Mismatch trials', '(no change)'});

h_ax                     = subplot(2, 4, 4);
running_trace_plot(h_ax, display_t, store_running_matched)
title(h_ax, 'Matched trials');


h_ax                     = subplot(2, 4, 5);
spiking_trace_plot(h_ax, display_t, store_spikes_mm_all)
ylabel(h_ax, '\Delta Hz');

h_ax                     = subplot(2, 4, 6);
spiking_trace_plot(h_ax, display_t, store_spikes_mm_change)

h_ax                     = subplot(2, 4, 7);
spiking_trace_plot(h_ax, display_t, store_spikes_mm_no_change)

h_ax                     = subplot(2, 4, 8);
spiking_trace_plot(h_ax, display_t, store_spikes_matched)

FigureTitle(gcf, 'Average running speed around MM onset');

date_time = datestr(now, 'yyyymmdd_HHMMSS');
full_fname_root = sprintf('averages_with_spikes_%s', date_time);

save(fullfile(figs.curr_dir, full_fname_root), ...
     'store_running_matched', ...
     'store_running_mm_change', ...
     'store_spikes_mm_change', ...
     'store_spikes_matched', ...
     'display_t');
 
figs.save_fig([full_fname_root, '.pdf']);



%% print p-value
baseline_idx = display_t >= -0.4 & display_t < 0;
response_idx = display_t >= 0 & display_t < 0.4;

baseline_mm = mean(store_spikes_mm_change(baseline_idx, :), 1);
response_mm = mean(store_spikes_mm_change(response_idx, :), 1);

baseline_matched = mean(store_spikes_matched(baseline_idx, :), 1);
response_matched = mean(store_spikes_matched(response_idx, :), 1);

p_mm = ranksum(baseline_mm, response_mm);
p_matched = ranksum(baseline_matched, response_matched);

fprintf('p-value spiking baseline vs. response period, mismatch trials: %.2e\n', p_mm);
fprintf('p-value spiking baseline vs. response period, matched trials: %.2f\n', p_matched);


%% print p-value
baseline_idx = display_t >= -0.4 & display_t < 0;
response_idx = display_t >= 0 & display_t < 0.4;

baseline_mm = mean(store_running_mm_change(baseline_idx, :), 1);
response_mm = mean(store_running_mm_change(response_idx, :), 1);

baseline_matched = mean(store_running_matched(baseline_idx, :), 1);
response_matched = mean(store_running_matched(response_idx, :), 1);

p_mm = ranksum(baseline_mm, baseline_matched);
p_matched = ranksum(response_mm, response_matched);

fprintf('p-value running baseline, median, mismatch vs matched: %.2f\n', p_mm);
fprintf('p-value running response, median, mismatch vs matched: %.2f\n', p_matched);


[~, p_mm] = vartest2(baseline_mm, baseline_matched);
[~, p_matched] = vartest2(response_mm, response_matched);

fprintf('p-value running baseline, variance, mismatch vs matched: %.2e\n', p_mm);
fprintf('p-value running response, variance, mismatch vs matched: %.2e\n', p_matched);




