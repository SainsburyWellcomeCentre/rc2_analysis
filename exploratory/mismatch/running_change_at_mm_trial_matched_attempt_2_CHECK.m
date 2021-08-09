% For each mismatch trial, takes the running speed around the mismatch and
% looks for periods most similar to that running speed in the period before
%  and after the mismatch.

% Second attempt:
%   Allow the search for each mismatch period to be over ALL trials instead 
%   of the same trial (then take the # trials to compare best matches)

% This script checks the main script by allowing the search to include the
% mismatch period (error should be 0) and this should recover the exact
% shape.

clearvars -except data

n_sds                   = 2;
window_1                = [-0.15, 0.05];
window_2                = [0.1, 0.3];
display_window          = [-1, 1];
protocol                = 2;



%%
config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('mismatch', 'running_change_at_mm', 'trial_matched');

recording_ids            = experiment_details('mismatch_nov20', 'protocol');

store_running_mm        = [];
store_running_matched   = [];
store_spikes_mm         = [];
store_spikes_matched    = [];

fname                   = '1s';


% average trace is a 1s long trace from -1 to 1s around mismatch onset
load( 'avg_to_compare_1s_incl_baseline.mat', 'avg_to_compare', 't_avg');

% template velocity
avg_to_take_start_t     = -0.5;
avg_to_take_end_t       = 0.5;
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
        trials_in_protocol = exp_obj.trials_of_type(prot_i);
        n_trials_in_protocol = length(trials_in_protocol);
        
        
        % get running traces for each trial to display
        [running_traces_around_mismatch_for_protocol, display_t] = exp_obj.running_around_mismatch_by_protocol(prot_i, display_window);
        
        n_samples_to_display = length(display_t);
        timebase_common_to_all_display_traces = display_window(1) + (0:n_samples_to_display-1)*(1/trials_in_protocol(1).fs);
        
        % preallocate cell array for spike rates of each cluster
        mm_spike_rate = cell(1, length(clusters));
        matched_spike_rate = cell(1, length(clusters));
        
        for cluster_i = 1 : length(clusters)
            
            % for this cluster, get the spike rate around the mismatch
            mm_spike_rate{cluster_i} = exp_obj.firing_around_mismatch_by_protocol(clusters(cluster_i), prot_i, display_window);
            
            % make sure it is the same number of trials_in_protocol
            assert(n_trials_in_protocol == size(running_traces_around_mismatch_for_protocol, 2), ...
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
        search_window_trial_number_id = [];
        check_template_velocity = {};
        
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
                                      mismatch_onset_sample + n_samples_after_mismatch_to_show + 1000;
            
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
            
            search_window_trial_number_id(end+1) = trials_to_search(trial_i).id;
            search_window_trial_number_id(end+1) = trials_to_search(trial_i).id;
            
            template_idx = mismatch_onset_sample - n_baseline_samples_in_avg_template + (0:n_samples_to_compare-1);
            check_template_velocity{end+1} = this_trial.velocity(template_idx);
            check_template_velocity{end+1} = this_trial.velocity(template_idx);
        end
        
        include_trial = false(1, n_trials_in_protocol);
        
        for trial_i = 1 : n_trials_in_protocol
            
            % was there a change in velocity for this trial?
            include_trial(trial_i) = exp_obj.trial_changed_velocity(trials_in_protocol(trial_i).id, window_1, window_2, n_sds);
            
            % skip trial if there was no change in velocity
            if include_trial(trial_i)
                store_running_mm = [store_running_mm, running_traces_around_mismatch_for_protocol(:, trial_i)];
            end
        end
        
        % only take search_windows with samples more than the number of
        % samples to compare
        invalid_search_window_mask = cellfun(@(x)(length(x) < n_samples_to_compare), search_window_velocities);
        search_window_velocities(invalid_search_window_mask) = [];
        search_window_probe_timebase(invalid_search_window_mask) = [];
        search_window_trial_sample_points(invalid_search_window_mask) = [];
        search_window_trial_number(invalid_search_window_mask) = [];
        search_window_trial_number_id(invalid_search_window_mask) = [];
        check_template_velocity(invalid_search_window_mask) = [];
        
        % number of search windows
        n_search_windows = length(search_window_velocities);
        
        error_at_sample_point_of_search_window = cell(n_search_windows, 1);
        
        % find corresponding window
        for search_i = 1 : n_search_windows
            
            %
            n_samples_in_search_window = length(search_window_velocities{search_i});
            
            % number of samples we have to search
            n_samples_to_search = n_samples_in_search_window - n_samples_to_compare + 1;
            
            % preallocate error
            error_at_sample_point_of_search_window{search_i} = nan(n_samples_to_search, 1);
        end
        
        
        % do the search
        for trial_i = 1 : n_trials_in_protocol
            
            if include_trial(trial_i)
                
                trial_id = trials_in_protocol(trial_i).id;
                
                % find search windows for this trial
                search_window_idx = find(search_window_trial_number_id == trial_id);
                
                % find corresponding window
                for search_i = search_window_idx
                    
                    %
                    n_samples_in_search_window = length(search_window_velocities{search_i});
                    
                    % number of samples we have to search
                    n_samples_to_search = n_samples_in_search_window - n_samples_to_compare + 1;
                    
                    % preallocate error
                    error_at_sample_point_of_search_window{search_i} = nan(n_samples_to_search, 1);
                    
                    if search_window_trial_number_id(search_i) ~= trial_id
                        continue
                    end
                    
                    fprintf('Doing the search... ');
                    for sample_i = 1 : n_samples_to_search
                        
                        if mod(sample_i, 2000) == 1
                            str = sprintf('%i/%i', sample_i, n_samples_to_search);
                            fprintf('%s\n', str);
                        end
                        
                        velocity_to_compare = search_window_velocities{search_i}(sample_i + (0:n_samples_to_compare-1));
                        
                        error_at_sample_point_of_search_window{search_i}(sample_i) = sum((velocity_to_compare - check_template_velocity{search_i}).^2);
                    end
                    fprintf('done\n');
                end
                
            end
        end
        
        
        
        
        
        
        
        n_trials_to_match = sum(include_trial);
        
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
            
            time_base = timebase_common_to_all_display_traces + this_search_window_probe_timebase(sample_to_look_in_this_search_window + n_baseline_samples_in_avg_template);
            
            for cluster_i = 1 : length(clusters)
                
                fprintf('  Cluster %i/%i\n', cluster_i, length(clusters));
                matched_spike_rate{cluster_i}(:, end+1) = cluster_fr(cluster_i).get_convolution(time_base);
            end
        end
       
        
        for cluster_i = 1 : length(clusters)
            
            store_spikes_mm         = [store_spikes_mm, mean(mm_spike_rate{cluster_i}(:, include_trial), 2)];
            store_spikes_matched    = [store_spikes_matched, mean(matched_spike_rate{cluster_i}, 2)];
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


date_time = datestr(now, 'yyyymmdd_HHMMSS');
full_fname_root = sprintf('averages_with_spikes_%s_%s', fname, date_time);

save(full_fname_root, ...
     'store_running_matched', ...
     'store_running_mm', ...
     'store_spikes_mm', ...
     'store_spikes_matched', ...
     'display_t');
 
figs.save_fig([full_fname_root, '.pdf']);
