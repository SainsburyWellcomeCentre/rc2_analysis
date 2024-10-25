% For each mismatch trial, takes the running speed around the mismatch and
% looks for periods most similar to that running speed in the period before
%  and after the mismatch.

% Second attempt:
%   Allow the search for each mismatch period to be over ALL trials instead 
%   of the same trial (then take the # trials to compare best matches)


%%
experiment_groups       = {'mismatch_nov20', 'mismatch_jul21'};
trial_group_label       = 'RVT_gain_up';
n_sds                   = 2;
window_1                = [-0.15, 0.05];
window_2                = [0.1, 0.3];
template_limits         = [-0.5, 1];
limits                  = [-1, 1];
fs                      = 10e3;

save_figs               = false;
overwrite               = true;
figure_dir              = {'matched_mismatch'};


%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

ctl.setup_figures(figure_dir, save_figs);

traces = cell(1, length(probe_ids));
mismatch_times = cell(1, length(probe_ids));
running_decreases = cell(1, length(probe_ids));
n_change_trials = nan(1, length(probe_ids));

decrease_raster_data = {};
no_change_raster_data = {};

% 1. get average running trace at mismatch onset
for ii = 1 : length(probe_ids)
    
    data                    = ctl.load_formatted_data(probe_ids{ii});
    clusters                = data.VISp_clusters();
    common_t                = data.timebase(limits, fs);
    traces{ii}              = data.get_traces_around_mismatch_onset(trial_group_label, limits, fs); %#ok<*AGROW>
    mismatch_times{ii}      = data.get_mismatch_onset_times(trial_group_label);
    
    % determine whether each trace shows decrease in running
    running_decreases{ii}   = does_mismatch_running_decrease(common_t, traces{ii}, window_1, window_2, n_sds);
    n_change_trials(ii)     = sum(running_decreases{ii});
    
    % for each cluster get the raster data
    for jj = 1 : length(clusters)
        
        decrease_raster_data{end+1} = RasterData(clusters(jj), limits);
        decrease_raster_data{end}.trigger_times = mismatch_times{ii}(running_decreases{ii});
        decrease_raster_data{end}.fs = fs;
        
        no_change_raster_data{end+1} = RasterData(clusters(jj), limits);
        no_change_raster_data{end}.trigger_times = mismatch_times{ii}(~running_decreases{ii});
        no_change_raster_data{end}.fs = fs;
    end
end
%%
% take the average of all 'change' trials
all_traces = [traces{:}];
template_idx = common_t >= template_limits(1) & common_t <= template_limits(2);
template_trace = mean(all_traces(template_idx, logical([running_decreases{:}])), 2);

% the number of samples to skip at the start of each window
n_samples_to_skip = sum(common_t < template_limits(1));
n_samples_to_zero = sum(common_t >= template_limits(1) & common_t < 0);
after_mismatch_to_skip = 3;  % sec

matched_traces      = {};
matched_raster_data = {};

%%
% go through all traces and find matches
for ii = 1 : length(probe_ids)
    
    data                    = ctl.load_formatted_data(probe_ids{ii});
    trials                  = data.get_trials_with_trial_group_label({'RVT_gain_up', 'RVT_gain_down'});
    clusters                = data.VISp_clusters();
    
    velocities_to_search    = {};
    probe_t                 = {}; 
    store_trials            = {};
    
    % go through trials and collect velocities to search (excl. mismatch)
    for jj = 1 : length(trials)
        
        onset_idx = trials{jj}.mismatch_onset_sample();
        
        % we search between
        % 1. 'n_samples_before_mismatch_to_show' (as we will have to take this many samples before any match for display)
        % 2. 'mismatch_onset_sample - n_samples_after_mismatch_to_show'
        %   (as we have to show this many samples after a match before the true mismatch onset)
        idx_1 = (n_samples_to_skip + 1) : onset_idx;
        idx_2 = (onset_idx + after_mismatch_to_skip * trials{jj}.fs):find(trials{jj}.analysis_window, 1, 'last');
        
        velocities_to_search{end+1} = trials{jj}.velocity(idx_1); %#ok<*SAGROW>
        velocities_to_search{end+1} = trials{jj}.velocity(idx_2);
        
        probe_t{end+1} = trials{jj}.probe_t(idx_1);
        probe_t{end+1} = trials{jj}.probe_t(idx_2);
        
        store_trials{end+1} = trials{jj};
        store_trials{end+1} = trials{jj};
    end
    
    to_remove = cellfun(@length, velocities_to_search) < length(template_trace);
    velocities_to_search(to_remove) = [];
    probe_t(to_remove) = [];
    store_trials(to_remove) = [];
    
    best_error = nan(length(velocities_to_search), 1);
    best_sample = nan(length(velocities_to_search), 1);
    best_time = nan(length(velocities_to_search), 1);
    
    % for each available window, search for matches
    for jj = 1 : length(velocities_to_search)
        
        % display progress
        fprintf('%i/%i, %i/%i mice\n', jj, length(velocities_to_search), ii, length(probe_ids));
        
        % number of samples which we can align to the template trace start
        n_samples_to_search = length(velocities_to_search{jj}) - length(template_trace) + 1;
        
        % preallocate error between this trace and the template for each of those samples
        error = nan(n_samples_to_search, 1);
        
        for kk = 1 : n_samples_to_search
            
            % select the bit of the window and do the comparison
            cmp = velocities_to_search{jj}(kk + (0:length(template_trace)-1));
            error(kk) = sum((cmp - template_trace).^2);
        end
        
        % for this window, find the best error and location of the best error
        [best_error(jj), idx] = min(error);
        
        % location of the sample which is aligned to zero time in the
        % template trace
        best_sample = idx + n_samples_to_zero;
        
        % the time of that sample point
        best_time(jj) = probe_t{jj}(best_sample);
    end
    
    % sort the windows from best to worst. each window gets a 'best' score,
    % and then we take the best of the best
    [~, glob_idx] = sort(best_error, 'ascend');
    
    % the times of those best moments
    matched_times{ii} = best_time(glob_idx(1:n_change_trials(ii)));
    
    % the corresponding trial objects for those moments
    matched_trials = store_trials(glob_idx(1:n_change_trials(ii)));
    
    % the running traces around those moments
    matched_traces{ii} = data.get_traces_around_times(matched_trials, matched_times{ii}, limits, fs);
    
    % for each cluster get the raster data
    for jj = 1 : length(clusters)
        
        matched_raster_data{end+1} = RasterData(clusters(jj), limits);
        matched_raster_data{end}.trigger_times = matched_times{ii};
        matched_raster_data{end}.fs = fs;
    end
end


%% save the data for the future
decrease_traces = cellfun(@(x, y)(x(:, y)), traces, running_decreases, 'uniformoutput', false);
no_change_traces = cellfun(@(x, y)(x(:, ~y)), traces, running_decreases, 'uniformoutput', false);
matched_decrease_traces = matched_traces;

decrease_traces_fr = cellfun(@(x)(x.spike_convolution_avg), decrease_raster_data, 'UniformOutput', false);
no_change_traces_fr = cellfun(@(x)(x.spike_convolution_avg), no_change_raster_data, 'UniformOutput', false);
matched_decrease_traces_fr = cellfun(@(x)(x.spike_convolution_avg), matched_raster_data, 'UniformOutput', false);

% library git info
lib_git = ctl.save.git.info;

save('running_around_mismatch_matched_trials', 'lib_git', ...
                                               'common_t', ...
                                               'decrease_traces', ...
                                               'no_change_traces', ...
                                               'matched_decrease_traces', ...
                                               'decrease_traces_fr', ...
                                               'no_change_traces_fr', ...
                                               'matched_decrease_traces_fr');
