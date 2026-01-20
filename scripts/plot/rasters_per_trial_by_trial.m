close all;

% Modified version of rasters_per_trial.m that groups by trial instead of by cluster
%
%   This script creates raster plots where each page shows one trial,
%   with all clusters displayed together in the same raster plot.

%   The time axis is aligned to solenoid timing:
%       - t=0 is when solenoid goes DOWN
%       - Display range: from -5s (before solenoid down) to +10s after solenoid goes UP
%       - This means the displayed time range varies per trial
%
%   Specify options:
%
%       experiment_groups:      Will generate raster plots for all trials
%                               and all probe recordings 
%                               in the specified experiment group.
%
%       trial_group_labels:     Will generate a raster for all trials specified.
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%
%       common_fs:              sampling frequency to compute the
%                               firing rate convolutions for each cluster
%                               (e.g. 60Hz) 
%
%       save_figs:              true or false, whether to save the figures to pdf
%
%       overwrite:              true or false. If figure pdf's already exist,
%                               whether to overwrite 
%       
%       figure_dir:             cell array of strings specifying which
%                               directory to save pdf's.  

%%
experiment_groups       = {'ambient_light'};
trial_group_labels      = {'RT'};
common_fs               = 60;
save_figs               = true;
overwrite               = true;
figure_dir              = {'rasters_by_trial', 'ambient_light'};

%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.selected_clusters();
    
    % Load the spatial analysis cache to get cluster sorting order
    % This ensures rasters use the same order as spatial heatmaps
    spatial_cache_dir = fullfile(ctl.path_config.figure_dir, 'spatial_firing_rate', experiment_groups{:});
    cache_filename = sprintf('%s_spatial_analysis_cache.mat', probe_ids{ii});
    cache_filepath = fullfile(spatial_cache_dir, cache_filename);
    
    if exist(cache_filepath, 'file')
        fprintf('Loading spatial analysis cache to determine cluster sorting order...\n');
        cached = load(cache_filepath);
        
        % Compute the same sorting as in spatial_firing_rate_profile.m
        % Sort clusters by maximum peak positional firing rate for long trials
        long_rate_mat = cached.all_Q2_rate_smooth_by_group.long;  % Use median (Q2)
        bin_centers = cached.all_bin_centers_by_group.long;
        peak_positions = zeros(cached.n_clusters, 1);
        
        for c = 1:cached.n_clusters
            rate_profile = long_rate_mat(c, :);
            
            % Find the position of the maximum peak
            [~, max_idx] = max(rate_profile);
            peak_positions(c) = bin_centers(max_idx);
        end
        [~, sort_idx] = sort(peak_positions);
        
        % Reorder clusters to match the spatial heatmap order
        clusters = clusters(sort_idx);
        
        % Get significance information for each cluster
        cluster_significance_long = false(cached.n_clusters, 1);
        cluster_significance_short = false(cached.n_clusters, 1);
        
        if isfield(cached, 'spatial_tuning_stats') && ~isempty(cached.spatial_tuning_stats)
            stats = cached.spatial_tuning_stats;
            for c = 1:cached.n_clusters
                cluster_field = sprintf('cluster_%d', cached.cluster_ids(c));
                if isfield(stats, cluster_field)
                    cluster_stats = stats.(cluster_field);
                    if isfield(cluster_stats, 'long') && isfield(cluster_stats.long, 'isSpatiallyTuned')
                        cluster_significance_long(c) = cluster_stats.long.isSpatiallyTuned;
                    end
                    if isfield(cluster_stats, 'short') && isfield(cluster_stats.short, 'isSpatiallyTuned')
                        cluster_significance_short(c) = cluster_stats.short.isSpatiallyTuned;
                    end
                end
            end
        end
        
        % Reorder significance arrays to match sorted clusters
        cluster_significance_long = cluster_significance_long(sort_idx);
        cluster_significance_short = cluster_significance_short(sort_idx);
        
        fprintf('Clusters sorted by peak position in long trials (matching spatial heatmaps)\n');
        fprintf('  %d/%d clusters significant in long trials\n', sum(cluster_significance_long), cached.n_clusters);
        fprintf('  %d/%d clusters significant in short trials\n', sum(cluster_significance_short), cached.n_clusters);
    else
        warning('Spatial analysis cache not found at: %s\nClusters will not be sorted.', cache_filepath);
        % Default: all clusters non-significant
        cluster_significance_long = false(length(clusters), 1);
        cluster_significance_short = false(length(clusters), 1);
    end
    
    % Get all trials for the specified trial group
    all_trials = {};
    for kk = 1 : length(trial_group_labels)
        trials = data.get_trials_with_trial_group_label(trial_group_labels{kk});
        all_trials = [all_trials; trials];
    end
    
    % For each trial
    for trial_idx = 1 : length(all_trials)
        
        trial = all_trials{trial_idx};
        
        % Find solenoid timing for this trial
        % Solenoid goes down (high to low): diff(solenoid > 2.5) == -1
        % Solenoid goes up (low to high): diff(solenoid > 2.5) == 1
        sol_down_idx = find(diff(trial.solenoid > 2.5) == -1, 1);
        sol_up_idx = find(diff(trial.solenoid > 2.5) == 1, 1);
        
        if isempty(sol_down_idx) || isempty(sol_up_idx)
            warning('Trial %d: Could not find solenoid transitions, skipping...', trial.trial_id);
            continue;
        end
        
        % Get absolute times for solenoid events
        sol_down_time = trial.probe_t(sol_down_idx);
        sol_up_time = trial.probe_t(sol_up_idx);
        
        % Determine trial type based on maximum position
        max_position = max(trial.position);
        is_long_trial = max_position > 90;  % threshold at 90 cm
        
        % Set color based on trial type: blue for long, red for short
        if is_long_trial
            base_color = [0, 0, 1];  % blue
            trial_significance = cluster_significance_long;
        else
            base_color = [1, 0, 0];  % red
            trial_significance = cluster_significance_short;
        end
        
        % Get motion mask for this trial
        motion_mask = trial.motion_mask();
        
        % Use the entire trial time range, but make it relative to solenoid down
        % So t=0 is at solenoid down
        trial_start_time = trial.probe_t(1);
        trial_end_time = trial.probe_t(end);
        
        % Set limits relative to solenoid down (t=0)
        limits = [trial_start_time - sol_down_time, trial_end_time - sol_down_time];
        
        % Create one raster display object for this trial (with 1 section)
        r = RasterDisplayFigure(1);
        
        % Combine all clusters' spike data for this trial
        all_spike_times = {};
        all_spike_colors = {};  % Cell array of color arrays (one per cluster)
        all_spike_alphas = {};  % Cell array of alpha arrays (one per cluster)
        all_spike_rates = [];
        
        for cluster_idx = 1 : length(clusters)
            
            cluster = clusters(cluster_idx);
            
            % Create RasterData for this cluster and trial
            rd = RasterData(cluster, limits);
            rd.trigger_times = sol_down_time;  % Set t=0 at solenoid down
            
            % Get spike data for this cluster around solenoid down time
            spike_times = rd.spike_array();
            spike_rates = rd.spike_convolutions();
            
            % Determine which spikes occur during motion
            % Convert spike times (relative to sol_down_time) back to absolute probe time
            spike_times_absolute = spike_times{1} + sol_down_time;
            
            % Interpolate motion mask to spike times
            % motion_mask is sampled at trial.fs (usually 10000 Hz)
            spike_in_motion = interp1(trial.probe_t, double(motion_mask), spike_times_absolute, 'nearest', 0) > 0.5;
            
            % Create color array for each spike
            n_spikes = length(spike_times{1});
            spike_colors = zeros(n_spikes, 3);
            spike_alphas = zeros(n_spikes, 1);
            
            % Determine alpha for this cluster based on significance
            cluster_alpha = trial_significance(cluster_idx) * 0.5 + 0.5;  % 1.0 if significant, 0.5 if not
            
            for spike_i = 1:n_spikes
                if spike_in_motion(spike_i)
                    % During motion: use trial color (red/blue) with cluster alpha
                    spike_colors(spike_i, :) = base_color;
                else
                    % Not during motion: use black with cluster alpha
                    spike_colors(spike_i, :) = [0, 0, 0];
                end
                spike_alphas(spike_i) = cluster_alpha;
            end
            
            % Add to combined data
            all_spike_times = [all_spike_times; spike_times];
            all_spike_colors = [all_spike_colors; {spike_colors}];
            all_spike_alphas = [all_spike_alphas; {spike_alphas}];
            all_spike_rates = [all_spike_rates, spike_rates];
        end
        
        % Get velocity trace for this trial
        % Extract velocity from trial and interpolate to common time points
        trial_velocity = trial.velocity;
        trial_times = trial.probe_t;
        
        velocity_trace = nan(length(rd.common_t), 1);  % Single column vector
        
        % Interpolate entire trial velocity to common time points (relative to solenoid down)
        common_times_absolute = rd.common_t + sol_down_time;
        velocity_trace(:, 1) = interp1(trial_times, trial_velocity, common_times_absolute, 'linear', nan);
        
        % fill the section with the combined data
        sol_duration = sol_up_time - sol_down_time;
        if is_long_trial
            trial_type_str = 'Long';
        else
            trial_type_str = 'Short';
        end
        title_str = sprintf('Trial %i (%s trial) - All Clusters (Solenoid down at t=0, up at t=%.2fs)', trial.trial_id, trial_type_str, sol_duration);
        r.fill_data(1, 1, all_spike_times, velocity_trace, all_spike_rates, rd.common_t, title_str, all_spike_colors, all_spike_alphas, base_color);
        
        % set the limits and synchronize the axes of all the sections
        r.x_lim(limits);
        r.sync_sections();
        
        % give the page a title
        FigureTitle(r.h_fig, sprintf('%s, Trial %i (%s)', probe_ids{ii}, trial.trial_id, trial.trial_group_label));
        
        ctl.figs.save_fig_to_join(true, 400);
    end
    
    ctl.figs.join_figs(sprintf('%s.pdf', probe_ids{ii}), overwrite);
    ctl.figs.clear_figs();
end