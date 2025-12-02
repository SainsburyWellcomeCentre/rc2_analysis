% Plot spatial firing rate profile for each cluster using position as x-axis
% For experiment group 'training_running', bin position into 2-cm bins,
% normalize spike counts by occupancy, and smooth with a Gaussian kernel (8-cm s.d.)
%
% SMOOTHING APPROACH:
%   Spike counts and occupancy are smoothed SEPARATELY before computing rates.
%   This is mathematically cleaner than smoothing rates directly:
%       rate_smooth = smooth(counts) / smooth(occupancy)
%   rather than:
%       rate_smooth = smooth(counts / occupancy)
%
%   Specify options:
%
%       experiment_groups:      Will generate plots for all trials
%                               and all probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21',
%                                   'training_running'
%                               Should be a cell array of strings with each
%                               entry an experiment group
%
%       save_figs:              true or false, whether to save the figures to pdf
%
%       overwrite:              true or false. If figure pdf's already exist,
%                               whether to overwrite 
%       
%       figure_dir:             cell array of strings specifying which
%                               directory to save pdf's. The directory will
%                               be relative to the directory specified by
%                               path_config.figure_dir (in
%                               `path_config.m`), so that {'one', 'two',
%                               'three'} will save .pdfs to:
%                               <path_config.figure_dir>\one\two\three\      
%
% If `save_figs` is true, one pdf will be created for each probe recording,
% and contain spatial firing rate profiles for all clusters in that probe.
%
% Requires functions from lib/spatial_analysis/:
%   - compute_trial_firing_rate.m
%   - smooth_spatial_rate.m
%   - compute_quartiles_and_smooth.m
%   - collect_and_group_trials.m

%%
% Configuration
experiment_groups       = {'ambient_light'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'spatial_firing_rate', 'ambient_light'};
plot_single_cluster_fig = true;
plot_heatmp_cluster_fig = true;
re_run_analysis         = false;  % Set to true to recompute all metrics, false to load cached data


% Analysis parameters
bin_size_cm = 2;
gauss_sigma_cm = 8; % standard deviation for smoothing (cm)
position_threshold_cm = 90; % threshold for long/short trial classification


% Initialize controller and get probe IDs
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));


for pid = 1:length(probe_ids)
    fprintf('\nProcessing probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    
    % Define cache file path
    cache_filename = sprintf('%s_spatial_analysis_cache.mat', probe_ids{pid});
    cache_filepath = fullfile(ctl.figs.curr_dir, cache_filename);
    
    % Try to load cached data if not re-running analysis
    if ~re_run_analysis && exist(cache_filepath, 'file')
        fprintf('  Loading cached analysis data from: %s\n', cache_filename);
        load(cache_filepath, 'all_rate_smooth_by_group', 'all_rate_unsmooth_by_group', ...
             'all_bin_centers_by_group', 'all_avg_velocity_by_group', ...
             'all_Q1_rate_smooth_by_group', 'all_Q2_rate_smooth_by_group', ...
             'all_Q3_rate_smooth_by_group', 'all_occ_by_group', ...
             'all_rate_per_trial_by_group', 'all_rate_per_trial_smooth_by_group', ...
             'all_spike_positions_by_group', 'all_global_trial_indices_by_group', ...
             'all_rate_pooled_by_group', ...
             'cluster_ids', 'n_clusters', 'group_names', 'group_labels');
        fprintf('  Loaded cached data for %d clusters\n', n_clusters);
        
        % Still need to load data object for clusters info
        data = ctl.load_formatted_data(probe_ids{pid});
        clusters = data.selected_clusters();
        sessions = data.motion_sessions();
    else
        if ~re_run_analysis
            fprintf('  No cache found, computing analysis...\n');
        else
            fprintf('  Re-running analysis (re_run_analysis=true)...\n');
        end
        
        data = ctl.load_formatted_data(probe_ids{pid});
        clusters = data.selected_clusters();
        sessions = data.motion_sessions();
        
        fprintf('  Found %d clusters and %d sessions\n', length(clusters), length(sessions));
        
        % --- Collect and group trials using helper function ---
        [trial_groups, group_names, group_labels] = collect_and_group_trials(sessions, position_threshold_cm);
        fprintf('  Total long trials: %d\n', numel(trial_groups.long));
        fprintf('  Total short trials: %d\n', numel(trial_groups.short));

    % --- Initialize storage structures ---
    all_rate_smooth_by_group = struct();
    all_rate_unsmooth_by_group = struct();
    all_bin_centers_by_group = struct();
    all_avg_velocity_by_group = struct();
    all_Q1_rate_unsmooth_by_group = struct();
    all_Q1_rate_smooth_by_group = struct();
    all_Q2_rate_unsmooth_by_group = struct();
    all_Q2_rate_smooth_by_group = struct();
    all_Q3_rate_unsmooth_by_group = struct();
    all_Q3_rate_smooth_by_group = struct();
    all_occ_by_group = struct();
    
    % Storage for per-trial data (for plotting individual traces)
    all_rate_per_trial_by_group = struct();
    all_rate_per_trial_smooth_by_group = struct();
    all_spike_positions_by_group = struct();  % For position-based raster plots
    all_global_trial_indices_by_group = struct();  % For preserving trial order
    
    % Storage for pooled trial rates (all trials combined before smoothing)
    all_rate_pooled_by_group = struct();
    
    cluster_ids = [clusters.id];
    n_clusters = length(clusters);

    % --- Compute spatial firing rate profiles for all clusters ---
    fprintf('  Computing spatial firing rate profiles for %d clusters...\n', n_clusters);
    
    for c = 1:n_clusters
        cluster = clusters(c);
        
        for g = 1:2
            group = group_names{g};
            trials_struct = trial_groups.(group);
            
            % Set bin edges explicitly for each group
            if strcmp(group, 'long')
                edges = 0:bin_size_cm:120;
            else
                edges = 0:bin_size_cm:60;
            end
            bin_centers = edges(1:end-1) + bin_size_cm/2;
            n_bins = length(bin_centers);
            
            % Adjust bin centers for plotting (short trials: 60-120 cm)
            if strcmp(group, 'short')
                plot_bin_centers = bin_centers + 60;
            else
                plot_bin_centers = bin_centers;
            end
            
            % --- Compute trial-by-trial firing rates using helper function ---
            n_trials = length(trials_struct);
            rate_per_trial = nan(n_trials, n_bins);
            occ_per_trial = nan(n_trials, n_bins);
            vel_per_trial = nan(n_trials, n_bins);
            count_per_trial = nan(n_trials, n_bins);  % Store spike counts
            spike_positions_per_trial = cell(n_trials, 1);  % For raster plots
            global_trial_indices = zeros(n_trials, 1);  % Store global trial order
            
            for k = 1:n_trials
                trial = trials_struct(k).trial;
                global_trial_indices(k) = trials_struct(k).global_trial_idx;  % Store global index
                [rate, occ, vel, counts] = compute_trial_firing_rate(trial, cluster, edges);
                rate_per_trial(k, :) = rate;
                occ_per_trial(k, :) = occ;
                vel_per_trial(k, :) = vel;
                count_per_trial(k, :) = counts;
                
                % Get spike positions for position-based raster plot (motion only)
                motion_mask = trial.motion_mask();
                pos = trial.position(motion_mask);
                tvec = trial.probe_t(motion_mask);
                if ~isempty(pos) && ~isempty(tvec)
                    st = cluster.spike_times;
                    spike_mask = st >= tvec(1) & st <= tvec(end);
                    st_in_trial = st(spike_mask);
                    if ~isempty(st_in_trial)
                        spike_pos = interp1(tvec, pos, st_in_trial, 'linear', 'extrap');
                        % Adjust short trial positions to 60-120 range
                        if strcmp(group, 'short')
                            spike_pos = spike_pos + 60;
                        end
                        spike_positions_per_trial{k} = spike_pos(:)';
                    else
                        spike_positions_per_trial{k} = [];
                    end
                else
                    spike_positions_per_trial{k} = [];
                end
            end
            
            % Smooth spike counts and occupancy separately, then compute rates
            % This is more mathematically sound than smoothing rates directly
            count_per_trial_smooth = nan(n_trials, n_bins);
            occ_per_trial_smooth = nan(n_trials, n_bins);
            rate_per_trial_smooth = nan(n_trials, n_bins);
            
            for k = 1:n_trials
                % Smooth counts and occupancy
                count_per_trial_smooth(k, :) = smooth_spatial_rate(count_per_trial(k, :), bin_size_cm, gauss_sigma_cm);
                occ_per_trial_smooth(k, :) = smooth_spatial_rate(occ_per_trial(k, :), bin_size_cm, gauss_sigma_cm);
                
                % Compute smoothed rate from smoothed counts and occupancy
                rate_per_trial_smooth(k, :) = count_per_trial_smooth(k, :) ./ occ_per_trial_smooth(k, :);
                rate_per_trial_smooth(k, isnan(rate_per_trial_smooth(k, :)) | isinf(rate_per_trial_smooth(k, :))) = 0;
            end
            
            % Compute quartiles from SMOOTHED trial data
            Q1_from_smoothed = prctile(rate_per_trial_smooth, 25, 1);
            Q2_from_smoothed = prctile(rate_per_trial_smooth, 50, 1);  % Median
            Q3_from_smoothed = prctile(rate_per_trial_smooth, 75, 1);
            
            % Compute mean across trials using smoothed counts/occupancy approach
            mean_count_smooth = nanmean(count_per_trial_smooth, 1);
            mean_occ_smooth = nanmean(occ_per_trial_smooth, 1);
            mean_rate_smooth = mean_count_smooth ./ mean_occ_smooth;
            mean_rate_smooth(isnan(mean_rate_smooth) | isinf(mean_rate_smooth)) = 0;
            
            % Also compute unsmoothed mean for comparison
            mean_rate_unsmooth = nanmean(rate_per_trial, 1);
            
            % Compute POOLED trial rate: sum all counts/occupancy, then smooth
            % This is for comparison with the _identification script approach
            pooled_count = nansum(count_per_trial, 1);
            pooled_occ = nansum(occ_per_trial, 1);
            pooled_count_smooth = smooth_spatial_rate(pooled_count, bin_size_cm, gauss_sigma_cm);
            pooled_occ_smooth = smooth_spatial_rate(pooled_occ, bin_size_cm, gauss_sigma_cm);
            rate_pooled = pooled_count_smooth ./ pooled_occ_smooth;
            rate_pooled(isnan(rate_pooled) | isinf(rate_pooled)) = 0;
            
            % --- Print statistics ---
            fprintf('\n    Cluster %d, %s trials (n=%d):\n', cluster.id, group, n_trials);
            fprintf('      Smoothed mean rate: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(mean_rate_smooth), min(mean_rate_smooth), max(mean_rate_smooth));
            fprintf('      Smoothed Q2 (median): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(Q2_from_smoothed), min(Q2_from_smoothed), max(Q2_from_smoothed));
            fprintf('      Pooled rate: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(rate_pooled), min(rate_pooled), max(rate_pooled));
            
            % Compute mean occupancy and velocity across trials
            mean_occupancy = nanmean(occ_per_trial, 1);
            mean_velocity = nanmean(vel_per_trial, 1);
            
            % --- Initialize storage if needed ---
            if ~isfield(all_rate_smooth_by_group, group)
                all_rate_smooth_by_group.(group) = nan(n_clusters, n_bins);
                all_rate_unsmooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q1_rate_unsmooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q1_rate_smooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q2_rate_unsmooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q2_rate_smooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q3_rate_unsmooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q3_rate_smooth_by_group.(group) = nan(n_clusters, n_bins);
                all_bin_centers_by_group.(group) = plot_bin_centers;
                all_avg_velocity_by_group.(group) = nan(n_clusters, n_bins);
                all_occ_by_group.(group) = nan(n_clusters, n_bins);
                all_rate_per_trial_by_group.(group) = cell(n_clusters, 1);
                all_rate_per_trial_smooth_by_group.(group) = cell(n_clusters, 1);
                all_spike_positions_by_group.(group) = cell(n_clusters, 1);
                all_global_trial_indices_by_group.(group) = cell(n_clusters, 1);
                all_rate_pooled_by_group.(group) = nan(n_clusters, n_bins);
            end
            
            % --- Store results ---
            all_rate_smooth_by_group.(group)(c, :) = mean_rate_smooth;
            all_rate_unsmooth_by_group.(group)(c, :) = mean_rate_unsmooth;
            % Store unsmoothed quartiles (computed from unsmoothed rates)
            all_Q1_rate_unsmooth_by_group.(group)(c, :) = prctile(rate_per_trial, 25, 1);
            all_Q2_rate_unsmooth_by_group.(group)(c, :) = prctile(rate_per_trial, 50, 1);
            all_Q3_rate_unsmooth_by_group.(group)(c, :) = prctile(rate_per_trial, 75, 1);
            % Store smoothed quartiles (computed from smoothed rates)
            all_Q1_rate_smooth_by_group.(group)(c, :) = Q1_from_smoothed;
            all_Q2_rate_smooth_by_group.(group)(c, :) = Q2_from_smoothed;
            all_Q3_rate_smooth_by_group.(group)(c, :) = Q3_from_smoothed;
            all_bin_centers_by_group.(group) = plot_bin_centers;
            all_avg_velocity_by_group.(group)(c, :) = mean_velocity;
            all_occ_by_group.(group)(c, :) = mean_occupancy;
            all_rate_per_trial_by_group.(group){c} = rate_per_trial;
            all_rate_per_trial_smooth_by_group.(group){c} = rate_per_trial_smooth;
            all_spike_positions_by_group.(group){c} = spike_positions_per_trial;
            all_global_trial_indices_by_group.(group){c} = global_trial_indices;
            all_rate_pooled_by_group.(group)(c, :) = rate_pooled;
        end
    end
    
        % --- Save computed data to cache ---
        fprintf('  Saving analysis data to cache: %s\n', cache_filename);
        save(cache_filepath, 'all_rate_smooth_by_group', 'all_rate_unsmooth_by_group', ...
             'all_bin_centers_by_group', 'all_avg_velocity_by_group', ...
             'all_Q1_rate_smooth_by_group', 'all_Q2_rate_smooth_by_group', ...
             'all_Q3_rate_smooth_by_group', 'all_occ_by_group', ...
             'all_rate_per_trial_by_group', 'all_rate_per_trial_smooth_by_group', ...
             'all_spike_positions_by_group', 'all_global_trial_indices_by_group', ...
             'all_rate_pooled_by_group', ...
             'cluster_ids', 'n_clusters', 'group_names', 'group_labels', '-v7.3');
        fprintf('  Cache saved successfully\n');
    end  % End of analysis computation block

    % --- Plot average running and occupancy profiles for each group ---
    if save_figs
        fprintf('  Creating average running and occupancy profiles...\n');
        fig_profiles = ctl.figs.a4figure('landscape');
        % Average velocity
        subplot(2, 1, 1, 'Parent', fig_profiles);
        hold on;
        for g = 1:2
            group = group_names{g};
            bin_centers = all_bin_centers_by_group.(group);
            avg_vel = nanmean(all_avg_velocity_by_group.(group), 1);
            plot(bin_centers, avg_vel, 'LineWidth', 2, 'DisplayName', group_labels{g});
        end
        hold off;
        xlabel('Position (cm)');
        ylabel('Avg Velocity (cm/s)');
        legend('show', 'Location', 'northeastoutside');
        % Removed title to avoid overlap with FigureTitle
        grid on;
        % Average occupancy
        subplot(2, 1, 2, 'Parent', fig_profiles);
        hold on;
        for g = 1:2
            group = group_names{g};
            bin_centers = all_bin_centers_by_group.(group);
            avg_occ = nanmean(all_occ_by_group.(group), 1);
            plot(bin_centers, avg_occ, 'LineWidth', 2, 'DisplayName', group_labels{g});
        end
        hold off;
        xlabel('Position (cm)');
        ylabel('Occupancy (s)');
        legend('show', 'Location', 'northeastoutside');
        % Removed title to avoid overlap with FigureTitle
        grid on;
        
        % Add figure title
        FigureTitle(fig_profiles, sprintf('Average Profiles - %s', probe_ids{pid}));
        ctl.figs.save_fig_to_join();
    end

    % Determine session types based on actual trial data
    session_types = cell(length(sessions), 1);
    session_max_positions = zeros(length(sessions), 1);
    
    for s = 1:length(sessions)
        session = sessions{s};
        trials_in_session = session.trials;
        session_max_pos = 0;
        
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            if ~isempty(pos)
                max_pos = max(pos);
                session_max_pos = max(session_max_pos, max_pos);
            end
        end
        
        session_max_positions(s) = session_max_pos;
        if session_max_pos > 90
            session_types{s} = 'long';
        else
            session_types{s} = 'short';
        end
        
        fprintf('  Session %d: max position = %.2f cm, type = %s\n', s, session_max_pos, session_types{s});
    end
    
    % Define colors for different sessions
    % Long sessions (120cm): cold colors (blue, cyan)
    % Short sessions (60cm): warm colors (red, orange)
    cold_colors = {[0 0 0.8], [0 0.8 0.8]}; % dark blue, cyan
    warm_colors = {[0.8 0 0], [1 0.5 0]}; % dark red, orange
    
    session_colors = cell(length(sessions), 1);
    session_labels = {};
    long_count = 0;
    short_count = 0;
    
    for s = 1:length(sessions)
        if strcmp(session_types{s}, 'long')
            long_count = long_count + 1;
            session_colors{s} = cold_colors{mod(long_count-1, length(cold_colors)) + 1};
            session_labels{s} = sprintf('120cm Session %d', s);
        else
            short_count = short_count + 1;
            session_colors{s} = warm_colors{mod(short_count-1, length(warm_colors)) + 1};
            session_labels{s} = sprintf('60cm Session %d', s);
        end
    end

    % For each cluster, create a combined figure with spatial profile and x-normalized comparison
    if plot_single_cluster_fig
        fprintf('  Creating individual cluster plots...\n');
        
        % Define colors for the two groups
        group_colors = struct('long', [0 0 0.8], 'short', [0.8 0 0]); % Blue for long, red for short
        
        for c = 1:n_clusters
            cluster = clusters(c);
            
            % Prepare rate data structure for this cluster
            rate_data = struct();
            for g = 1:2
                group = group_names{g};
                rate_data.(group).Q1_smooth = all_Q1_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).Q2_smooth = all_Q2_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).Q3_smooth = all_Q3_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).rate_per_trial = all_rate_per_trial_by_group.(group){c};
                rate_data.(group).rate_per_trial_smooth = all_rate_per_trial_smooth_by_group.(group){c};
                rate_data.(group).spike_positions = all_spike_positions_by_group.(group){c};
                rate_data.(group).global_trial_indices = all_global_trial_indices_by_group.(group){c};
                rate_data.(group).rate_pooled = all_rate_pooled_by_group.(group)(c, :);
            end
            
            % Create combined figure using helper function
            fig_cluster = plot_cluster_spatial_profile(cluster.id, all_bin_centers_by_group, ...
                                                       rate_data, group_names, group_labels, ...
                                                       group_colors, probe_ids{pid});
            
            % Save using the figure management system
             if save_figs
                % Convert to a4figure format for PDF joining
                set(fig_cluster, 'PaperOrientation', 'landscape');
                set(fig_cluster, 'PaperUnits', 'normalized');
                set(fig_cluster, 'PaperPosition', [0 0 1 1]);
                ctl.figs.save_fig_to_join();
            end
            
            if ~save_figs
                close(fig_cluster);
            end
        end
    end
    
    % Plot heatmaps for each group in a single figure with four subplots
    if plot_heatmp_cluster_fig

        fprintf('  Creating heatmap plots...\n');
        fig_heatmaps = ctl.figs.a4figure('landscape');
    
        % Sort clusters by maximum peak positional firing rate for long trials
        long_rate_mat = all_rate_smooth_by_group.('long');
        peak_positions = zeros(n_clusters, 1);
        for c = 1:n_clusters
            rate_profile = long_rate_mat(c, :);
            bin_centers = all_bin_centers_by_group.('long');
            
            % Find the position of the maximum peak
            [~, max_idx] = max(rate_profile);
            peak_positions(c) = bin_centers(max_idx);
        end
        [~, sort_idx] = sort(peak_positions);
        sorted_cluster_ids = cluster_ids(sort_idx);
    
        for g = 1:2
            group = group_names{g};
            rate_mat = all_rate_smooth_by_group.(group);
            bin_centers = all_bin_centers_by_group.(group);
            sorted_rate_mat = rate_mat(sort_idx, :);
            
            % Raw heatmap
            subplot(2, 2, g, 'Parent', fig_heatmaps);
            imagesc(bin_centers, 1:n_clusters, sorted_rate_mat);
            colormap('jet'); % Changed to blue-to-yellow colormap
            colorbar;
            xlabel('Position (cm)');
            ylabel('Cluster (sorted by maximum peak position)');
            % Removed title to avoid overlap with FigureTitle
            set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sorted_cluster_ids);
            
            % Set x-axis limits to make short trials visually half width
            if strcmp(group, 'short')
                xlim([60, 120]); % Short trials span 60-120 cm
            else
                xlim([0, 120]); % Long trials span 0-120 cm
            end
            
            % Normalized heatmap
            norm_rate_mat = sorted_rate_mat ./ max(sorted_rate_mat, [], 2);
            norm_rate_mat(isnan(norm_rate_mat)) = 0;
            subplot(2, 2, g+2, 'Parent', fig_heatmaps);
            imagesc(bin_centers, 1:n_clusters, norm_rate_mat);
            colormap('jet'); % Changed to blue-to-yellow colormap
            colorbar;
            xlabel('Position (cm)');
            ylabel('Cluster (sorted by maximum peak position)');
            % Removed title to avoid overlap with FigureTitle
            set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sorted_cluster_ids);
            
            % Set x-axis limits to make short trials visually half width
            if strcmp(group, 'short')
                xlim([60, 120]); % Short trials span 60-120 cm
            else
                xlim([0, 120]); % Long trials span 0-120 cm
            end
        end
        % Removed sgtitle to avoid overlap with FigureTitle
        
        % Add figure title with cleaner formatting
        FigureTitle(fig_heatmaps, sprintf('Heatmaps - %s', probe_ids{pid}));
        ctl.figs.save_fig_to_join();
        
        % Join all plots for this probe into one PDF
        fprintf('  Saving PDF for probe %s...\n', probe_ids{pid});
        fname = sprintf('%s.pdf', probe_ids{pid});
        ctl.figs.join_figs(fname, overwrite);
        ctl.figs.clear_figs();
        
        fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
    end
end

fprintf('\nSpatial firing rate profile analysis completed for %d probe(s)\n', length(probe_ids)); 