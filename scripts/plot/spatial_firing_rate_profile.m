% Plot spatial firing rate profile for each cluster using position as x-axis
% For experiment group 'training_running', bin position into 2-cm bins,
% normalize spike counts by occupancy, and smooth with a Gaussian kernel (8-cm s.d.)
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
save_figs               = false;
overwrite               = false;
figure_dir              = {'spatial_firing_rate', 'ambient_light'};
plot_single_cluster_fig = false;
plot_heatmp_cluster_fig = false;


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
            
            for k = 1:n_trials
                trial = trials_struct(k).trial;
                [rate, occ, vel] = compute_trial_firing_rate(trial, cluster, edges);
                rate_per_trial(k, :) = rate;
                occ_per_trial(k, :) = occ;
                vel_per_trial(k, :) = vel;
            end
            
            % --- Compute quartiles and smooth using helper function ---
            fprintf('\n    Cluster %d, %s trials:\n', cluster.id, group);
            quartile_results = compute_quartiles_and_smooth(rate_per_trial, bin_size_cm, gauss_sigma_cm, true);
            
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
            end
            
            % --- Store results ---
            all_rate_smooth_by_group.(group)(c, :) = quartile_results.mean_smooth;
            all_rate_unsmooth_by_group.(group)(c, :) = quartile_results.mean_unsmooth;
            all_Q1_rate_unsmooth_by_group.(group)(c, :) = quartile_results.Q1_unsmooth;
            all_Q1_rate_smooth_by_group.(group)(c, :) = quartile_results.Q1_smooth;
            all_Q2_rate_unsmooth_by_group.(group)(c, :) = quartile_results.Q2_unsmooth;
            all_Q2_rate_smooth_by_group.(group)(c, :) = quartile_results.Q2_smooth;
            all_Q3_rate_unsmooth_by_group.(group)(c, :) = quartile_results.Q3_unsmooth;
            all_Q3_rate_smooth_by_group.(group)(c, :) = quartile_results.Q3_smooth;
            all_bin_centers_by_group.(group) = plot_bin_centers;
            all_avg_velocity_by_group.(group)(c, :) = mean_velocity;
            all_occ_by_group.(group)(c, :) = mean_occupancy;
        end
    end

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
        legend('show');
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
        legend('show');
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
                rate_data.(group).mean_smooth = all_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).Q1_smooth = all_Q1_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).Q2_smooth = all_Q2_rate_smooth_by_group.(group)(c, :);
                rate_data.(group).Q3_smooth = all_Q3_rate_smooth_by_group.(group)(c, :);
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