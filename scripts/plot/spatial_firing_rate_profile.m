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
    
    % Session analysis and trial counting
    for s = 1:length(sessions)
        session = sessions{s};
        fprintf('  Session %d/%d: session_id = %s\n', s, length(sessions), session.session_id);
        trials_in_session = session.trials;
        n_in_motion_trials = 0;
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            if ~isempty(pos)
                n_in_motion_trials = n_in_motion_trials + 1;
                max_pos = max(pos);
                fprintf('    Trial %d: max(position in-motion) = %.2f cm\n', t, max_pos);
            else
                fprintf('    Trial %d: [no in-motion samples]\n', t);
            end
        end
        fprintf('    Total in-motion trials in session: %d\n', n_in_motion_trials);
    end
 

    % --- Collect all trials with their session index and max in-motion position ---
    all_trials_info = struct('trial', {}, 'session_idx', {}, 'trial_idx', {}, 'max_pos', {});
    for s = 1:length(sessions)
        session = sessions{s};
        trials_in_session = session.trials;
        for t = 1:length(trials_in_session)
            trial = trials_in_session{t}.to_aligned;
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            if ~isempty(pos)
                max_pos = max(pos);
                all_trials_info(end+1) = struct('trial', trial, 'session_idx', s, 'trial_idx', t, 'max_pos', max_pos); %#ok<AGROW>
            end
        end
    end
    
    % --- Group trials by max position (long/short) ---
    long_trials = all_trials_info([all_trials_info.max_pos] > 90);
    short_trials = all_trials_info([all_trials_info.max_pos] <= 90);
    fprintf('  Total long trials: %d\n', numel(long_trials));
    fprintf('  Total short trials: %d\n', numel(short_trials));

    % --- Prepare for cluster analysis as before, but using long_trials and short_trials ---
    trial_groups = struct('long', long_trials, 'short', short_trials);
    group_names = {'long', 'short'};
    group_labels = {'Long (0-120 cm)', 'Short (60-120 cm)'};
    all_rate_smooth_by_group = struct();
    all_rate_unsmooth_by_group = struct(); % NEW: Store unsmoothed rates
    all_bin_centers_by_group = struct();
    all_avg_velocity_by_group = struct();

    all_Q1_rate_unsmooth_by_group = struct();
    all_Q1_rate_smooth_by_group = struct();
    all_Q2_rate_unsmooth_by_group = struct();  % Median
    all_Q2_rate_smooth_by_group = struct();
    all_Q3_rate_unsmooth_by_group = struct();
    all_Q3_rate_smooth_by_group = struct();

    all_occ_by_group = struct();
    cluster_ids = [clusters.id];
    n_clusters = length(clusters);

    % Compute and fill group summary structures for all clusters (needed for later plots)
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
                edges = 0:bin_size_cm:60; % Short trials have positions 0-60 cm
            end
            bin_centers = edges(1:end-1) + bin_size_cm/2;
            n_bins = length(bin_centers);
            
            % For short trials, adjust bin centers for plotting to align with 60-120 cm
            if strcmp(group, 'short')
                plot_bin_centers = bin_centers + 60; % Shift to 60-120 cm for plotting
            else
                plot_bin_centers = bin_centers;
            end
            
            % --- Trial-by-trial firing rate computation ---
            n_trials = length(trials_struct);
            rate_per_trial = nan(n_trials, n_bins);  % Each row is one trial
            occ_per_trial = nan(n_trials, n_bins);
            vel_per_trial = nan(n_trials, n_bins);
            
            for k = 1:n_trials
                trial = trials_struct(k).trial;
                motion_mask = trial.motion_mask();
                pos = trial.position(motion_mask);
                tvec = trial.probe_t(motion_mask);
                vel = trial.velocity(motion_mask);
                
                if isempty(pos)
                    continue;
                end
                
                % Get spike times for this trial
                st = cluster.spike_times;
                mask = st >= tvec(1) & st <= tvec(end);
                st = st(mask);
                
                % Interpolate spike positions
                if ~isempty(st)
                    spike_pos = interp1(tvec, pos, st, 'linear', 'extrap');
                else
                    spike_pos = [];
                end
                
                % Compute spike count per bin for this trial
                spike_count = histcounts(spike_pos, edges);
                
                % Compute occupancy per bin for this trial
                dt_vec = [diff(tvec); 0];
                occupancy = zeros(1, n_bins);
                velocity_sum = zeros(1, n_bins);
                for i = 1:n_bins
                    in_bin = pos >= edges(i) & pos < edges(i+1);
                    occupancy(i) = sum(dt_vec(in_bin));
                    velocity_sum(i) = sum(vel(in_bin) .* dt_vec(in_bin));
                end
                
                % Compute firing rate for this trial
                rate = spike_count ./ occupancy;
                rate(isnan(rate) | isinf(rate)) = NaN;  % Keep NaN for bins with no occupancy
                
                avg_velocity = velocity_sum ./ occupancy;
                avg_velocity(isnan(avg_velocity) | isinf(avg_velocity)) = NaN;
                
                rate_per_trial(k, :) = rate;
                occ_per_trial(k, :) = occupancy;
                vel_per_trial(k, :) = avg_velocity;
            end
            
            % --- Compute quartiles across trials ---
            % Q1 (25th percentile), Q2 (median, 50th percentile), Q3 (75th percentile)
            Q1_rate = prctile(rate_per_trial, 25, 1);
            Q2_rate = prctile(rate_per_trial, 50, 1);  % Median
            Q3_rate = prctile(rate_per_trial, 75, 1);
            
            % Also compute mean for comparison
            mean_rate = nanmean(rate_per_trial, 1);
            
            % Compute mean occupancy and velocity across trials
            mean_occupancy = nanmean(occ_per_trial, 1);
            mean_velocity = nanmean(vel_per_trial, 1);
            
            % Replace NaN with 0 for smoothing
            Q1_rate(isnan(Q1_rate)) = 0;
            Q2_rate(isnan(Q2_rate)) = 0;
            Q3_rate(isnan(Q3_rate)) = 0;
            mean_rate(isnan(mean_rate)) = 0;
            
            % --- Print statistics for this cluster and group ---
            fprintf('\n    Cluster %d, %s trials (n=%d):\n', cluster.id, group, n_trials);
            fprintf('      Firing rate Q1 (25%%): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(Q1_rate), min(Q1_rate), max(Q1_rate));
            fprintf('      Firing rate Q2 (50%%, median): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(Q2_rate), min(Q2_rate), max(Q2_rate));
            fprintf('      Firing rate Q3 (75%%): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(Q3_rate), min(Q3_rate), max(Q3_rate));
            fprintf('      Mean firing rate: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(mean_rate), min(mean_rate), max(mean_rate));
            
            % --- Create proper spatial Gaussian kernel for smoothing ---
            kernel_size_cm = 8 * gauss_sigma_cm; % 64 cm total kernel size
            kernel_samples = round(kernel_size_cm / bin_size_cm); 
            x = linspace(-kernel_samples/2, kernel_samples/2, kernel_samples);
            kernel = exp(-x.^2 / (2 * (gauss_sigma_cm/bin_size_cm)^2));
            kernel = kernel / sum(kernel); % normalize
            
            % Pad and smooth each quartile
            pad_size = floor(length(kernel) / 2);
            
            % Smooth Q1
            Q1_padded = [fliplr(Q1_rate(1:pad_size)), Q1_rate, fliplr(Q1_rate(end-pad_size+1:end))];
            Q1_smooth_padded = conv(Q1_padded, kernel, 'same');
            Q1_rate_smooth = Q1_smooth_padded(pad_size + 1 : end - pad_size);
            
            % Smooth Q2 (median)
            Q2_padded = [fliplr(Q2_rate(1:pad_size)), Q2_rate, fliplr(Q2_rate(end-pad_size+1:end))];
            Q2_smooth_padded = conv(Q2_padded, kernel, 'same');
            Q2_rate_smooth = Q2_smooth_padded(pad_size + 1 : end - pad_size);
            
            % Smooth Q3
            Q3_padded = [fliplr(Q3_rate(1:pad_size)), Q3_rate, fliplr(Q3_rate(end-pad_size+1:end))];
            Q3_smooth_padded = conv(Q3_padded, kernel, 'same');
            Q3_rate_smooth = Q3_smooth_padded(pad_size + 1 : end - pad_size);
            
            % Smooth mean rate (for backward compatibility)
            mean_padded = [fliplr(mean_rate(1:pad_size)), mean_rate, fliplr(mean_rate(end-pad_size+1:end))];
            mean_smooth_padded = conv(mean_padded, kernel, 'same');
            rate_smooth = mean_smooth_padded(pad_size + 1 : end - pad_size);
            
            % --- Print smoothed statistics ---
            fprintf('      After smoothing (sigma=%.1f cm):\n', gauss_sigma_cm);
            fprintf('        Smoothed Q1: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(Q1_rate_smooth), min(Q1_rate_smooth), max(Q1_rate_smooth));
            fprintf('        Smoothed Q2 (median): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(Q2_rate_smooth), min(Q2_rate_smooth), max(Q2_rate_smooth));
            fprintf('        Smoothed Q3: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(Q3_rate_smooth), min(Q3_rate_smooth), max(Q3_rate_smooth));
            fprintf('        Smoothed mean: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
                    nanmean(rate_smooth), min(rate_smooth), max(rate_smooth));
            
            % --- Initialize storage structures if not already done ---
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
            all_rate_smooth_by_group.(group)(c, :) = rate_smooth;
            all_rate_unsmooth_by_group.(group)(c, :) = mean_rate;
            all_Q1_rate_unsmooth_by_group.(group)(c, :) = Q1_rate;
            all_Q1_rate_smooth_by_group.(group)(c, :) = Q1_rate_smooth;
            all_Q2_rate_unsmooth_by_group.(group)(c, :) = Q2_rate;
            all_Q2_rate_smooth_by_group.(group)(c, :) = Q2_rate_smooth;
            all_Q3_rate_unsmooth_by_group.(group)(c, :) = Q3_rate;
            all_Q3_rate_smooth_by_group.(group)(c, :) = Q3_rate_smooth;
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

    % For each cluster, create a single figure with both group firing rates overlaid
    if plot_single_cluster_fig
        fprintf('  Creating individual cluster plots with quartile distributions...\n');
        
        % Define colors for the two groups
        group_colors = struct('long', [0 0 0.8], 'short', [0.8 0 0]); % Blue for long, red for short
        
        for c = 1:n_clusters
            cluster = clusters(c);
            fig_cluster = ctl.figs.a4figure('landscape');
            
            hold on;
            
            % Plot each group (long and short)
            for g = 1:2
                group = group_names{g};
                
                % Get the precomputed data for this cluster and group
                bin_centers = all_bin_centers_by_group.(group);
                rate_smooth = all_rate_smooth_by_group.(group)(c, :);  % Mean (smoothed)
                Q1_smooth = all_Q1_rate_smooth_by_group.(group)(c, :);
                Q2_smooth = all_Q2_rate_smooth_by_group.(group)(c, :);  % Median (smoothed)
                Q3_smooth = all_Q3_rate_smooth_by_group.(group)(c, :);
                
                % Skip if all NaN
                if all(isnan(rate_smooth))
                    continue;
                end
                
                % Plot Q1-Q3 shaded area (interquartile range)
                fill_x = [bin_centers, fliplr(bin_centers)];
                fill_y = [Q1_smooth, fliplr(Q3_smooth)];
                fill(fill_x, fill_y, group_colors.(group), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                     'HandleVisibility', 'off');
                
                % Plot median (Q2) as main line
                plot(bin_centers, Q2_smooth, '-', 'LineWidth', 2, 'Color', group_colors.(group), ...
                     'DisplayName', sprintf('%s median', group_labels{g}));
                
                % Plot mean as dashed line for comparison
                plot(bin_centers, rate_smooth, '--', 'LineWidth', 1.5, 'Color', group_colors.(group), ...
                     'DisplayName', sprintf('%s mean', group_labels{g}));
            end
            
            hold off;
            xlabel('Position (cm)');
            ylabel('Firing rate (Hz)');
            legend('show', 'Location', 'best');
            grid on;
            
            % Add figure title
            FigureTitle(fig_cluster, sprintf('Cluster %d - %s (shaded: IQR)', cluster.id, probe_ids{pid}));
            ctl.figs.save_fig_to_join();
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

    %  Normalized x-position with respect to firing rate
    % === Make sure required variables exist (from your main script) ===
    if ~exist('all_rate_smooth_by_group','var') || ...
       ~isfield(all_rate_smooth_by_group,'long') || ...
       ~isfield(all_rate_smooth_by_group,'short') || ...
       ~exist('all_bin_centers_by_group','var') || ...
       ~exist('clusters','var')
        error('Required variables not found. Run your main script first.');
    end
    
    % Number of clusters from the long group matrix
    n_clusters = size(all_rate_smooth_by_group.long, 1);
    
    % Decide output PDF name based on the last processed probe
    if exist('probe_ids','var') && exist('pid','var')
        out_pdf = sprintf('%s_shape_comparison_xnorm.pdf', probe_ids{pid});
    else
        out_pdf = 'shape_comparison_xnorm.pdf';
    end
    
    % If file already exists, delete so we can recreate it
    if exist(out_pdf, 'file')
        delete(out_pdf);
    end
    
    fprintf('Saving x-normalized shape comparison for %d clusters to %s\n', ...
            n_clusters, out_pdf);
    
    % Loop over all clusters
    for c = 1:n_clusters
        cluster_id = clusters(c).id;
    
        % Extract x and y for this cluster
        pos_long   = all_bin_centers_by_group.long(:);            % 0–120
        rate_long  = all_rate_smooth_by_group.long(c,:).';        % row -> column
    
        pos_short  = all_bin_centers_by_group.short(:);           % 60–120 (shifted)
        rate_short = all_rate_smooth_by_group.short(c,:).';       % row -> column
    
        % Skip if this cluster has no data in one of the groups
        if all(isnan(rate_long)) || all(isnan(rate_short))
            fprintf('Cluster %d: skipped (no data in one group)\n', cluster_id);
            continue;
        end
    
        % Normalized position (start=0, end=1) for each condition
        pos_long_norm  = (pos_long  - min(pos_long))  / (max(pos_long)  - min(pos_long));
        pos_short_norm = (pos_short - min(pos_short)) / (max(pos_short) - min(pos_short));
    
        % Keep only finite points
        idxL = isfinite(pos_long_norm)  & isfinite(rate_long);
        idxS = isfinite(pos_short_norm) & isfinite(rate_short);
    
        pos_long_n  = pos_long_norm(idxL);
        rate_long_n = rate_long(idxL);
    
        pos_short_n  = pos_short_norm(idxS);
        rate_short_n = rate_short(idxS);
    
        % Need at least two points in each group for interpolation
        if numel(pos_long_n) < 2 || numel(pos_short_n) < 2
            fprintf('Cluster %d: skipped (not enough data after filtering)\n', cluster_id);
            continue;
        end
    
        % Common normalized grid (x-axis only)
        x_norm = linspace(0,1,100);
    
        % Interpolate both curves onto the same normalized grid
        rate_long_i  = interp1(pos_long_n,  rate_long_n,  x_norm, 'linear', 'extrap');
        rate_short_i = interp1(pos_short_n, rate_short_n, x_norm, 'linear', 'extrap');
    
        % IMPORTANT: keep firing rates RAW (no amplitude normalization here)
        % If you do want shape-only correlation later, you can z-score just for r.
    
        % Compute correlation as a shape/similarity metric on RAW FR
        r = corr(rate_long_i(:), rate_short_i(:), 'rows', 'complete');
    
        % --- Plot and append to the PDF ---
        fig = figure('Visible','off');  % or 'on' if you want to see them
        hold on;
        plot(x_norm, rate_long_i,  'b', 'LineWidth', 2);
        plot(x_norm, rate_short_i, 'r', 'LineWidth', 2);
        hold off;
    
        xlabel('Normalized track position (start \rightarrow end)');
        ylabel('Firing rate (Hz)');   % raw firing rate
        legend({'Long (0–120 cm)','Short (60–120 cm)'}, 'Location', 'best');
        title(sprintf('Cluster %d – x-normalized comparison (r = %.2f)', cluster_id, r));
        grid on;
    
        % Append this figure as a new page in the PDF
        exportgraphics(fig, out_pdf, 'Append', true);
        close(fig);
    
        fprintf('Cluster %d: added to PDF (r = %.2f)\n', cluster_id, r);
    end
end

fprintf('\nSpatial firing rate profile analysis completed for %d probe(s)\n', length(probe_ids)); 