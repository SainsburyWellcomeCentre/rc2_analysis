% ============================
% PLACE CELL IDENTIFICATION (integrated)
% Neuropixels – RC2Analysis pipeline (per-cluster, per-group classification)
% =============================

clear all; 
clc;

output_base_dir = '/Users/randybentbarker/Documents/SWC/spatial_firing_rate_profile_identification';


% Configuration
experiment_groups       = {'ambient_light'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'spatial_firing_rate', 'ambient_light'};


% Analysis parameters (you can tweak)
bin_size_cm       = 2;
gauss_sigma_cm    = 8;    % standard deviation for smoothing (cm)
minSpikes         = 25;   % minimum spikes per cluster per group to attempt classification
minPeakRate       = 0.5;  % Hz minimum peak for field consideration
fieldFrac         = 0.2;  % field threshold fraction of peak
minFieldBins      = 3;    % contiguous bins for a field
nShuf             = 100; % number of circular-shift shuffles
pThresh           = 0.05; % significance threshold for Skaggs shuffle
% stabilityMinR     = 0.25; % minimum split-half correlation required

% Initialize controller and get probe IDs
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));

% Prepare place_cell_results container (per probe)
place_cell_results = struct('probe_id', {}, 'clusters', {});

for pid = 1:length(probe_ids)
    fprintf('\nProcessing probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    
    data = ctl.load_formatted_data(probe_ids{pid});
    clusters = data.selected_clusters();
    sessions = data.motion_sessions();
    
    fprintf('  Found %d clusters and %d sessions\n', length(clusters), length(sessions));
    
    % Session analysis and trial counting (unchanged)
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
    all_rate_unsmooth_by_group = struct(); % Store unsmoothed rates
    all_bin_centers_by_group = struct();
    all_avg_velocity_by_group = struct();
    all_occ_by_group = struct();
    cluster_ids = [clusters.id];
    n_clusters = length(clusters);

    % Initialize place cell results for this probe
    place_cell_results(pid).probe_id  = probe_ids{pid};
    place_cell_results(pid).clusters  = repmat(struct(), n_clusters, 1);
    % Compute and fill group summary structures for all clusters
    fprintf('  Computing spatial firing rate profiles for %d clusters...\n', n_clusters);
    for c = 1:n_clusters
        cluster = clusters(c);
        % Initialize per-cluster place_cell_results entries
        place_cell_results(pid).clusters(c).cluster_id = cluster.id;
        place_cell_results(pid).clusters(c).group = struct(); % will contain .long and .short

        for g = 1:2
            group = group_names{g};
            trials_struct = trial_groups.(group);
            all_spike_positions = [];
            all_positions = [];
            all_times = [];
            all_velocities = [];
            all_spike_times = []; % collect spike times used for shuffles / stability
            for k = 1:length(trials_struct)
                trial = trials_struct(k).trial;
                motion_mask = trial.motion_mask();
                pos = trial.position(motion_mask);
                tvec = trial.probe_t(motion_mask);
                vel = trial.velocity(motion_mask);
                all_positions = [all_positions; pos(:)];
                all_times = [all_times; tvec(:)];
                all_velocities = [all_velocities; vel(:)];
                st = cluster.spike_times;
                mask = st >= tvec(1) & st <= tvec(end);
                st_in = st(mask);
                if ~isempty(st_in)
                    spike_pos = interp1(tvec, pos, st_in, 'linear', 'extrap');
                    all_spike_positions = [all_spike_positions; spike_pos];
                    all_spike_times = [all_spike_times; st_in(:)]; %#ok<AGROW>
                end
            end
            % If no position samples for this group, continue
            if isempty(all_positions)
                % store NaNs to keep matrix shapes consistent later
                if ~isfield(all_rate_smooth_by_group, group)
                    % set up shapes with n_clusters x 0 placeholder; will fill later when data exists
                    all_rate_smooth_by_group.(group) = nan(n_clusters, 0);
                    all_rate_unsmooth_by_group.(group) = nan(n_clusters, 0);
                    all_bin_centers_by_group.(group) = [];
                    all_avg_velocity_by_group.(group) = nan(n_clusters, 0);
                    all_occ_by_group.(group) = nan(n_clusters, 0);
                end
                % Mark place result as no_data
                place_cell_results(pid).clusters(c).group.(group) = struct( ...
                    'isPlace', false, 'reason', 'no_position_samples', 'info', NaN, 'pVal', NaN, ...
                    'peakRate', NaN, 'fieldBins', 0, 'stabilityR', NaN);
                continue;
            end

            % Set bin edges explicitly for each group
            if strcmp(group, 'long')
                edges = 0:bin_size_cm:120;
            else
                edges = 0:bin_size_cm:60; % Short trials have positions 0-60 cm
            end
            bin_centers = edges(1:end-1) + bin_size_cm/2;
            
            % For short trials, adjust bin centers for plotting to align with 60-120 cm
            if strcmp(group, 'short')
                plot_bin_centers = bin_centers + 60; % Shift to 60-120 cm for plotting
            else
                plot_bin_centers = bin_centers;
            end
            
            % Compute spike_count and occupancy (as before)
            spike_count = histcounts(all_spike_positions, edges);
            dt_vec = [diff(all_times); 0];
            occupancy = zeros(1, length(edges)-1);
            velocity_sum = zeros(1, length(edges)-1);
            for i = 1:length(occupancy)
                in_bin = all_positions >= edges(i) & all_positions < edges(i+1);
                occupancy(i) = sum(dt_vec(in_bin));
                velocity_sum(i) = sum(all_velocities(in_bin) .* dt_vec(in_bin));
            end
            avg_velocity = velocity_sum ./ occupancy;
            avg_velocity(isnan(avg_velocity) | isinf(avg_velocity)) = 0;
            rate = spike_count ./ occupancy;
            rate(isnan(rate) | isinf(rate)) = 0;

            % Create proper spatial Gaussian kernel (keep your original kernel construction)
            kernel_size_cm = 8 * gauss_sigma_cm; % 64 cm total kernel size
            kernel_samples = round(kernel_size_cm / bin_size_cm); 
            if kernel_samples < 3
                kernel_samples = 3; 
            end
            x = linspace(-kernel_samples/2, kernel_samples/2, kernel_samples);
            kernel = exp(-x.^2 / (2 * (gauss_sigma_cm/bin_size_cm)^2));
            kernel = kernel / sum(kernel); % normalize

            % Pad the rate data to avoid edge effects using mirror padding
            pad_size = floor(length(kernel) / 2);
            % Make sure rate has length > pad_size; otherwise pad more conservatively
            if length(rate) <= pad_size
                pad_size = max(1, floor(length(rate)/2));
            end
            % --- Smooth occupancy (shared by obs + shuffles) ---
            occ_padded        = [fliplr(occupancy(1:pad_size)), occupancy, fliplr(occupancy(end-pad_size+1:end))];
            smoothOcc_padded  = conv(occ_padded, kernel, 'same');
            smoothOcc         = smoothOcc_padded(pad_size + 1 : end - pad_size);
            smoothOcc(smoothOcc==0) = NaN;
            
            % --- Smooth spike COUNTS and divide by smoothOcc: OBSERVED rate_smooth ---
            counts_padded       = [fliplr(spike_count(1:pad_size)), spike_count, fliplr(spike_count(end-pad_size+1:end))];
            smoothCounts_padded = conv(counts_padded, kernel, 'same');
            smoothCounts        = smoothCounts_padded(pad_size + 1 : end - pad_size);
            rate_smooth         = smoothCounts ./ smoothOcc;
            rate_smooth(isnan(rate_smooth)) = 0;
            
            % Store for heatmap and group plots
            if ~isfield(all_rate_smooth_by_group, group)
                all_rate_smooth_by_group.(group) = nan(n_clusters, length(rate_smooth));
                all_rate_unsmooth_by_group.(group) = nan(n_clusters, length(rate)); % Store unsmoothed rates
                all_bin_centers_by_group.(group) = plot_bin_centers; % Use plot_bin_centers for display
                all_avg_velocity_by_group.(group) = nan(n_clusters, length(avg_velocity));
                all_occ_by_group.(group) = nan(n_clusters, length(occupancy));
            end
            all_rate_smooth_by_group.(group)(c, 1:length(rate_smooth)) = rate_smooth;
            all_rate_unsmooth_by_group.(group)(c, 1:length(rate)) = rate; % Store unsmoothed rates
            all_bin_centers_by_group.(group) = plot_bin_centers; % Use plot_bin_centers for display
            all_avg_velocity_by_group.(group)(c, 1:length(avg_velocity)) = avg_velocity;
            all_occ_by_group.(group)(c, 1:length(occupancy)) = occupancy;


            % -----------------------------
            % PLACE CELL IDENTIFICATION
            % -----------------------------
            % Build Skaggs information for this cluster and group
            infoObs = skaggs_info(rate_smooth, smoothOcc);

            % If too few spikes in this group, skip classification
            n_spikes_group = numel(all_spike_times);
            if n_spikes_group < minSpikes
                place_cell_results(pid).clusters(c).group.(group) = struct( ...
                    'isPlace', false, 'reason', 'low_spikes', 'info', infoObs, 'pVal', NaN, ...
                    'peakRate', max(rate_smooth), 'fieldBins', 0, 'stabilityR', NaN);
                continue;
            end

            % Build shuffle null distribution (circular shift)
            % Define time window for circular shift using the union of all_times
            t_min = min(all_times);
            t_max = max(all_times);
            sessionDur = t_max - t_min;
            if sessionDur <= 0
                % cannot shuffle; mark as error
                place_cell_results(pid).clusters(c).group.(group) = struct( ...
                    'isPlace', false, 'reason', 'zero_time_span', 'info', infoObs, 'pVal', NaN, ...
                    'peakRate', max(rate_smooth), 'fieldBins', 0, 'stabilityR', NaN);
                continue;
            end

            infoNull = nan(nShuf,1);
           

            for s = 1:nShuf
                shift = rand * sessionDur;
                spkSh = mod(all_spike_times - t_min + shift, sessionDur) + t_min;
                % map shifted spike times back to positions using same all_times/all_positions mapping
                % interp1 requires monotonic all_times. If repeating times occur, use 'previous'/'next' would be problematic.
                % Use linear interpolation
                spkPosSh = interp1(all_times, all_positions, spkSh, 'linear', 'extrap');
                spike_count_sh = histcounts(spkPosSh, edges);
                % smooth counts & divide by smooth occupancy
                % pad counts, convolve, unpad
                counts_padded = [fliplr(spike_count_sh(1:pad_size)), spike_count_sh, fliplr(spike_count_sh(end-pad_size+1:end))];
                smoothCounts_padded = conv(counts_padded, kernel, 'same');
                smoothCounts = smoothCounts_padded(pad_size + 1 : end - pad_size);
                rate_sh = smoothCounts ./ smoothOcc;
                rate_sh(isnan(rate_sh)) = 0;
                infoNull(s) = skaggs_info(rate_sh, smoothOcc);
            end

            pVal = (sum(infoNull >= infoObs) + 1) / (nShuf + 1);

            % Peak and field detection (on smoothed rate)
            peakRate = max(rate_smooth);
            if peakRate < minPeakRate
                place_cell_results(pid).clusters(c).group.(group) = struct( ...
                    'isPlace', false, 'reason', 'low_peak', 'info', infoObs, 'pVal', pVal, ...
                    'peakRate', peakRate, 'fieldBins', 0, 'stabilityR', NaN);
                continue;
            end

            thresh = fieldFrac * peakRate;
            above = find(rate_smooth >= thresh);
            if isempty(above)
                place_cell_results(pid).clusters(c).group.(group) = struct( ...
                    'isPlace', false, 'reason', 'no_field', 'info', infoObs, 'pVal', pVal, ...
                    'peakRate', peakRate, 'fieldBins', 0, 'stabilityR', NaN);
                continue;
            end

            % find contiguous segments
            d = diff(above);
            breaks = [0, find(d > 1), length(above)];
            segLens = zeros(length(breaks)-1,1);
            for k = 1:length(segLens)
                idx = (breaks(k)+1) : breaks(k+1);
                seg = above(idx);
                segLens(k) = length(seg);
            end
            longestField = max(segLens);
            if longestField < minFieldBins
                place_cell_results(pid).clusters(c).group.(group) = struct( ...
                    'isPlace', false, 'reason', 'small_field', 'info', infoObs, 'pVal', pVal, ...
                    'peakRate', peakRate, 'fieldBins', longestField, 'stabilityR', NaN);
                continue;
            end

            % Split-half stability: split by median time of all_times
            mid = (t_min + t_max) / 2;
            % split spike times into two sets relative to mid
            spk1 = all_spike_times(all_spike_times < mid);
            spk2 = all_spike_times(all_spike_times >= mid);
            % compute rates for spk1/spk2 with same binning & smoothing
            rm1 = compute_rate_from_spikes(spk1, all_times, all_positions, edges, kernel, pad_size, smoothOcc);
            rm2 = compute_rate_from_spikes(spk2, all_times, all_positions, edges, kernel, pad_size, smoothOcc);

            validBins = ~isnan(rm1) & ~isnan(rm2);
            if sum(validBins) < 5
                stabilityR = NaN;
            else
                stabilityR = corr(rm1(validBins)', rm2(validBins)');
            end

            % Final decision
            % isPlace = (pVal < pThresh) && (longestField >= minFieldBins) && (~isnan(stabilityR) && (stabilityR > stabilityMinR));

            isPlace = (pVal < pThresh) && (longestField >= minFieldBins);

            reasonStr = 'passes';
            if ~isPlace
                if ~(pVal < pThresh); reasonStr = 'info_not_significant'; end
                if longestField < minFieldBins; reasonStr = 'small_field'; end
                % if isnan(stabilityR) || (stabilityR <= stabilityMinR); reasonStr = 'low_stability'; end
            end

            place_cell_results(pid).clusters(c).group.(group) = struct( ...
                'isPlace', isPlace, 'reason', reasonStr, 'info', infoObs, 'pVal', pVal, ...
                'peakRate', peakRate, 'fieldBins', longestField, 'stabilityR', stabilityR);
        end % group loop
    end % cluster loop

    % --- Plot average running and occupancy profiles for each group (unchanged) ---
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
    grid on;
    
    % Add figure title
    FigureTitle(fig_profiles, sprintf('Average Profiles - %s', probe_ids{pid}));
    ctl.figs.save_fig_to_join();

    % Determine session types based on actual trial data (unchanged)
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
    
    % Define colors for different sessions (unchanged)
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

    % For each cluster, create a single figure with both group firing rates overlaid (unchanged)
    fprintf('  Creating individual cluster plots...\n');
    
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
            rate_smooth = all_rate_smooth_by_group.(group)(c, :);
            rate_unsmooth = all_rate_unsmooth_by_group.(group)(c, :);
            
            % Skip if all NaN
            if isempty(rate_smooth) || all(isnan(rate_smooth))
                continue;
            end
            
            % Plot unsmoothed rate (lighter)
            plot(bin_centers, rate_unsmooth, '-', 'LineWidth', 1, 'Color', [group_colors.(group), 0.5], 'DisplayName', sprintf('%s (unsmoothed)', group_labels{g}));
            % Plot smoothed rate (solid line)
            plot(bin_centers, rate_smooth, '-', 'LineWidth', 2, 'Color', group_colors.(group), 'DisplayName', group_labels{g});
            
            % Annotate place cell result on the plot
            res = place_cell_results(pid).clusters(c).group.(group);
            if isfield(res, 'isPlace') && res.isPlace
                % mark peak position
                [~, peak_idx] = max(rate_smooth);
                pkx = bin_centers(peak_idx);
                pky = rate_smooth(peak_idx);
                plot(pkx, pky, 'v', 'MarkerSize', 8, 'MarkerFaceColor', group_colors.(group), 'MarkerEdgeColor', 'k');
                text(pkx, pky, sprintf('  place'), 'Color', group_colors.(group), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
            end
        end
        
        hold off;
        xlabel('Position (cm)');
        ylabel('Firing rate (Hz)');
        legend('show');
        grid on;
        
        % Add figure title
        FigureTitle(fig_cluster, sprintf('Cluster %d - %s', cluster.id, probe_ids{pid}));
        ctl.figs.save_fig_to_join();
    end

    % Plot heatmaps for each group in a single figure with four subplots (unchanged)
    fprintf('  Creating heatmap plots...\n');
    fig_heatmaps = ctl.figs.a4figure('landscape');

    % Sort clusters by maximum peak positional firing rate for long trials
    long_rate_mat = all_rate_smooth_by_group.('long');
    peak_positions = zeros(n_clusters, 1);
    for c = 1:n_clusters
        rate_profile = long_rate_mat(c, :);
        bin_centers = all_bin_centers_by_group.('long');
        if isempty(rate_profile) || all(isnan(rate_profile))
            peak_positions(c) = Inf;
            continue;
        end
        [~, max_idx] = max(rate_profile);
        peak_positions(c) = bin_centers(max_idx);
    end
    [~, sort_idx] = sort(peak_positions);
    sorted_cluster_ids = cluster_ids(sort_idx);

    for g = 1:2
        group = group_names{g};
        rate_mat = all_rate_smooth_by_group.(group);
        bin_centers = all_bin_centers_by_group.(group);
        if isempty(rate_mat)
            continue;
        end
        sorted_rate_mat = rate_mat(sort_idx, :);
        
        % Raw heatmap
        subplot(2, 2, g, 'Parent', fig_heatmaps);
        imagesc(bin_centers, 1:n_clusters, sorted_rate_mat);
        colormap('jet');
        colorbar;
        xlabel('Position (cm)');
        ylabel('Cluster (sorted by maximum peak position)');
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
        colormap('jet');
        colorbar;
        xlabel('Position (cm)');
        ylabel('Cluster (sorted by maximum peak position)');
        set(gca, 'YTick', 1:n_clusters, 'YTickLabel', sorted_cluster_ids);
        
        if strcmp(group, 'short')
            xlim([60, 120]);
        else
            xlim([0, 120]);
        end
    end
    
    FigureTitle(fig_heatmaps, sprintf('Heatmaps - %s', probe_ids{pid}));
    ctl.figs.save_fig_to_join();
    
    % Join all plots for this probe into one PDF
    fprintf('  Saving PDF for probe %s...\n', probe_ids{pid});
    fname = sprintf('%s.pdf', probe_ids{pid});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
    
    %save to folder
    save_results_dir = fullfile(output_base_dir, 'place_cell_results');
    if ~exist(save_results_dir, 'dir')
        mkdir(save_results_dir);
    end
    
    save_fname = fullfile(save_results_dir, sprintf('%s_place_cells.mat', probe_ids{pid}));
    save(save_fname, 'place_cell_results');
    
    fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
      % === X-normalized comparison (fraction of path) for this probe ===
    % This should be placed INSIDE the probe loop, AFTER you have computed:
    %   all_rate_smooth_by_group.long / .short
    %   all_bin_centers_by_group.long / .short
    %   cluster_ids, clusters

    % Make output directory for these PDFs
    xnorm_dir = fullfile(output_base_dir, 'xnorm_shape_comparison');
    if ~exist(xnorm_dir, 'dir')
        mkdir(xnorm_dir);
    end

    % Number of clusters (rows of the long rate matrix)
    n_clusters = size(all_rate_smooth_by_group.long, 1);

    % Output file name for this probe
    out_pdf = fullfile(xnorm_dir, sprintf('%s_shape_comparison_xnorm.pdf', probe_ids{pid}));

    % If file already exists, delete so we can recreate it
    if exist(out_pdf, 'file')
        delete(out_pdf);
    end

    fprintf('Saving x-normalized shape comparison for %d clusters to %s\n', ...
            n_clusters, out_pdf);

    % Loop over clusters
    for c = 1:n_clusters
        cluster_id = cluster_ids(c);

        % Extract x and y for this cluster
        pos_long   = all_bin_centers_by_group.long(:);          % 0–120 (physical)
        rate_long  = all_rate_smooth_by_group.long(c,:).';      % row -> column

        pos_short  = all_bin_centers_by_group.short(:);         % 60–120 (shifted)
        rate_short = all_rate_smooth_by_group.short(c,:).';     % row -> column

        % Skip if this cluster has no data in one of the groups
        if all(isnan(rate_long)) || all(isnan(rate_short))
            fprintf('Cluster %d: skipped (no data in one group)\n', cluster_id);
            continue;
        end

        % Normalized position (0 = start, 1 = goal) for each condition
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

        % Common normalized grid
        x_norm = linspace(0,1,100);

        % Interpolate both curves onto the same normalized grid
        rate_long_i  = interp1(pos_long_n,  rate_long_n,  x_norm, 'linear', 'extrap');
        rate_short_i = interp1(pos_short_n, rate_short_n, x_norm, 'linear', 'extrap');

        % Correlation of raw firing rates across normalized position
        r = corr(rate_long_i(:), rate_short_i(:), 'rows', 'complete');

        % Plot and append to PDF
        fig = figure('Visible','off');
        hold on;
        plot(x_norm, rate_long_i,  'b', 'LineWidth', 2);
        plot(x_norm, rate_short_i, 'r', 'LineWidth', 2);
        hold off;

        xlabel('Normalized track position (start \rightarrow goal)');
        ylabel('Firing rate (Hz)');   % raw firing rate
        legend({'Long (0–120 cm)','Short (60–120 cm)'}, 'Location', 'best');
        title(sprintf('Cluster %d – x-normalized comparison (r = %.2f)', cluster_id, r));
        grid on;

        exportgraphics(fig, out_pdf, 'Append', true);
        close(fig);

        fprintf('Cluster %d: added to PDF (r = %.2f)\n', cluster_id, r);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions (defined after main script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = skaggs_info(rateMap, occupancy)
    % compute Skaggs spatial information (bits/spike)
    % rateMap: vector of rates (Hz)
    % occupancy: vector of occupancy times (s) for same bins
    valid = ~isnan(rateMap) & ~isnan(occupancy) & (occupancy > 0);
    if sum(valid) == 0
        info = NaN; return;
    end
    p_i = occupancy(valid) ./ nansum(occupancy(valid));
    lam_i = rateMap(valid);
    lam = nansum(lam_i .* p_i);
    if lam == 0
        info = 0; return;
    end
    mask = (lam_i > 0) & (p_i > 0);
    info = nansum( p_i(mask) .* (lam_i(mask) ./ lam) .* log2(lam_i(mask) ./ lam) );
end

function rm = compute_rate_from_spikes(spk_times, pos_time_vec, pos_vec, edges, kernel, pad_size, smoothOcc)
    % compute smoothed rate map for a vector of spike times using the
    % provided position time series (pos_time_vec,pos_vec) and precomputed kernel.
    if isempty(spk_times)
        rm = nan(1, length(edges)-1); return;
    end
    % map spike times to positions
    spkPos = interp1(pos_time_vec, pos_vec, spk_times, 'linear', 'extrap');
    spike_count = histcounts(spkPos, edges);
    % pad counts, smooth by conv with kernel, unpad
    counts_padded = [fliplr(spike_count(1:pad_size)), spike_count, fliplr(spike_count(end-pad_size+1:end))];
    smoothCounts_padded = conv(counts_padded, kernel, 'same');
    smoothCounts = smoothCounts_padded(pad_size + 1 : end - pad_size);
    rm = smoothCounts ./ smoothOcc;
    % replace NaNs by zeros for correlation compatibility
    rm(isnan(rm)) = 0;
end
