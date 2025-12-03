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

% Statistical analysis parameters (set run_statistics=false to skip slow shuffle tests)
%
% IMPORTANT: If you change any of the statistical parameters below (minSpikes, minPeakRate,
% fieldFrac, minFieldBins, maxNumFields, nShuf, pThresh), you MUST set re_run_analysis=true
% to recompute the cached data. The cache does not track parameter changes, so old cached
% statistics will be invalid if you change classification criteria.
%
run_statistics = true;      % Set to false to skip statistical analysis and speed up processing
minSpikes = 25;             % minimum spikes per cluster per group to attempt classification
minPeakRate = 1.0;          % Hz minimum peak for field consideration
fieldFrac = 0.7;            % field threshold fraction of peak
minFieldBins = 5;           % contiguous bins for a field
maxNumFields = 1;           % maximum number of fields for uniqueness criterion
nShuf = 100;                % number of circular-shift shuffles
pThresh = 0.05;             % significance threshold for Skaggs shuffle
use_parallel = true;        % Use parallel processing for shuffles (requires Parallel Computing Toolbox)


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
        loaded_vars = load(cache_filepath);
        
        % Load core variables
        all_rate_smooth_by_group = loaded_vars.all_rate_smooth_by_group;
        all_bin_centers_by_group = loaded_vars.all_bin_centers_by_group;
        all_avg_velocity_by_group = loaded_vars.all_avg_velocity_by_group;
        all_Q1_rate_smooth_by_group = loaded_vars.all_Q1_rate_smooth_by_group;
        all_Q2_rate_smooth_by_group = loaded_vars.all_Q2_rate_smooth_by_group;
        all_Q3_rate_smooth_by_group = loaded_vars.all_Q3_rate_smooth_by_group;
        all_occ_by_group = loaded_vars.all_occ_by_group;
        all_rate_per_trial_by_group = loaded_vars.all_rate_per_trial_by_group;
        all_rate_per_trial_smooth_by_group = loaded_vars.all_rate_per_trial_smooth_by_group;
        all_spike_positions_by_group = loaded_vars.all_spike_positions_by_group;
        all_global_trial_indices_by_group = loaded_vars.all_global_trial_indices_by_group;
        all_rate_pooled_by_group = loaded_vars.all_rate_pooled_by_group;
        cluster_ids = loaded_vars.cluster_ids;
        n_clusters = loaded_vars.n_clusters;
        group_names = loaded_vars.group_names;
        group_labels = loaded_vars.group_labels;
        
        % Load statistical results if available (for backward compatibility)
        if isfield(loaded_vars, 'spatial_tuning_stats')
            spatial_tuning_stats = loaded_vars.spatial_tuning_stats;
            fprintf('  Loaded cached data for %d clusters (including statistical analysis)\n', n_clusters);
        elseif isfield(loaded_vars, 'place_cell_stats')
            % Legacy support for old cache files
            spatial_tuning_stats = loaded_vars.place_cell_stats;
            fprintf('  Loaded cached data for %d clusters (including statistical analysis - legacy format)\n', n_clusters);
        else
            fprintf('  Loaded cached data for %d clusters (no statistical analysis found - use re_run_analysis=true)\n', n_clusters);
        end
        
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
    all_bin_centers_by_group = struct();
    all_avg_velocity_by_group = struct();
    all_Q1_rate_smooth_by_group = struct();
    all_Q2_rate_smooth_by_group = struct();
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
                    endlll
                    
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
            
            % --- Pre-compute trial data for statistical testing (if enabled) ---
            trial_data_cache = cell(n_trials, 1);  % Cache trial data to avoid reloading
            
            if run_statistics
                for k = 1:n_trials
                    trial = trials_struct(k).trial;
                    motion_mask = trial.motion_mask();
                    pos = trial.position(motion_mask);
                    tvec = trial.probe_t(motion_mask);
                    
                    if ~isempty(pos) && ~isempty(tvec)
                        % Get spike times within this trial
                        st = cluster.spike_times;
                        spike_mask = st >= tvec(1) & st <= tvec(end);
                        st_in_trial = st(spike_mask);
                        
                        % Cache trial-specific data for shuffle test
                        trial_data_cache{k} = struct( ...
                            'pos', pos, 'tvec', tvec, 'st_in_trial', st_in_trial, ...
                            'duration', tvec(end) - tvec(1), ...
                            'occ_smooth', occ_per_trial_smooth(k, :));
                    end
                end
            end
            
            % --- Initialize storage if needed ---
            if ~isfield(all_rate_smooth_by_group, group)
                all_rate_smooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q1_rate_smooth_by_group.(group) = nan(n_clusters, n_bins);
                all_Q2_rate_smooth_by_group.(group) = nan(n_clusters, n_bins);
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
            
            % --- Perform statistical analysis on median firing rate ---
            if run_statistics
                % Count total spikes across all trials in this group
                n_spikes_group = 0;
                for k = 1:n_trials
                    if ~isempty(trial_data_cache{k})
                        n_spikes_group = n_spikes_group + length(trial_data_cache{k}.st_in_trial);
                    end
                end
                
                % Check minimum spike count
                if n_spikes_group < minSpikes
                    spatial_tuning_stats.(sprintf('cluster_%d', cluster.id)).(group) = struct( ...
                        'isSpatiallyTuned', false, 'reason', 'insufficient_spikes', ...
                        'n_spikes', n_spikes_group, 'info', NaN, 'pVal', NaN, ...
                        'peakRate', NaN, 'peakPosition', NaN, 'fieldBins', 0, 'numFields', 0);
                    continue;
                end
                
                % Use median of smoothed trials as the rate map for statistical testing
                median_rate_smooth = Q2_from_smoothed;
                median_occ_smooth = nanmean(occ_per_trial_smooth, 1);
                
                % Compute Skaggs information on median rate
                infoObs = skaggs_info(median_rate_smooth, median_occ_smooth);
                
                % Circular shuffle test (with parallel processing)
                infoNull = zeros(nShuf, 1);
                
                % Pre-generate random shifts for all shuffles (vectorized)
                shift_offsets = rand(nShuf, n_trials);
                
                % Choose acceleration method: Parallel > Serial
                if use_parallel && nShuf > 20
                    % Use parfor for parallel CPU processing
                    parfor s = 1:nShuf
                        rate_per_trial_shuf = compute_shuffled_rates(trial_data_cache, shift_offsets(s, :), ...
                                                                     edges, bin_size_cm, gauss_sigma_cm, n_bins);
                        median_rate_shuf = nanmedian(rate_per_trial_shuf, 1);
                        infoNull(s) = skaggs_info(median_rate_shuf, median_occ_smooth);
                    end
                else
                    % Use regular serial for loop
                    for s = 1:nShuf
                        rate_per_trial_shuf = compute_shuffled_rates(trial_data_cache, shift_offsets(s, :), ...
                                                                     edges, bin_size_cm, gauss_sigma_cm, n_bins);
                        median_rate_shuf = nanmedian(rate_per_trial_shuf, 1);
                        infoNull(s) = skaggs_info(median_rate_shuf, median_occ_smooth);
                    end
                end
                
                % Compute p-value
                pVal = (sum(infoNull >= infoObs) + 1) / (nShuf + 1);
            else
                % Skip statistical tests
                n_spikes_group = NaN;
                infoObs = NaN;
                pVal = NaN;
            end
            
            % Peak detection and field analysis on median rate
            if run_statistics
                median_rate_smooth = Q2_from_smoothed;
                [peakRate, peakIdx] = max(median_rate_smooth);
                peakPosition = bin_centers(peakIdx);
        
                % Field detection: count how many spatial bins belong to firing fields
                % 
                % FIELD DETECTION ALGORITHM:
                % 1. Set threshold = fieldFrac * peakRate (e.g., 0.2 * peak = 20% of max)
                % 2. Find all bins where firing rate >= threshold
                % 3. Group consecutive bins into "fields" (contiguous regions above threshold)
                % 4. Count the number of bins in each field
                % 5. Record the longest field and total number of separate fields
                %
                % EXAMPLE: If bins [5,6,7,8] and [15,16] are above threshold:
                %   - numFields = 2 (two separate regions)
                %   - fieldBins for field 1 = 4 bins (covering 8 cm with 2cm bins)
                %   - fieldBins for field 2 = 2 bins (covering 4 cm with 2cm bins)
                %   - longestField = 4 bins (the larger of the two)
                %
                % For spatial tuning classification, we require:
                %   - longestField >= minFieldBins (e.g., >= 5 bins = 10 cm)
                %   - numFields <= maxNumFields (typically 1, for single-field criterion)
                
                thresh = fieldFrac * peakRate;
                above = find(median_rate_smooth >= thresh);
                
                if isempty(above)
                    longestField = 0;
                    numFields = 0;
                    fieldStart = NaN;
                    fieldEnd = NaN;
                else
                    % Find contiguous segments by looking for breaks in the sequence
                    % diff(above) > 1 indicates a gap between consecutive bins
                    d = diff(above);
                    breaks = [0, find(d > 1), length(above)];
                    numFields = length(breaks) - 1;
                    segLens = zeros(numFields, 1);
                    segStarts = zeros(numFields, 1);
                    segEnds = zeros(numFields, 1);
                    for k = 1:numFields
                        idx = (breaks(k)+1) : breaks(k+1);
                        seg = above(idx);
                        segLens(k) = length(seg);
                        segStarts(k) = seg(1);  % First bin index in this field
                        segEnds(k) = seg(end);  % Last bin index in this field
                    end
                    [longestField, longestIdx] = max(segLens);
                    % Store the boundaries of the longest field (in bin indices)
                    fieldStartIdx = segStarts(longestIdx);
                    fieldEndIdx = segEnds(longestIdx);
                    % Convert to position in cm (bin edges, not centers)
                    fieldStart = bin_centers(fieldStartIdx) - bin_size_cm/2;
                    fieldEnd = bin_centers(fieldEndIdx) + bin_size_cm/2;
                end
                
                % Classification logic
                isSpatiallyTuned = (pVal < pThresh) && (longestField >= minFieldBins) && ...
                          (numFields <= maxNumFields) && (peakRate >= minPeakRate);
                
                reasonStr = 'passes';
                if ~isSpatiallyTuned
                    if pVal >= pThresh
                        reasonStr = 'info_not_significant';
                    elseif peakRate < minPeakRate
                        reasonStr = 'low_peak';
                    elseif longestField < minFieldBins
                        reasonStr = 'small_field';
                    elseif numFields > maxNumFields
                        reasonStr = 'multiple_fields';
                    end
                end
                
                % Store statistical results (including field boundaries for plotting)
                spatial_tuning_stats.(sprintf('cluster_%d', cluster.id)).(group) = struct( ...
                    'isSpatiallyTuned', isSpatiallyTuned, 'reason', reasonStr, ...
                    'n_spikes', n_spikes_group, 'info', infoObs, 'pVal', pVal, ...
                    'peakRate', peakRate, 'peakPosition', peakPosition, ...
                    'fieldBins', longestField, 'numFields', numFields, ...
                    'fieldStart', fieldStart, 'fieldEnd', fieldEnd);
                
                fprintf('      Spatial tuning analysis: isSpatiallyTuned=%d, reason=%s, info=%.3f, pVal=%.3f, peak=%.2f Hz at %.1f cm\n', ...
                        isSpatiallyTuned, reasonStr, infoObs, pVal, peakRate, peakPosition);
            end
        end
    end
    
        % --- Save computed data to cache ---
        fprintf('  Saving analysis data to cache: %s\n', cache_filename);
        save(cache_filepath, 'all_rate_smooth_by_group', ...
             'all_bin_centers_by_group', 'all_avg_velocity_by_group', ...
             'all_Q1_rate_smooth_by_group', 'all_Q2_rate_smooth_by_group', ...
             'all_Q3_rate_smooth_by_group', 'all_occ_by_group', ...
             'all_rate_per_trial_by_group', 'all_rate_per_trial_smooth_by_group', ...
             'all_spike_positions_by_group', 'all_global_trial_indices_by_group', ...
             'all_rate_pooled_by_group', ...
             'spatial_tuning_stats', ...
             'cluster_ids', 'n_clusters', 'group_names', 'group_labels', ...
             'bin_size_cm', 'gauss_sigma_cm', '-v7.3');
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
            
            % Get statistics for this cluster (if available)
            cluster_stats = [];
            if exist('spatial_tuning_stats', 'var')
                cluster_field = sprintf('cluster_%d', cluster.id);
                if isfield(spatial_tuning_stats, cluster_field)
                    cluster_stats = spatial_tuning_stats.(cluster_field);
                end
            end
            
            % Create combined figure using helper function
            fig_cluster = plot_cluster_spatial_profile(cluster.id, all_bin_centers_by_group, ...
                                                       rate_data, group_names, group_labels, ...
                                                       group_colors, probe_ids{pid}, cluster_stats);
            
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
            rate_mat = all_Q2_rate_smooth_by_group.(group);  % Use median (Q2) instead of mean
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
            
            % Normalized heatmap (min-max normalization: (x - min) / (max - min))
            min_vals = min(sorted_rate_mat, [], 2);
            max_vals = max(sorted_rate_mat, [], 2);
            norm_rate_mat = (sorted_rate_mat - min_vals) ./ (max_vals - min_vals);
            norm_rate_mat(isnan(norm_rate_mat) | isinf(norm_rate_mat)) = 0;
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
        
        % --- Plot spatially tuned clusters with Gaussian fits ---
        if exist('spatial_tuning_stats', 'var')
            fprintf('  Creating spatially tuned clusters plot with Gaussian fits...\n');
            
            % Identify clusters that are spatially tuned in BOTH long AND short trials
            spatially_tuned_both = [];
            
            for c = 1:n_clusters
                cluster_field = sprintf('cluster_%d', cluster_ids(c));
                if isfield(spatial_tuning_stats, cluster_field)
                    stats = spatial_tuning_stats.(cluster_field);
                    % Must be tuned in BOTH groups
                    if isfield(stats, 'long') && stats.long.isSpatiallyTuned && ...
                       isfield(stats, 'short') && stats.short.isSpatiallyTuned
                        spatially_tuned_both = [spatially_tuned_both; c];
                    end
                end
            end
            
            n_tuned_both = length(spatially_tuned_both);
            fprintf('    Found %d clusters spatially tuned in BOTH long and short trials\n', n_tuned_both);
            
            if n_tuned_both > 0
                % Sort clusters by peak position in long trials
                peak_pos_long = zeros(n_tuned_both, 1);
                for i = 1:n_tuned_both
                    c = spatially_tuned_both(i);
                    cluster_field = sprintf('cluster_%d', cluster_ids(c));
                    peak_pos_long(i) = spatial_tuning_stats.(cluster_field).long.peakPosition;
                end
                [~, sort_idx] = sort(peak_pos_long);
                spatially_tuned_both_sorted = spatially_tuned_both(sort_idx);
                
                % Generate pastel colors for all clusters
                pastel_colors = generate_pastel_colors(n_tuned_both);
                
                % --- Create 3D plot with clusters along z-axis ---
                fprintf('  Creating 3D visualization of spatially tuned clusters...\n');
                fig_3d = ctl.figs.a4figure('landscape');
                
                % Left subplot: Long trials in 3D
                subplot(1, 2, 1, 'Parent', fig_3d);
                hold on;
                
                for i = 1:n_tuned_both
                    c = spatially_tuned_both_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    % Get median firing rate
                    bin_centers_long = all_bin_centers_by_group.long;
                    median_rate_long = all_Q2_rate_smooth_by_group.long(c, :);
                    
                    % Fit Gaussian
                    peak_pos = spatial_tuning_stats.(cluster_field).long.peakPosition;
                    peak_rate = spatial_tuning_stats.(cluster_field).long.peakRate;
                    sigma_guess = gauss_sigma_cm;
                    p0 = [peak_rate, peak_pos, sigma_guess];
                    
                    % Create z-coordinate (cluster depth) - convert to double
                    z_position = double(cluster_id);
                    
                    % Fit and plot in 3D
                    try
                        gauss_fit = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2 * p(3)^2));
                        options = optimset('Display', 'off');
                        p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_long, median_rate_long, ...
                                           [0, 0, 1], [inf, 120, 50], options);
                        
                        x_fit = linspace(0, 120, 200);
                        y_fit = gauss_fit(p_fit, x_fit);
                        z_fit = ones(size(x_fit)) * z_position;
                        
                        % Plot filled area using patch (y=cluster, z=firing rate)
                        x_patch = [x_fit, fliplr(x_fit)];
                        y_patch = [z_fit, fliplr(z_fit)];  % cluster ID on y-axis
                        z_patch = [y_fit, zeros(size(y_fit))];  % firing rate on z-axis
                        patch(x_patch, y_patch, z_patch, color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                        
                        % Plot line on top
                        plot3(x_fit, z_fit, y_fit, 'Color', color, 'LineWidth', 2);
                    catch
                        % If fit fails, just plot the median
                        y_data = ones(size(bin_centers_long)) * z_position;
                        plot3(bin_centers_long, y_data, median_rate_long, 'Color', color, 'LineWidth', 2);
                    end
                end
                
                xlabel('Position (cm)');
                ylabel('Cluster ID');
                zlabel('Firing Rate (Hz)');
                title('Long Trials (0-120 cm)');
                grid on;
                view(3);
                hold off;
                
                % Right subplot: Short trials in 3D
                subplot(1, 2, 2, 'Parent', fig_3d);
                hold on;
                
                for i = 1:n_tuned_both
                    c = spatially_tuned_both_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    % Get median firing rate
                    bin_centers_short = all_bin_centers_by_group.short;
                    median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
                    
                    % Fit Gaussian
                    peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
                    peak_rate = spatial_tuning_stats.(cluster_field).short.peakRate;
                    p0 = [peak_rate, peak_pos, sigma_guess];
                    
                    % Create z-coordinate (cluster depth) - convert to double
                    z_position = double(cluster_id);
                    
                    % Fit and plot in 3D
                    try
                        gauss_fit = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2 * p(3)^2));
                        options = optimset('Display', 'off');
                        p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_short, ...
                                           [0, 60, 1], [inf, 120, 50], options);
                        
                        x_fit = linspace(60, 120, 200);
                        y_fit = gauss_fit(p_fit, x_fit);
                        z_fit = ones(size(x_fit)) * z_position;
                        
                        % Plot filled area using patch (y=cluster, z=firing rate)
                        x_patch = [x_fit, fliplr(x_fit)];
                        y_patch = [z_fit, fliplr(z_fit)];  % cluster ID on y-axis
                        z_patch = [y_fit, zeros(size(y_fit))];  % firing rate on z-axis
                        patch(x_patch, y_patch, z_patch, color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                        
                        % Plot line on top
                        plot3(x_fit, z_fit, y_fit, 'Color', color, 'LineWidth', 2);
                    catch
                        % If fit fails, just plot the median
                        y_data = ones(size(bin_centers_short)) * z_position;
                        plot3(bin_centers_short, y_data, median_rate_short, 'Color', color, 'LineWidth', 2);
                    end
                end
                
                xlabel('Position (cm)');
                ylabel('Cluster ID');
                zlabel('Firing Rate (Hz)');
                title('Short Trials (60-120 cm)');
                grid on;
                view(3);
                hold off;
                
                FigureTitle(fig_3d, sprintf('3D Spatial Tuning - %s', probe_ids{pid}));
                ctl.figs.save_fig_to_join();
            else
                fprintf('    No clusters found that are spatially tuned in both groups. Skipping plot.\n');
            end
        end
        
        % Save statistics text summary
        %
        % NOTE: Only a text summary is saved for reference. The actual caching 
        % (including all firing rate data and statistics) is done with 
        % *_spatial_analysis_cache.mat
        %
        if run_statistics && exist('spatial_tuning_stats', 'var')
            % Prepare analysis parameters structure
            analysis_params = struct( ...
                'bin_size_cm', bin_size_cm, ...
                'gauss_sigma_cm', gauss_sigma_cm, ...
                'minSpikes', minSpikes, ...
                'minPeakRate', minPeakRate, ...
                'fieldFrac', fieldFrac, ...
                'minFieldBins', minFieldBins, ...
                'maxNumFields', maxNumFields, ...
                'nShuf', nShuf, ...
                'pThresh', pThresh);
            
            % Call external function to save statistics
            save_spatial_statistics(spatial_tuning_stats, probe_ids{pid}, ctl.figs.curr_dir, ...
                                   analysis_params, cluster_ids, group_names, group_labels);
        end
        
        % Join all plots for this probe into one PDF
        fprintf('  Saving PDF for probe %s...\n', probe_ids{pid});
        fname = sprintf('%s.pdf', probe_ids{pid});
        ctl.figs.join_figs(fname, overwrite);
        ctl.figs.clear_figs();
        
        fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
    end
end

fprintf('\nSpatial firing rate profile analysis completed for %d probe(s)\n', length(probe_ids));

% --- Summary of spatial tuning classification ---
if exist('spatial_tuning_stats', 'var')
    fprintf('\n=== SPATIAL TUNING CLASSIFICATION SUMMARY ===\n');
    cluster_names = fieldnames(spatial_tuning_stats);
    for c = 1:length(cluster_names)
        cluster_name = cluster_names{c};
        fprintf('\n%s:\n', cluster_name);
        stats = spatial_tuning_stats.(cluster_name);
        group_fields = fieldnames(stats);
        for g = 1:length(group_fields)
            group = group_fields{g};
            s = stats.(group);
            fprintf('  %s: isSpatiallyTuned=%d, reason=%s, spikes=%d, info=%.3f, pVal=%.3f, peak=%.2f Hz @ %.1f cm, fields=%d\n', ...
                    group, s.isSpatiallyTuned, s.reason, s.n_spikes, s.info, s.pVal, s.peakRate, s.peakPosition, s.numFields);
        end
    end
end

% --- Helper function to generate pastel colors ---
function colors = generate_pastel_colors(n)
    % Generate n distinct pastel colors using HSV color space
    % Pastel colors have high value (brightness) and low-to-medium saturation
    hues = linspace(0, 1, n+1);
    hues = hues(1:end-1); % Remove duplicate at end
    colors = zeros(n, 3);
    for i = 1:n
        % Use moderate saturation (0.4-0.6) and high value (0.85-0.95) for pastel effect
        hsv_color = [hues(i), 0.5, 0.9];
        colors(i, :) = hsv2rgb(hsv_color);
    end
end