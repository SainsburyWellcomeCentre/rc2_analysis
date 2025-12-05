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
% STATISTICAL TESTING:
%   Spatial tuning significance is tested using a bin-shuffle approach analogous
%   to velocity tuning in ShuffleTuning.m. The rate matrix (n_trials × n_bins) is
%   completely shuffled to break trial-bin associations, and Skaggs information
%   is computed on the shuffled data to create a null distribution.
%
% Requires functions from lib/spatial_analysis/:
%   - compute_trial_firing_rate.m
%   - smooth_spatial_rate.m
%   - compute_quartiles_and_smooth.m
%   - collect_and_group_trials.m
%   - compute_shuffled_spatial_info.m

%%
% Configuration
experiment_groups       = {'ambient_light'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'spatial_firing_rate', 'ambient_light'};
plot_single_cluster_fig = true;
plot_heatmp_cluster_fig = true;
re_run_analysis         = false;  % Set to true to recompute all metrics, false to load cached data

% Parallel processing configuration
use_parallel            = true;   % Set to false to disable parallel processing entirely
max_workers             = 4;      % Maximum number of parallel workers (reduce if running out of memory)
                                  % Recommended: 2-4 for laptops, 4-8 for desktops with 16GB+ RAM


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
fieldFrac = 0.7;            % field threshold fraction of rate range (min to max)
minFieldBins = 5;           % contiguous bins for a field
maxNumFields = 1;           % maximum number of fields for uniqueness criterion
nShuf = 100;                % number of bin-shuffles for significance testing
pThresh = 0.05;             % significance threshold for Skaggs shuffle


% Initialize controller and get probe IDs
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));


for pid = 5:length(probe_ids)
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
        
        % Still need to load data object for clusters info and recreate trial_groups
        data = ctl.load_formatted_data(probe_ids{pid});
        clusters = data.selected_clusters();
        sessions = data.motion_sessions();
        
        % Recreate trial_groups for plotting (needed even when loading from cache)
        [trial_groups, ~, ~] = collect_and_group_trials(sessions, position_threshold_cm);
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

        % --- Pre-calculate bin information for both groups ---
        % Long trials
        edges_long = 0:bin_size_cm:120;
        bin_centers_long = edges_long(1:end-1) + bin_size_cm/2;
        n_bins_long = length(bin_centers_long);
        plot_bin_centers_long = bin_centers_long;
        
        % Short trials
        edges_short = 0:bin_size_cm:60;
        bin_centers_short = edges_short(1:end-1) + bin_size_cm/2;
        n_bins_short = length(bin_centers_short);
        plot_bin_centers_short = bin_centers_short + 60;
        
        % Pre-allocate plain arrays for parfor (not structures)
        rate_smooth_long = nan(n_clusters, n_bins_long);
        rate_smooth_short = nan(n_clusters, n_bins_short);
        Q1_rate_smooth_long = nan(n_clusters, n_bins_long);
        Q1_rate_smooth_short = nan(n_clusters, n_bins_short);
        Q2_rate_smooth_long = nan(n_clusters, n_bins_long);
        Q2_rate_smooth_short = nan(n_clusters, n_bins_short);
        Q3_rate_smooth_long = nan(n_clusters, n_bins_long);
        Q3_rate_smooth_short = nan(n_clusters, n_bins_short);
        avg_velocity_long = nan(n_clusters, n_bins_long);
        avg_velocity_short = nan(n_clusters, n_bins_short);
        occ_long = nan(n_clusters, n_bins_long);
        occ_short = nan(n_clusters, n_bins_short);
        rate_pooled_long = nan(n_clusters, n_bins_long);
        rate_pooled_short = nan(n_clusters, n_bins_short);
        rate_per_trial_long = cell(n_clusters, 1);
        rate_per_trial_short = cell(n_clusters, 1);
        rate_per_trial_smooth_long = cell(n_clusters, 1);
        rate_per_trial_smooth_short = cell(n_clusters, 1);
        spike_positions_long = cell(n_clusters, 1);
        spike_positions_short = cell(n_clusters, 1);
        global_trial_indices_long = cell(n_clusters, 1);
        global_trial_indices_short = cell(n_clusters, 1);
        
        % Pre-allocate arrays for statistical results (to avoid structure issues in parfor)
        if run_statistics
            stats_isSpatiallyTuned_long = false(n_clusters, 1);
            stats_isSpatiallyTuned_short = false(n_clusters, 1);
            stats_reason_long = cell(n_clusters, 1);
            stats_reason_short = cell(n_clusters, 1);
            stats_n_spikes_long = nan(n_clusters, 1);
            stats_n_spikes_short = nan(n_clusters, 1);
            stats_info_long = nan(n_clusters, 1);
            stats_info_short = nan(n_clusters, 1);
            stats_pVal_long = nan(n_clusters, 1);
            stats_pVal_short = nan(n_clusters, 1);
            stats_infoNull_long = cell(n_clusters, 1);
            stats_infoNull_short = cell(n_clusters, 1);
            stats_peakRate_long = nan(n_clusters, 1);
            stats_peakRate_short = nan(n_clusters, 1);
            stats_peakPosition_long = nan(n_clusters, 1);
            stats_peakPosition_short = nan(n_clusters, 1);
            stats_fieldBins_long = nan(n_clusters, 1);
            stats_fieldBins_short = nan(n_clusters, 1);
            stats_numFields_long = nan(n_clusters, 1);
            stats_numFields_short = nan(n_clusters, 1);
            stats_fieldStart_long = nan(n_clusters, 1);
            stats_fieldStart_short = nan(n_clusters, 1);
            stats_fieldEnd_long = nan(n_clusters, 1);
            stats_fieldEnd_short = nan(n_clusters, 1);
        end

        % --- Compute spatial firing rate profiles for all clusters ---
        if use_parallel
            % Initialize or get existing parallel pool with limited workers
            current_pool = gcp('nocreate');
            if isempty(current_pool)
                fprintf('  Starting parallel pool with %d workers...\n', max_workers);
                parpool('local', max_workers);
            elseif current_pool.NumWorkers ~= max_workers
                fprintf('  Restarting parallel pool with %d workers (was %d)...\n', max_workers, current_pool.NumWorkers);
                delete(current_pool);
                parpool('local', max_workers);
            else
                fprintf('  Using existing parallel pool with %d workers\n', current_pool.NumWorkers);
            end
            fprintf('  Computing spatial firing rate profiles for %d clusters (parallelized)...\n', n_clusters);
        else
            fprintf('  Computing spatial firing rate profiles for %d clusters (sequential - parallel disabled)...\n', n_clusters);
        end
        
        if use_parallel
            parfor c = 1:n_clusters
            cluster = clusters(c);
            
            for g = 1:2
                group = group_names{g};
                trials_struct = trial_groups.(group);
            
            % Use pre-calculated bin information
            if strcmp(group, 'long')
                edges = edges_long;
                bin_centers = bin_centers_long;
                n_bins = n_bins_long;
                plot_bin_centers = plot_bin_centers_long;
            else
                edges = edges_short;
                bin_centers = bin_centers_short;
                n_bins = n_bins_short;
                plot_bin_centers = plot_bin_centers_short;
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
            
            % Compute POOLED trial rate: sum all counts/occupancy, then smooth
            % This is for comparison with the _identification script approach
            pooled_count = nansum(count_per_trial, 1);
            pooled_occ = nansum(occ_per_trial, 1);
            pooled_count_smooth = smooth_spatial_rate(pooled_count, bin_size_cm, gauss_sigma_cm);
            pooled_occ_smooth = smooth_spatial_rate(pooled_occ, bin_size_cm, gauss_sigma_cm);
            rate_pooled = pooled_count_smooth ./ pooled_occ_smooth;
            rate_pooled(isnan(rate_pooled) | isinf(rate_pooled)) = 0;
            
            % --- Print statistics (commented out for parfor compatibility) ---
            % fprintf('\n    Cluster %d, %s trials (n=%d):\n', cluster.id, group, n_trials);
            % fprintf('      Smoothed mean rate: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
            %         nanmean(mean_rate_smooth), min(mean_rate_smooth), max(mean_rate_smooth));
            % fprintf('      Smoothed Q2 (median): mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
            %         nanmean(Q2_from_smoothed), min(Q2_from_smoothed), max(Q2_from_smoothed));
            % fprintf('      Pooled rate: mean=%.2f Hz, range=[%.2f, %.2f] Hz\n', ...
            %         nanmean(rate_pooled), min(rate_pooled), max(rate_pooled));
            
            % Compute mean occupancy and velocity across trials
            mean_occupancy = nanmean(occ_per_trial, 1);
            mean_velocity = nanmean(vel_per_trial, 1);
            
            % --- Store results in plain arrays (not structures, for parfor compatibility) ---
            if strcmp(group, 'long')
                rate_smooth_long(c, :) = mean_rate_smooth;
                Q1_rate_smooth_long(c, :) = Q1_from_smoothed;
                Q2_rate_smooth_long(c, :) = Q2_from_smoothed;
                Q3_rate_smooth_long(c, :) = Q3_from_smoothed;
                avg_velocity_long(c, :) = mean_velocity;
                occ_long(c, :) = mean_occupancy;
                rate_per_trial_long{c} = rate_per_trial;
                rate_per_trial_smooth_long{c} = rate_per_trial_smooth;
                spike_positions_long{c} = spike_positions_per_trial;
                global_trial_indices_long{c} = global_trial_indices;
                rate_pooled_long(c, :) = rate_pooled;
            else
                rate_smooth_short(c, :) = mean_rate_smooth;
                Q1_rate_smooth_short(c, :) = Q1_from_smoothed;
                Q2_rate_smooth_short(c, :) = Q2_from_smoothed;
                Q3_rate_smooth_short(c, :) = Q3_from_smoothed;
                avg_velocity_short(c, :) = mean_velocity;
                occ_short(c, :) = mean_occupancy;
                rate_per_trial_short{c} = rate_per_trial;
                rate_per_trial_smooth_short{c} = rate_per_trial_smooth;
                spike_positions_short{c} = spike_positions_per_trial;
                global_trial_indices_short{c} = global_trial_indices;
                rate_pooled_short(c, :) = rate_pooled;
            end
            
            % --- Perform statistical analysis on median firing rate ---
            if run_statistics
                % Count total spikes across all trials in this group
                n_spikes_group = nansum(count_per_trial(:));
                
                % Check minimum spike count
                if n_spikes_group < minSpikes
                    if strcmp(group, 'long')
                        stats_isSpatiallyTuned_long(c) = false;
                        stats_reason_long{c} = 'insufficient_spikes';
                        stats_n_spikes_long(c) = n_spikes_group;
                        stats_info_long(c) = NaN;
                        stats_pVal_long(c) = NaN;
                        stats_peakRate_long(c) = NaN;
                        stats_peakPosition_long(c) = NaN;
                        stats_fieldBins_long(c) = 0;
                        stats_numFields_long(c) = 0;
                        stats_fieldStart_long(c) = NaN;
                        stats_fieldEnd_long(c) = NaN;
                    else
                        stats_isSpatiallyTuned_short(c) = false;
                        stats_reason_short{c} = 'insufficient_spikes';
                        stats_n_spikes_short(c) = n_spikes_group;
                        stats_info_short(c) = NaN;
                        stats_pVal_short(c) = NaN;
                        stats_peakRate_short(c) = NaN;
                        stats_peakPosition_short(c) = NaN;
                        stats_fieldBins_short(c) = 0;
                        stats_numFields_short(c) = 0;
                        stats_fieldStart_short(c) = NaN;
                        stats_fieldEnd_short(c) = NaN;
                    end
                    continue;
                end
                
                % Use median of smoothed trials as the rate map for statistical testing
                median_rate_smooth = Q2_from_smoothed;
                median_occ_smooth = nanmean(occ_per_trial_smooth, 1);
                
                % Compute Skaggs information on median rate
                infoObs = skaggs_info(median_rate_smooth, median_occ_smooth);
                
                % Bin-shuffle test (analogous to velocity tuning in ShuffleTuning.m)
                % Shuffles all rate matrix values to break trial-bin associations
                infoNull = compute_shuffled_spatial_info(rate_per_trial_smooth, median_occ_smooth, nShuf);
                
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
                % 1. Set threshold = min_rate + fieldFrac * (max_rate - min_rate)
                %    This adapts to the modulation depth rather than using absolute peak
                %    Example: rates 10-20 Hz with fieldFrac=0.7 → threshold = 10 + 0.7*10 = 17 Hz
                %             rates 0-20 Hz with fieldFrac=0.7 → threshold = 0 + 0.7*20 = 14 Hz
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
                
                minRate = min(median_rate_smooth);
                maxRate = max(median_rate_smooth);
                thresh = minRate + fieldFrac * (maxRate - minRate);
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
                
                % Store statistical results in plain arrays (not structures, for parfor compatibility)
                if strcmp(group, 'long')
                    stats_isSpatiallyTuned_long(c) = isSpatiallyTuned;
                    stats_reason_long{c} = reasonStr;
                    stats_n_spikes_long(c) = n_spikes_group;
                    stats_info_long(c) = infoObs;
                    stats_pVal_long(c) = pVal;
                    stats_infoNull_long{c} = infoNull;
                    stats_peakRate_long(c) = peakRate;
                    stats_peakPosition_long(c) = peakPosition;
                    stats_fieldBins_long(c) = longestField;
                    stats_numFields_long(c) = numFields;
                    stats_fieldStart_long(c) = fieldStart;
                    stats_fieldEnd_long(c) = fieldEnd;
                else
                    stats_isSpatiallyTuned_short(c) = isSpatiallyTuned;
                    stats_reason_short{c} = reasonStr;
                    stats_n_spikes_short(c) = n_spikes_group;
                    stats_info_short(c) = infoObs;
                    stats_pVal_short(c) = pVal;
                    stats_infoNull_short{c} = infoNull;
                    stats_peakRate_short(c) = peakRate;
                    stats_peakPosition_short(c) = peakPosition;
                    stats_fieldBins_short(c) = longestField;
                    stats_numFields_short(c) = numFields;
                    stats_fieldStart_short(c) = fieldStart;
                    stats_fieldEnd_short(c) = fieldEnd;
                end
                
                % fprintf('      Spatial tuning analysis: isSpatiallyTuned=%d, reason=%s, info=%.3f, pVal=%.3f, peak=%.2f Hz at %.1f cm\n', ...
                %         isSpatiallyTuned, reasonStr, infoObs, pVal, peakRate, peakPosition);
            end
            end
        end
        else
            % Sequential processing (no parallel pool)
            for c = 1:n_clusters
            cluster = clusters(c);
            
            for g = 1:2
                group = group_names{g};
                trials_struct = trial_groups.(group);
            
            % Use pre-calculated bin information
            if strcmp(group, 'long')
                edges = edges_long;
                bin_centers = bin_centers_long;
                n_bins = n_bins_long;
                plot_bin_centers = plot_bin_centers_long;
            else
                edges = edges_short;
                bin_centers = bin_centers_short;
                n_bins = n_bins_short;
                plot_bin_centers = plot_bin_centers_short;
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
                        spike_pos = interp1(tvec, pos, st_in_trial, 'linear');
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
            
            % Compute POOLED trial rate: sum all counts/occupancy, then smooth
            % This is for comparison with the _identification script approach
            pooled_count = nansum(count_per_trial, 1);
            pooled_occ = nansum(occ_per_trial, 1);
            pooled_count_smooth = smooth_spatial_rate(pooled_count, bin_size_cm, gauss_sigma_cm);
            pooled_occ_smooth = smooth_spatial_rate(pooled_occ, bin_size_cm, gauss_sigma_cm);
            rate_pooled = pooled_count_smooth ./ pooled_occ_smooth;
            rate_pooled(isnan(rate_pooled) | isinf(rate_pooled)) = 0;
            
            % Compute mean occupancy and velocity across trials
            mean_occupancy = nanmean(occ_per_trial, 1);
            mean_velocity = nanmean(vel_per_trial, 1);
            
            % --- Store results in plain arrays (not structures, for compatibility) ---
            if strcmp(group, 'long')
                rate_smooth_long(c, :) = mean_rate_smooth;
                Q1_rate_smooth_long(c, :) = Q1_from_smoothed;
                Q2_rate_smooth_long(c, :) = Q2_from_smoothed;
                Q3_rate_smooth_long(c, :) = Q3_from_smoothed;
                avg_velocity_long(c, :) = mean_velocity;
                occ_long(c, :) = mean_occupancy;
                rate_per_trial_long{c} = rate_per_trial;
                rate_per_trial_smooth_long{c} = rate_per_trial_smooth;
                spike_positions_long{c} = spike_positions_per_trial;
                global_trial_indices_long{c} = global_trial_indices;
                rate_pooled_long(c, :) = rate_pooled;
            else
                rate_smooth_short(c, :) = mean_rate_smooth;
                Q1_rate_smooth_short(c, :) = Q1_from_smoothed;
                Q2_rate_smooth_short(c, :) = Q2_from_smoothed;
                Q3_rate_smooth_short(c, :) = Q3_from_smoothed;
                avg_velocity_short(c, :) = mean_velocity;
                occ_short(c, :) = mean_occupancy;
                rate_per_trial_short{c} = rate_per_trial;
                rate_per_trial_smooth_short{c} = rate_per_trial_smooth;
                spike_positions_short{c} = spike_positions_per_trial;
                global_trial_indices_short{c} = global_trial_indices;
                rate_pooled_short(c, :) = rate_pooled;
            end
            
            % --- Perform statistical analysis on median firing rate ---
            if run_statistics
                % Count total spikes across all trials in this group
                n_spikes_group = nansum(count_per_trial(:));
                
                % Check minimum spike count
                if n_spikes_group < minSpikes
                    if strcmp(group, 'long')
                        stats_isSpatiallyTuned_long(c) = false;
                        stats_reason_long{c} = 'insufficient_spikes';
                        stats_n_spikes_long(c) = n_spikes_group;
                        stats_info_long(c) = NaN;
                        stats_pVal_long(c) = NaN;
                        stats_infoNull_long{c} = [];
                        stats_peakRate_long(c) = NaN;
                        stats_peakPosition_long(c) = NaN;
                        stats_fieldBins_long(c) = NaN;
                        stats_numFields_long(c) = NaN;
                        stats_fieldStart_long(c) = NaN;
                        stats_fieldEnd_long(c) = NaN;
                    else
                        stats_isSpatiallyTuned_short(c) = false;
                        stats_reason_short{c} = 'insufficient_spikes';
                        stats_n_spikes_short(c) = n_spikes_group;
                        stats_info_short(c) = NaN;
                        stats_pVal_short(c) = NaN;
                        stats_infoNull_short{c} = [];
                        stats_peakRate_short(c) = NaN;
                        stats_peakPosition_short(c) = NaN;
                        stats_fieldBins_short(c) = NaN;
                        stats_numFields_short(c) = NaN;
                        stats_fieldStart_short(c) = NaN;
                        stats_fieldEnd_short(c) = NaN;
                    end
                    continue;
                end
                
                % Use median of smoothed trials as the rate map for statistical testing
                median_rate_smooth = Q2_from_smoothed;
                median_occ_smooth = nanmean(occ_per_trial_smooth, 1);
                
                % Compute Skaggs information on median rate
                infoObs = skaggs_info(median_rate_smooth, median_occ_smooth);
                
                % Bin-shuffle test (analogous to velocity tuning in ShuffleTuning.m)
                % Shuffles all rate matrix values to break trial-bin associations
                infoNull = compute_shuffled_spatial_info(rate_per_trial_smooth, median_occ_smooth, nShuf);
                
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
        
                minRate = min(median_rate_smooth);
                maxRate = max(median_rate_smooth);
                thresh = minRate + fieldFrac * (maxRate - minRate);
                above = find(median_rate_smooth >= thresh);
                
                if isempty(above)
                    longestField = 0;
                    numFields = 0;
                    fieldStart = NaN;
                    fieldEnd = NaN;
                else
                    d = diff(above);
                    breaks = [0, find(d > 1), length(above)];
                    numFields = length(breaks) - 1;
                    segLens = zeros(numFields, 1);
                    segStarts = zeros(numFields, 1);
                    segEnds = zeros(numFields, 1);
                    for k = 1:numFields
                        seg = above(breaks(k)+1:breaks(k+1));
                        segLens(k) = length(seg);
                        segStarts(k) = seg(1);
                        segEnds(k) = seg(end);
                    end
                    [longestField, longestIdx] = max(segLens);
                    fieldStartIdx = segStarts(longestIdx);
                    fieldEndIdx = segEnds(longestIdx);
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
                
                % Store statistical results
                if strcmp(group, 'long')
                    stats_isSpatiallyTuned_long(c) = isSpatiallyTuned;
                    stats_reason_long{c} = reasonStr;
                    stats_n_spikes_long(c) = n_spikes_group;
                    stats_info_long(c) = infoObs;
                    stats_pVal_long(c) = pVal;
                    stats_infoNull_long{c} = infoNull;
                    stats_peakRate_long(c) = peakRate;
                    stats_peakPosition_long(c) = peakPosition;
                    stats_fieldBins_long(c) = longestField;
                    stats_numFields_long(c) = numFields;
                    stats_fieldStart_long(c) = fieldStart;
                    stats_fieldEnd_long(c) = fieldEnd;
                else
                    stats_isSpatiallyTuned_short(c) = isSpatiallyTuned;
                    stats_reason_short{c} = reasonStr;
                    stats_n_spikes_short(c) = n_spikes_group;
                    stats_info_short(c) = infoObs;
                    stats_pVal_short(c) = pVal;
                    stats_infoNull_short{c} = infoNull;
                    stats_peakRate_short(c) = peakRate;
                    stats_peakPosition_short(c) = peakPosition;
                    stats_fieldBins_short(c) = longestField;
                    stats_numFields_short(c) = numFields;
                    stats_fieldStart_short(c) = fieldStart;
                    stats_fieldEnd_short(c) = fieldEnd;
                end
                
                fprintf('    Cluster %d/%d, %s: isSpatiallyTuned=%d, reason=%s\n', c, n_clusters, group, isSpatiallyTuned, reasonStr);
            end
            end
        end
        end
        
        % --- Reassemble structures from plain arrays after parfor ---
        all_rate_smooth_by_group.long = rate_smooth_long;
        all_rate_smooth_by_group.short = rate_smooth_short;
        all_Q1_rate_smooth_by_group.long = Q1_rate_smooth_long;
        all_Q1_rate_smooth_by_group.short = Q1_rate_smooth_short;
        all_Q2_rate_smooth_by_group.long = Q2_rate_smooth_long;
        all_Q2_rate_smooth_by_group.short = Q2_rate_smooth_short;
        all_Q3_rate_smooth_by_group.long = Q3_rate_smooth_long;
        all_Q3_rate_smooth_by_group.short = Q3_rate_smooth_short;
        all_bin_centers_by_group.long = plot_bin_centers_long;
        all_bin_centers_by_group.short = plot_bin_centers_short;
        all_avg_velocity_by_group.long = avg_velocity_long;
        all_avg_velocity_by_group.short = avg_velocity_short;
        all_occ_by_group.long = occ_long;
        all_occ_by_group.short = occ_short;
        all_rate_per_trial_by_group.long = rate_per_trial_long;
        all_rate_per_trial_by_group.short = rate_per_trial_short;
        all_rate_per_trial_smooth_by_group.long = rate_per_trial_smooth_long;
        all_rate_per_trial_smooth_by_group.short = rate_per_trial_smooth_short;
        all_spike_positions_by_group.long = spike_positions_long;
        all_spike_positions_by_group.short = spike_positions_short;
        all_global_trial_indices_by_group.long = global_trial_indices_long;
        all_global_trial_indices_by_group.short = global_trial_indices_short;
        all_rate_pooled_by_group.long = rate_pooled_long;
        all_rate_pooled_by_group.short = rate_pooled_short;
        
        % --- Reassemble spatial_tuning_stats structure from plain arrays ---
        if run_statistics
            spatial_tuning_stats = struct();
            for c = 1:n_clusters
                cluster_field = sprintf('cluster_%d', cluster_ids(c));
                spatial_tuning_stats.(cluster_field).long = struct( ...
                    'isSpatiallyTuned', stats_isSpatiallyTuned_long(c), ...
                    'reason', stats_reason_long{c}, ...
                    'n_spikes', stats_n_spikes_long(c), ...
                    'info', stats_info_long(c), ...
                    'pVal', stats_pVal_long(c), ...
                    'infoNull', stats_infoNull_long{c}, ...
                    'peakRate', stats_peakRate_long(c), ...
                    'peakPosition', stats_peakPosition_long(c), ...
                    'fieldBins', stats_fieldBins_long(c), ...
                    'numFields', stats_numFields_long(c), ...
                    'fieldStart', stats_fieldStart_long(c), ...
                    'fieldEnd', stats_fieldEnd_long(c));
                spatial_tuning_stats.(cluster_field).short = struct( ...
                    'isSpatiallyTuned', stats_isSpatiallyTuned_short(c), ...
                    'reason', stats_reason_short{c}, ...
                    'n_spikes', stats_n_spikes_short(c), ...
                    'info', stats_info_short(c), ...
                    'pVal', stats_pVal_short(c), ...
                    'infoNull', stats_infoNull_short{c}, ...
                    'peakRate', stats_peakRate_short(c), ...
                    'peakPosition', stats_peakPosition_short(c), ...
                    'fieldBins', stats_fieldBins_short(c), ...
                    'numFields', stats_numFields_short(c), ...
                    'fieldStart', stats_fieldStart_short(c), ...
                    'fieldEnd', stats_fieldEnd_short(c));
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
        
        % Reconstruct per-trial velocity and occupancy data for plotting
        trial_vel_by_group = struct();
        trial_occ_by_group = struct();
        
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
            
            n_trials = length(trials_struct);
            vel_per_trial = nan(n_trials, n_bins);
            occ_per_trial = nan(n_trials, n_bins);
            
            for k = 1:n_trials
                trial = trials_struct(k).trial;
                [~, occ, vel, ~] = compute_trial_firing_rate(trial, clusters(1), edges);
                vel_per_trial(k, :) = vel;
                occ_per_trial(k, :) = occ;
            end
            
            trial_vel_by_group.(group) = vel_per_trial;
            trial_occ_by_group.(group) = occ_per_trial;
        end
        
        % Average velocity plot with individual trials
        subplot(2, 1, 1, 'Parent', fig_profiles);
        hold on;
        for g = 1:2
            group = group_names{g};
            bin_centers = all_bin_centers_by_group.(group);
            
            % Plot individual trials with thin, semi-transparent lines
            vel_trials = trial_vel_by_group.(group);
            for k = 1:size(vel_trials, 1)
                if g == 1  % long trials - blue
                    plot(bin_centers, vel_trials(k, :), 'Color', [0.5 0.5 1], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                else  % short trials - red
                    plot(bin_centers, vel_trials(k, :), 'Color', [1 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                end
            end
            
            % Plot average with thick line
            avg_vel = nanmean(all_avg_velocity_by_group.(group), 1);
            if g == 1  % long trials - blue
                plot(bin_centers, avg_vel, 'Color', [0 0 0.8], 'LineWidth', 2, 'DisplayName', group_labels{g});
            else  % short trials - red
                plot(bin_centers, avg_vel, 'Color', [0.8 0 0], 'LineWidth', 2, 'DisplayName', group_labels{g});
            end
        end
        hold off;
        xlabel('Position (cm)');
        ylabel('Avg Velocity (cm/s)');
        legend('show', 'Location', 'northeastoutside');
        grid on;
        
        % Average occupancy plot with individual trials
        subplot(2, 1, 2, 'Parent', fig_profiles);
        hold on;
        for g = 1:2
            group = group_names{g};
            bin_centers = all_bin_centers_by_group.(group);
            
            % Plot individual trials with thin, semi-transparent lines
            occ_trials = trial_occ_by_group.(group);
            for k = 1:size(occ_trials, 1)
                if g == 1  % long trials - blue
                    plot(bin_centers, occ_trials(k, :), 'Color', [0.5 0.5 1], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                else  % short trials - red
                    plot(bin_centers, occ_trials(k, :), 'Color', [1 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                end
            end
            
            % Plot average with thick line
            avg_occ = nanmean(all_occ_by_group.(group), 1);
            if g == 1  % long trials - blue
                plot(bin_centers, avg_occ, 'Color', [0 0 0.8], 'LineWidth', 2, 'DisplayName', group_labels{g});
            else  % short trials - red
                plot(bin_centers, avg_occ, 'Color', [0.8 0 0], 'LineWidth', 2, 'DisplayName', group_labels{g});
            end
        end
        hold off;
        xlabel('Position (cm)');
        ylabel('Occupancy (s)');
        legend('show', 'Location', 'northeastoutside');
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
        
        % Adjust subplot positions BEFORE adding title to create significant space at top
        all_axes = findall(fig_heatmaps, 'Type', 'axes');
        for ax_idx = 1:length(all_axes)
            ax = all_axes(ax_idx);
            pos = get(ax, 'Position');
            % Shift all plots down to make room at top for title
            pos(4) = pos(4) * 0.82;  % Reduce height to 82% 
            pos(2) = max(0.08, pos(2) - 0.03);  % Lower the bottom position
            set(ax, 'Position', pos);
        end
        
        % Add figure title after adjusting subplot positions
        FigureTitle(fig_heatmaps, sprintf('Heatmaps - %s', probe_ids{pid}));
        
        ctl.figs.save_fig_to_join();
        
        % --- Plot spatially tuned clusters with Gaussian fits in 2x3 grid ---
        if exist('spatial_tuning_stats', 'var')
            fprintf('  Creating spatially tuned clusters plot with Gaussian fits...\n');
            
            % Identify clusters that are spatially tuned in AT LEAST ONE condition
            spatially_tuned_any = [];
            
            for c = 1:n_clusters
                cluster_field = sprintf('cluster_%d', cluster_ids(c));
                if isfield(spatial_tuning_stats, cluster_field)
                    stats = spatial_tuning_stats.(cluster_field);
                    % Must be tuned in at least one group
                    if (isfield(stats, 'long') && stats.long.isSpatiallyTuned) || ...
                       (isfield(stats, 'short') && stats.short.isSpatiallyTuned)
                        spatially_tuned_any = [spatially_tuned_any; c];
                    end
                end
            end
            
            n_tuned = length(spatially_tuned_any);
            fprintf('    Found %d clusters spatially tuned in at least one condition\n', n_tuned);
            
            if n_tuned > 0
                % Sort clusters by peak position in long trials (or short if not tuned in long)
                peak_pos = zeros(n_tuned, 1);
                for i = 1:n_tuned
                    c = spatially_tuned_any(i);
                    cluster_field = sprintf('cluster_%d', cluster_ids(c));
                    stats = spatial_tuning_stats.(cluster_field);
                    if isfield(stats, 'long') && stats.long.isSpatiallyTuned
                        peak_pos(i) = stats.long.peakPosition;
                    elseif isfield(stats, 'short') && stats.short.isSpatiallyTuned
                        peak_pos(i) = stats.short.peakPosition;
                    end
                end
                [~, sort_idx] = sort(peak_pos);
                spatially_tuned_sorted = spatially_tuned_any(sort_idx);
                
                % Generate pastel colors for all clusters (ensure unique colors)
                pastel_colors = generate_pastel_colors(n_tuned);
                
                % --- Create 2x3 grid plot ---
                fprintf('  Creating 2x3 grid visualization of spatially tuned clusters...\n');
                fig_grid = figure('Position', [100, 100, 1200, 900], 'PaperPositionMode', 'auto');
                
                % Define Gaussian fit function
                gauss_fit = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2 * p(3)^2));
                options = optimset('Display', 'off');
                
                % Row 1, Column 1: Long trials (unnormalized)
                subplot('Position', [0.08, 0.72, 0.38, 0.20], 'Parent', fig_grid);  % Row 1, left half
                hold on;
                for i = 1:n_tuned
                    c = spatially_tuned_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    bin_centers_long = all_bin_centers_by_group.long;
                    median_rate_long = all_Q2_rate_smooth_by_group.long(c, :);
                    
                    % Plot original median with reduced alpha
                    h = plot(bin_centers_long, median_rate_long, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('C%d', cluster_id));
                    h.Color(4) = 0.6;  % Set alpha to 0.6
                    
                    % Fit and plot Gaussian if tuned in long
                    if isfield(spatial_tuning_stats.(cluster_field), 'long') && spatial_tuning_stats.(cluster_field).long.isSpatiallyTuned
                        peak_pos = spatial_tuning_stats.(cluster_field).long.peakPosition;
                        peak_rate = spatial_tuning_stats.(cluster_field).long.peakRate;
                        p0 = [peak_rate, peak_pos, gauss_sigma_cm];
                        try
                            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_long, median_rate_long, ...
                                               [0, 0, 1], [inf, 120, 50], options);
                            x_fit = linspace(0, 120, 300);
                            y_fit = gauss_fit(p_fit, x_fit);
                            % Fill area under Gaussian with pale color and low alpha
                            fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, ...
                                 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        catch
                            % Fit failed, median already plotted
                        end
                    end
                end
                xlabel('Position (cm)');
                ylabel('Firing Rate (Hz)');
                title('Long Trials (0-120 cm)');
                grid on;
                xlim([0, 120]);
                hold off;
                
                % Row 1, Column 2: Long trials (normalized)
                subplot('Position', [0.54, 0.72, 0.38, 0.20], 'Parent', fig_grid);  % Row 1, right half
                hold on;
                for i = 1:n_tuned
                    c = spatially_tuned_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    bin_centers_long = all_bin_centers_by_group.long;
                    median_rate_long = all_Q2_rate_smooth_by_group.long(c, :);
                    
                    % Normalize to [0, 1]
                    min_rate = min(median_rate_long);
                    max_rate = max(median_rate_long);
                    if max_rate > min_rate
                        median_rate_norm = (median_rate_long - min_rate) / (max_rate - min_rate);
                    else
                        median_rate_norm = zeros(size(median_rate_long));
                    end
                    
                    % Plot normalized median with reduced alpha
                    h = plot(bin_centers_long, median_rate_norm, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('C%d', cluster_id));
                    h.Color(4) = 0.6;  % Set alpha to 0.6
                    
                    % Fit and plot normalized Gaussian if tuned in long
                    if isfield(spatial_tuning_stats.(cluster_field), 'long') && spatial_tuning_stats.(cluster_field).long.isSpatiallyTuned
                        peak_pos = spatial_tuning_stats.(cluster_field).long.peakPosition;
                        peak_rate_norm = 1.0;  % Normalized peak is always 1
                        p0 = [peak_rate_norm, peak_pos, gauss_sigma_cm];
                        try
                            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_long, median_rate_norm, ...
                                               [0, 0, 1], [2, 120, 50], options);
                            x_fit = linspace(0, 120, 300);
                            y_fit = gauss_fit(p_fit, x_fit);
                            % Fill area under Gaussian with pale color and low alpha
                            fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, ...
                                 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        catch
                            % Fit failed, median already plotted
                        end
                    end
                end
                xlabel('Position (cm)');
                ylabel('Normalized Firing Rate');
                title('Long Trials (normalized)');
                grid on;
                xlim([0, 120]);
                ylim([0, 1]);
                hold off;
                
                % Row 2, Column 2: Short trials (unnormalized, 60-120 cm) - 60% width in position 2
                subplot('Position', [0.23, 0.42, 0.23, 0.20], 'Parent', fig_grid);  % Row 2, position 2 (40% gap + 60% plot)
                hold on;
                for i = 1:n_tuned
                    c = spatially_tuned_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    bin_centers_short = all_bin_centers_by_group.short;
                    median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
                    
                    % Plot original median with reduced alpha
                    h = plot(bin_centers_short, median_rate_short, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('C%d', cluster_id));
                    h.Color(4) = 0.6;  % Set alpha to 0.6
                    
                    % Fit and plot Gaussian if tuned in short
                    if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
                        peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
                        peak_rate = spatial_tuning_stats.(cluster_field).short.peakRate;
                        p0 = [peak_rate, peak_pos, gauss_sigma_cm];
                        try
                            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_short, ...
                                               [0, 60, 1], [inf, 120, 50], options);
                            x_fit = linspace(60, 120, 300);
                            y_fit = gauss_fit(p_fit, x_fit);
                            % Fill area under Gaussian with pale color and low alpha
                            fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, ...
                                 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        catch
                            % Fit failed, median already plotted
                        end
                    end
                end
                xlabel('Position (cm)');
                ylabel('Firing Rate (Hz)');
                title('Short Trials (60-120 cm)');
                grid on;
                xlim([60, 120]);
                hold off;
                
                % Row 2, Column 4: Short trials (normalized, 60-120 cm) - 60% width in position 4
                subplot('Position', [0.69, 0.42, 0.23, 0.20], 'Parent', fig_grid);  % Row 2, position 4 (40% gap + 60% plot)
                hold on;
                for i = 1:n_tuned
                    c = spatially_tuned_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    bin_centers_short = all_bin_centers_by_group.short;
                    median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
                    
                    % Normalize to [0, 1]
                    min_rate = min(median_rate_short);
                    max_rate = max(median_rate_short);
                    if max_rate > min_rate
                        median_rate_norm = (median_rate_short - min_rate) / (max_rate - min_rate);
                    else
                        median_rate_norm = zeros(size(median_rate_short));
                    end
                    
                    % Plot normalized median with reduced alpha
                    h = plot(bin_centers_short, median_rate_norm, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('C%d', cluster_id));
                    h.Color(4) = 0.6;  % Set alpha to 0.6
                    
                    % Fit and plot normalized Gaussian if tuned in short
                    if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
                        peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
                        peak_rate_norm = 1.0;
                        p0 = [peak_rate_norm, peak_pos, gauss_sigma_cm];
                        try
                            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_norm, ...
                                               [0, 60, 1], [2, 120, 50], options);
                            x_fit = linspace(60, 120, 300);
                            y_fit = gauss_fit(p_fit, x_fit);
                            % Fill area under Gaussian with pale color and low alpha
                            fill([x_fit, fliplr(x_fit)], [y_fit, zeros(size(y_fit))], color, ...
                                 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        catch
                            % Fit failed, median already plotted
                        end
                    end
                end
                xlabel('Position (cm)');
                ylabel('Normalized Firing Rate');
                title('Short Trials (normalized, 60-120 cm)');
                grid on;
                xlim([60, 120]);
                ylim([0, 1]);
                hold off;
                
                % Row 3, Column 1: Short trials with percentage x-axis (data unchanged, spans full width)
                subplot('Position', [0.08, 0.12, 0.38, 0.20], 'Parent', fig_grid);  % Row 3, left half
                hold on;
                for i = 1:n_tuned
                    c = spatially_tuned_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    bin_centers_short = all_bin_centers_by_group.short;  % Original 60-120 cm
                    median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
                    
                    % Convert position to percentage (60cm=0%, 120cm=100%)
                    bin_centers_percent = (bin_centers_short - 60) / 60 * 100;
                    
                    % Plot original data with percentage x-axis and reduced alpha
                    h = plot(bin_centers_percent, median_rate_short, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('C%d', cluster_id));
                    h.Color(4) = 0.6;  % Set alpha to 0.6
                    
                    % Plot Gaussian with SAME fit as row 2, just converted to percentage x-axis
                    if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
                        peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
                        peak_rate = spatial_tuning_stats.(cluster_field).short.peakRate;
                        p0 = [peak_rate, peak_pos, gauss_sigma_cm];
                        try
                            % Fit using ORIGINAL position coordinates (same as row 2)
                            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_short, ...
                                               [0, 60, 1], [inf, 120, 50], options);
                            % Generate x values in original coordinates then convert to percentage
                            x_fit_cm = linspace(60, 120, 300);
                            x_fit_percent = (x_fit_cm - 60) / 60 * 100;
                            y_fit = gauss_fit(p_fit, x_fit_cm);
                            % Fill area under Gaussian with pale color and low alpha
                            fill([x_fit_percent, fliplr(x_fit_percent)], [y_fit, zeros(size(y_fit))], color, ...
                                 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        catch
                            % Fit failed, median already plotted
                        end
                    end
                end
                xlabel('Track Position (%)');
                ylabel('Firing Rate (Hz)');
                title('Short Trials (0-100% of track)');
                grid on;
                xlim([0, 100]);
                hold off;
                
                % Row 3, Column 2: Short trials normalized with percentage x-axis (data unchanged, spans full width)
                subplot('Position', [0.54, 0.12, 0.38, 0.20], 'Parent', fig_grid);  % Row 3, right half
                hold on;
                for i = 1:n_tuned
                    c = spatially_tuned_sorted(i);
                    cluster_id = cluster_ids(c);
                    cluster_field = sprintf('cluster_%d', cluster_id);
                    color = pastel_colors(i, :);
                    
                    bin_centers_short = all_bin_centers_by_group.short;  % Original 60-120 cm
                    median_rate_short = all_Q2_rate_smooth_by_group.short(c, :);
                    
                    % Convert position to percentage (60cm=0%, 120cm=100%)
                    bin_centers_percent = (bin_centers_short - 60) / 60 * 100;
                    
                    % Normalize to [0, 1]
                    min_rate = min(median_rate_short);
                    max_rate = max(median_rate_short);
                    if max_rate > min_rate
                        median_rate_norm = (median_rate_short - min_rate) / (max_rate - min_rate);
                    else
                        median_rate_norm = zeros(size(median_rate_short));
                    end
                    
                    % Plot normalized data with percentage x-axis and reduced alpha
                    h = plot(bin_centers_percent, median_rate_norm, 'o-', 'Color', color, 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('C%d', cluster_id));
                    h.Color(4) = 0.6;  % Set alpha to 0.6
                    
                    % Plot Gaussian with SAME fit as row 2, just normalized and converted to percentage
                    if isfield(spatial_tuning_stats.(cluster_field), 'short') && spatial_tuning_stats.(cluster_field).short.isSpatiallyTuned
                        peak_pos = spatial_tuning_stats.(cluster_field).short.peakPosition;
                        peak_rate_norm = 1.0;
                        p0 = [peak_rate_norm, peak_pos, gauss_sigma_cm];
                        try
                            % Fit using ORIGINAL position coordinates (same as row 2)
                            p_fit = lsqcurvefit(gauss_fit, p0, bin_centers_short, median_rate_norm, ...
                                               [0, 60, 1], [2, 120, 50], options);
                            % Generate x values in original coordinates then convert to percentage
                            x_fit_cm = linspace(60, 120, 300);
                            x_fit_percent = (x_fit_cm - 60) / 60 * 100;
                            y_fit = gauss_fit(p_fit, x_fit_cm);
                            % Fill area under Gaussian with pale color and low alpha
                            fill([x_fit_percent, fliplr(x_fit_percent)], [y_fit, zeros(size(y_fit))], color, ...
                                 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        catch
                            % Fit failed, median already plotted
                        end
                    end
                end
                xlabel('Track Position (%)');
                ylabel('Normalized Firing Rate');
                title('Short Trials (normalized, 0-100% of track)');
                grid on;
                xlim([0, 100]);
                ylim([0, 1]);
                hold off;
                
                % Adjust subplot positions BEFORE adding title
                % All subplots were positioned manually, so adjust them to create top space
                all_axes = findall(fig_grid, 'Type', 'axes');
                for ax_idx = 1:length(all_axes)
                    ax = all_axes(ax_idx);
                    pos = get(ax, 'Position');
                    % Move all plots down by reducing top position
                    pos(4) = pos(4) * 0.88;  % Reduce height 
                    pos(2) = pos(2) * 0.85;  % Scale down vertical position to compress towards bottom
                    set(ax, 'Position', pos);
                end
                
                % Add figure title after adjusting positions
                FigureTitle(fig_grid, sprintf('Spatial Tuning Grid - %s', probe_ids{pid}));
                
                ctl.figs.save_fig_to_join();
            else
                fprintf('    No spatially tuned clusters found. Skipping plot.\n');
            end
            
            % --- Generate 3D visualization ---
            fprintf('  Generating 3D spatial tuning visualization...\n');
            spatial_firing_rate_profile_3d;
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
    
    % --- Clean up memory after processing each probe ---
    % Clear large probe-specific variables to free RAM
    clear data clusters sessions trial_groups
    clear all_rate_smooth_by_group all_bin_centers_by_group all_avg_velocity_by_group
    clear all_Q1_rate_smooth_by_group all_Q2_rate_smooth_by_group all_Q3_rate_smooth_by_group
    clear all_occ_by_group all_rate_per_trial_by_group all_rate_per_trial_smooth_by_group
    clear all_spike_positions_by_group all_global_trial_indices_by_group all_rate_pooled_by_group
    clear spatial_tuning_stats
    
    % Clear intermediate parfor arrays if they exist
    clear rate_smooth_long rate_smooth_short
    clear Q1_rate_smooth_long Q1_rate_smooth_short
    clear Q2_rate_smooth_long Q2_rate_smooth_short
    clear Q3_rate_smooth_long Q3_rate_smooth_short
    clear avg_velocity_long avg_velocity_short
    clear occ_long occ_short
    clear rate_pooled_long rate_pooled_short
    clear rate_per_trial_long rate_per_trial_short
    clear rate_per_trial_smooth_long rate_per_trial_smooth_short
    clear spike_positions_long spike_positions_short
    clear global_trial_indices_long global_trial_indices_short
    
    % Clear statistical arrays if they exist
    clear stats_isSpatiallyTuned_long stats_isSpatiallyTuned_short
    clear stats_reason_long stats_reason_short
    clear stats_n_spikes_long stats_n_spikes_short
    clear stats_info_long stats_info_short
    clear stats_pVal_long stats_pVal_short
    clear stats_infoNull_long stats_infoNull_short
    clear stats_peakRate_long stats_peakRate_short
    clear stats_peakPosition_long stats_peakPosition_short
    clear stats_fieldBins_long stats_fieldBins_short
    clear stats_numFields_long stats_numFields_short
    clear stats_fieldStart_long stats_fieldStart_short
    clear stats_fieldEnd_long stats_fieldEnd_short
    
    fprintf('  Memory cleaned up after probe %d/%d\n', pid, length(probe_ids));
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