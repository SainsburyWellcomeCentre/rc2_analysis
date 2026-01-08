classdef SpatialTuningAnalyzer < handle
    % SPATIALTUN INGANALYZER Analyzes spatial tuning of neural clusters
    %
    %   Encapsulates spatial firing rate analysis for neural clusters across
    %   different trial groups (long/short). Handles computation, statistical
    %   testing, and result management.
    %
    %   Properties:
    %       bin_size_cm         - Spatial bin size (cm)
    %       gauss_sigma_cm      - Gaussian smoothing kernel sigma (cm)
    %       position_threshold  - Threshold for long/short classification (cm)
    %       stats_params        - Statistical analysis parameters
    %       use_parallel        - Enable parallel processing
    %       max_workers         - Maximum parallel workers
    %
    %   Methods:
    %       analyze_all_clusters    - Compute spatial firing for all clusters
    %       classify_spatial_tuning - Perform statistical testing
    %       save_cache              - Save results to cache file
    %       load_cache              - Load results from cache file
    %       get_cluster_data        - Retrieve data for specific cluster
    %
    %   Example:
    %       analyzer = SpatialTuningAnalyzer(clusters, sessions);
    %       analyzer.bin_size_cm = 2;
    %       analyzer.gauss_sigma_cm = 8;
    %       analyzer.analyze_all_clusters();
    %       data = analyzer.get_cluster_data(cluster_id);
    
    properties
        % Configuration parameters
        bin_size_cm         = 2     % Spatial bin size (cm)
        gauss_sigma_cm      = 8     % Gaussian smoothing sigma (cm)
        position_threshold  = 90    % Long/short classification threshold (cm)
        
        % Statistical parameters
        stats_params = struct(...
            'minSpikes', 25, ...        % Minimum spikes for classification
            'minPeakRate', 1.0, ...     % Minimum peak rate (Hz)
            'fieldFrac', 0.7, ...       % Field threshold fraction
            'minFieldBins', 5, ...      % Minimum field size (bins)
            'maxNumFields', 1, ...      % Maximum number of fields
            'nShuf', 100, ...           % Number of shuffles
            'pThresh', 0.05, ...        % Significance threshold
            'run_statistics', true)     % Enable/disable statistics
        
        % Parallel processing
        use_parallel = true
        max_workers = 4
        
        % Input data
        clusters            % Array of cluster objects
        sessions            % Cell array of session objects
        
        % Processed data
        trial_groups        % Struct with 'long' and 'short' trial arrays
        group_names         % {'long', 'short'}
        group_labels        % Descriptive labels for plotting
        
        % Bin configuration (per group)
        bin_config          % Struct: long/short -> edges, centers, n_bins, plot_centers
        
        % Results (per cluster, per group)
        cluster_ids         % Array of cluster IDs
        firing_rates        % Struct: long/short -> matrices for all clusters
        tuning_stats        % Struct: cluster_id -> long/short -> stats
        ttg_analysis        % Struct: per-cluster TTG data and statistics
        distribution_comparisons % Cell array: per-cluster Kruskal-Wallis results
    end
    
    methods
        function obj = SpatialTuningAnalyzer(clusters, sessions)
            % Constructor
            %
            %   analyzer = SpatialTuningAnalyzer(clusters, sessions)
            %
            %   Inputs:
            %       clusters - Array of cluster objects
            %       sessions - Cell array of session objects
            
            if nargin > 0
                obj.clusters = clusters;
                obj.sessions = sessions;
                obj.cluster_ids = [clusters.id];
                obj.group_names = {'long', 'short'};
                obj.group_labels = {'Long (0-120 cm)', 'Short (60-120 cm)'};
            end
        end
        
        function setup_bins(obj)
            % Setup bin configuration for each trial group
            
            % Long trials: 0-120 cm
            edges_long = 0:obj.bin_size_cm:120;
            bin_centers_long = edges_long(1:end-1) + obj.bin_size_cm/2;
            n_bins_long = length(bin_centers_long);
            
            % Short trials: 0-60 cm (plotted at 60-120 cm)
            edges_short = 0:obj.bin_size_cm:60;
            bin_centers_short = edges_short(1:end-1) + obj.bin_size_cm/2;
            n_bins_short = length(bin_centers_short);
            
            obj.bin_config = struct(...
                'long', struct(...
                    'edges', edges_long, ...
                    'centers', bin_centers_long, ...
                    'n_bins', n_bins_long, ...
                    'plot_centers', bin_centers_long), ...
                'short', struct(...
                    'edges', edges_short, ...
                    'centers', bin_centers_short, ...
                    'n_bins', n_bins_short, ...
                    'plot_centers', bin_centers_short + 60));
        end
        
        function collect_trials(obj)
            % Collect and group trials by track length
            [obj.trial_groups, obj.group_names, obj.group_labels] = ...
                collect_and_group_trials(obj.sessions, obj.position_threshold);
        end
        
        function analyze_all_clusters(obj)
            % Compute spatial firing rate profiles for all clusters
            %
            %   analyzer.analyze_all_clusters()
            %
            %   Processes all clusters and computes spatial firing rates,
            %   quartiles, and optionally statistical significance.
            
            % Setup if not already done
            if isempty(obj.trial_groups)
                obj.collect_trials();
            end
            if isempty(obj.bin_config)
                obj.setup_bins();
            end
            
            n_clusters = length(obj.clusters);
            
            fprintf('  Computing spatial firing rate profiles for %d clusters (%s)...\n', ...
                n_clusters, obj.parallel_status());
            
            % Initialize result storage
            obj.initialize_storage();
            
            % Setup parallel pool if needed
            if obj.use_parallel
                obj.setup_parallel_pool();
            end
            
            % Process all clusters
            cluster_results = cell(n_clusters, 1);
            
            if obj.use_parallel
                parfor c = 1:n_clusters
                    cluster_results{c} = obj.compute_cluster_data(c);
                end
            else
                for c = 1:n_clusters
                    cluster_results{c} = obj.compute_cluster_data(c);
                    if mod(c, 10) == 0
                        fprintf('    Processed %d/%d clusters\n', c, n_clusters);
                    end
                end
            end
            
            % Store results
            obj.store_results(cluster_results);
            
            fprintf('  Analysis complete for %d clusters\n', n_clusters);
        end
        
        function result = compute_cluster_data(obj, cluster_idx)
            % Compute spatial data for a single cluster (parfor-compatible)
            %
            %   This method is designed to be called within parfor loops
            %   and does not modify object properties directly.
            
            cluster = obj.clusters(cluster_idx);
            result = struct();
            result.cluster_id = cluster.id;
            result.cluster_idx = cluster_idx;
            
            % Process each group (long/short)
            for g = 1:2
                group = obj.group_names{g};
                trials_struct = obj.trial_groups.(group);
                cfg = obj.bin_config.(group);
                
                % Compute trial-by-trial firing rates
                [rate_data, spike_data] = obj.compute_group_firing_rates(...
                    cluster, trials_struct, cfg);
                
                % Store rate data
                result.(group).rate_smooth = rate_data.mean_rate_smooth;
                result.(group).Q1_smooth = rate_data.Q1_smooth;
                result.(group).Q2_smooth = rate_data.Q2_smooth;
                result.(group).Q3_smooth = rate_data.Q3_smooth;
                result.(group).avg_velocity = rate_data.avg_velocity;
                result.(group).occupancy = rate_data.occupancy;
                result.(group).rate_per_trial = spike_data.rate_per_trial;
                result.(group).rate_per_trial_smooth = spike_data.rate_per_trial_smooth;
                result.(group).spike_positions = spike_data.spike_positions;
                result.(group).global_trial_indices = spike_data.global_trial_indices;
                
                % Statistical analysis
                if obj.stats_params.run_statistics
                    stats = obj.compute_spatial_statistics(...
                        rate_data, spike_data, cfg, group);
                    result.(group).stats = stats;
                else
                    result.(group).stats = [];
                end
            end
        end
        
        function [rate_data, spike_data] = compute_group_firing_rates(obj, cluster, trials_struct, cfg)
            % Compute firing rates for all trials in a group
            
            n_trials = length(trials_struct);
            n_bins = cfg.n_bins;
            edges = cfg.edges;
            
            % Pre-allocate
            rate_per_trial = nan(n_trials, n_bins);
            occ_per_trial = nan(n_trials, n_bins);
            vel_per_trial = nan(n_trials, n_bins);
            count_per_trial = nan(n_trials, n_bins);
            spike_positions_per_trial = cell(n_trials, 1);
            global_trial_indices = zeros(n_trials, 1);
            
            % Compute for each trial
            for k = 1:n_trials
                trial = trials_struct(k).trial;
                global_trial_indices(k) = trials_struct(k).global_trial_idx;
                
                [rate, occ, vel, counts] = compute_trial_firing_rate(trial, cluster, edges);
                rate_per_trial(k, :) = rate;
                occ_per_trial(k, :) = occ;
                vel_per_trial(k, :) = vel;
                count_per_trial(k, :) = counts;
                
                % Get spike positions for raster plots
                spike_positions_per_trial{k} = obj.extract_spike_positions(trial, cluster);
            end
            
            % Smooth counts and occupancy separately, then compute rates
            count_per_trial_smooth = nan(n_trials, n_bins);
            occ_per_trial_smooth = nan(n_trials, n_bins);
            rate_per_trial_smooth = nan(n_trials, n_bins);
            
            for k = 1:n_trials
                count_per_trial_smooth(k, :) = smooth_spatial_rate(...
                    count_per_trial(k, :), obj.bin_size_cm, obj.gauss_sigma_cm);
                occ_per_trial_smooth(k, :) = smooth_spatial_rate(...
                    occ_per_trial(k, :), obj.bin_size_cm, obj.gauss_sigma_cm);
                
                rate_per_trial_smooth(k, :) = count_per_trial_smooth(k, :) ./ occ_per_trial_smooth(k, :);
                rate_per_trial_smooth(k, isnan(rate_per_trial_smooth(k, :)) | isinf(rate_per_trial_smooth(k, :))) = 0;
            end
            
            % Compute quartiles from smoothed trial data
            Q1_smooth = prctile(rate_per_trial_smooth, 25, 1);
            Q2_smooth = prctile(rate_per_trial_smooth, 50, 1);
            Q3_smooth = prctile(rate_per_trial_smooth, 75, 1);
            
            % Compute mean across trials
            mean_count_smooth = nanmean(count_per_trial_smooth, 1);
            mean_occ_smooth = nanmean(occ_per_trial_smooth, 1);
            mean_rate_smooth = mean_count_smooth ./ mean_occ_smooth;
            mean_rate_smooth(isnan(mean_rate_smooth) | isinf(mean_rate_smooth)) = 0;
            
            % Compute mean occupancy and velocity
            mean_occupancy = nanmean(occ_per_trial, 1);
            mean_velocity = nanmean(vel_per_trial, 1);
            
            % Package results
            rate_data = struct(...
                'mean_rate_smooth', mean_rate_smooth, ...
                'Q1_smooth', Q1_smooth, ...
                'Q2_smooth', Q2_smooth, ...
                'Q3_smooth', Q3_smooth, ...
                'avg_velocity', mean_velocity, ...
                'occupancy', mean_occupancy);
            
            spike_data = struct(...
                'rate_per_trial', rate_per_trial, ...
                'rate_per_trial_smooth', rate_per_trial_smooth, ...
                'count_per_trial', count_per_trial, ...
                'occ_per_trial_smooth', occ_per_trial_smooth, ...
                'spike_positions', {spike_positions_per_trial}, ...
                'global_trial_indices', global_trial_indices);
        end
        
        function spike_pos = extract_spike_positions(obj, trial, cluster)
            % Extract spike positions for raster plots
            
            motion_mask = trial.motion_mask();
            pos = trial.position(motion_mask);
            tvec = trial.probe_t(motion_mask);
            
            if ~isempty(pos) && ~isempty(tvec)
                st = cluster.spike_times;
                spike_mask = st >= tvec(1) & st <= tvec(end);
                st_in_trial = st(spike_mask);
                
                if ~isempty(st_in_trial)
                    spike_pos = interp1(tvec, pos, st_in_trial, 'linear', 'extrap');
                    spike_pos = spike_pos(:)';
                else
                    spike_pos = [];
                end
            else
                spike_pos = [];
            end
        end
        
        function stats = compute_spatial_statistics(obj, rate_data, spike_data, cfg, group)
            % Perform statistical testing for spatial tuning
            
            % Count total spikes
            n_spikes = nansum(spike_data.count_per_trial(:));
            
            % Check minimum spike count
            if n_spikes < obj.stats_params.minSpikes
                stats = obj.create_insufficient_spikes_stats(n_spikes);
                return;
            end
            
            % Use median rate and total occupancy for testing
            median_rate_smooth = nanmedian(spike_data.rate_per_trial_smooth, 1);
            total_occ_smooth = nansum(spike_data.occ_per_trial_smooth, 1);
            
            % Compute Skaggs information
            infoObs = skaggs_info(median_rate_smooth, total_occ_smooth);
            
            % Bin-shuffle test
            infoNull = compute_shuffled_spatial_info(...
                spike_data.rate_per_trial_smooth, total_occ_smooth, obj.stats_params.nShuf);
            
            % Compute p-value
            pVal = (sum(infoNull >= infoObs) + 1) / (obj.stats_params.nShuf + 1);
            
            % Peak detection and field analysis
            [peakRate, peakIdx] = max(median_rate_smooth);
            peakPosition = cfg.centers(peakIdx);
            
            % Field detection
            [longestField, numFields, fieldStart, fieldEnd] = detect_spatial_fields(...
                median_rate_smooth, cfg.centers, obj.bin_size_cm, ...
                obj.stats_params.fieldFrac, obj.stats_params.minFieldBins);
            
            % Classification
            isSpatiallyTuned = (pVal < obj.stats_params.pThresh) && ...
                              (longestField >= obj.stats_params.minFieldBins) && ...
                              (numFields <= obj.stats_params.maxNumFields) && ...
                              (peakRate >= obj.stats_params.minPeakRate);
            
            reasonStr = 'passes';
            if ~isSpatiallyTuned
                if pVal >= obj.stats_params.pThresh
                    reasonStr = 'info_not_significant';
                elseif peakRate < obj.stats_params.minPeakRate
                    reasonStr = 'low_peak';
                elseif longestField < obj.stats_params.minFieldBins
                    reasonStr = 'small_field';
                elseif numFields > obj.stats_params.maxNumFields
                    reasonStr = 'multiple_fields';
                end
            end
            
            % Package results
            stats = struct(...
                'isSpatiallyTuned', isSpatiallyTuned, ...
                'reason', reasonStr, ...
                'n_spikes', n_spikes, ...
                'info', infoObs, ...
                'pVal', pVal, ...
                'infoNull', infoNull, ...
                'peakRate', peakRate, ...
                'peakPosition', peakPosition, ...
                'fieldBins', longestField, ...
                'numFields', numFields, ...
                'fieldStart', fieldStart, ...
                'fieldEnd', fieldEnd);
        end
        
        function stats = create_insufficient_spikes_stats(obj, n_spikes)
            % Create stats struct for clusters with insufficient spikes
            stats = struct(...
                'isSpatiallyTuned', false, ...
                'reason', 'insufficient_spikes', ...
                'n_spikes', n_spikes, ...
                'info', NaN, ...
                'pVal', NaN, ...
                'infoNull', [], ...
                'peakRate', NaN, ...
                'peakPosition', NaN, ...
                'fieldBins', 0, ...
                'numFields', 0, ...
                'fieldStart', NaN, ...
                'fieldEnd', NaN);
        end
        
        function initialize_storage(obj)
            % Initialize storage structures for results
            obj.firing_rates = struct();
            obj.tuning_stats = struct();
            
            for g = 1:2
                group = obj.group_names{g};
                n_bins = obj.bin_config.(group).n_bins;
                n_clusters = length(obj.clusters);
                
                obj.firing_rates.(group) = struct(...
                    'rate_smooth', nan(n_clusters, n_bins), ...
                    'Q1_smooth', nan(n_clusters, n_bins), ...
                    'Q2_smooth', nan(n_clusters, n_bins), ...
                    'Q3_smooth', nan(n_clusters, n_bins), ...
                    'avg_velocity', nan(n_clusters, n_bins), ...
                    'occupancy', nan(n_clusters, n_bins), ...
                    'rate_per_trial', {cell(n_clusters, 1)}, ...
                    'rate_per_trial_smooth', {cell(n_clusters, 1)}, ...
                    'spike_positions', {cell(n_clusters, 1)}, ...
                    'global_trial_indices', {cell(n_clusters, 1)});
            end
        end
        
        function store_results(obj, cluster_results)
            % Store results from parallel computation into object properties
            
            n_clusters = length(cluster_results);
            
            for c = 1:n_clusters
                result = cluster_results{c};
                cluster_id = result.cluster_id;
                cluster_field = sprintf('cluster_%d', cluster_id);
                
                for g = 1:2
                    group = obj.group_names{g};
                    
                    % Store firing rate matrices
                    obj.firing_rates.(group).rate_smooth(c, :) = result.(group).rate_smooth;
                    obj.firing_rates.(group).Q1_smooth(c, :) = result.(group).Q1_smooth;
                    obj.firing_rates.(group).Q2_smooth(c, :) = result.(group).Q2_smooth;
                    obj.firing_rates.(group).Q3_smooth(c, :) = result.(group).Q3_smooth;
                    obj.firing_rates.(group).avg_velocity(c, :) = result.(group).avg_velocity;
                    obj.firing_rates.(group).occupancy(c, :) = result.(group).occupancy;
                    
                    % Store per-trial data
                    obj.firing_rates.(group).rate_per_trial{c} = result.(group).rate_per_trial;
                    obj.firing_rates.(group).rate_per_trial_smooth{c} = result.(group).rate_per_trial_smooth;
                    obj.firing_rates.(group).spike_positions{c} = result.(group).spike_positions;
                    obj.firing_rates.(group).global_trial_indices{c} = result.(group).global_trial_indices;
                    
                    % Store statistics
                    if ~isempty(result.(group).stats)
                        obj.tuning_stats.(cluster_field).(group) = result.(group).stats;
                    end
                end
            end
        end
        
        function rate_data = get_cluster_data(obj, cluster_id, group)
            % Retrieve firing rate data for a specific cluster
            %
            %   rate_data = analyzer.get_cluster_data(cluster_id)
            %   rate_data = analyzer.get_cluster_data(cluster_id, 'long')
            %
            %   Returns a struct with firing rate data for the specified
            %   cluster. If group is specified, returns only that group's data.
            
            c = find(obj.cluster_ids == cluster_id, 1);
            if isempty(c)
                error('Cluster %d not found', cluster_id);
            end
            
            if nargin < 3
                % Return both groups
                rate_data = struct();
                for g = 1:2
                    grp = obj.group_names{g};
                    rate_data.(grp) = obj.extract_cluster_group_data(c, grp);
                end
            else
                % Return specific group
                rate_data = obj.extract_cluster_group_data(c, group);
            end
        end
        
        function data = extract_cluster_group_data(obj, cluster_idx, group)
            % Extract data for one cluster and one group
            data = struct(...
                'Q1_smooth', obj.firing_rates.(group).Q1_smooth(cluster_idx, :), ...
                'Q2_smooth', obj.firing_rates.(group).Q2_smooth(cluster_idx, :), ...
                'Q3_smooth', obj.firing_rates.(group).Q3_smooth(cluster_idx, :), ...
                'rate_per_trial', obj.firing_rates.(group).rate_per_trial{cluster_idx}, ...
                'rate_per_trial_smooth', obj.firing_rates.(group).rate_per_trial_smooth{cluster_idx}, ...
                'spike_positions', obj.firing_rates.(group).spike_positions{cluster_idx}, ...
                'global_trial_indices', obj.firing_rates.(group).global_trial_indices{cluster_idx});
        end
        
        function stats = get_cluster_stats(obj, cluster_id, group)
            % Retrieve statistics for a specific cluster
            %
            %   stats = analyzer.get_cluster_stats(cluster_id)
            %   stats = analyzer.get_cluster_stats(cluster_id, 'long')
            
            cluster_field = sprintf('cluster_%d', cluster_id);
            if ~isfield(obj.tuning_stats, cluster_field)
                stats = [];
                return;
            end
            
            if nargin < 3
                stats = obj.tuning_stats.(cluster_field);
            else
                stats = obj.tuning_stats.(cluster_field).(group);
            end
        end
        
        function [rate_map_2d, vel_bin_edges, trial_count_map] = compute_2d_position_velocity_tuning(obj, cluster, group, n_vel_bins, sigma_velocity)
            % Compute 2D position × velocity firing rate map with equal-occupancy velocity bins
            %
            %   [rate_map_2d, vel_bin_edges, trial_count_map] = analyzer.compute_2d_position_velocity_tuning(cluster, group, n_vel_bins, sigma_velocity)
            %
            %   Inputs:
            %       cluster         - Cluster object
            %       group           - 'long' or 'short'
            %       n_vel_bins      - Number of velocity bins (will have equal occupancy)
            %       sigma_velocity  - Smoothing sigma for velocity dimension (cm/s)
            %
            %   Outputs:
            %       rate_map_2d     - 2D firing rate map (position × velocity), median across trials
            %       vel_bin_edges   - Velocity bin edges (variable width for equal occupancy)
            %       trial_count_map - Number of trials that sampled each bin (for alpha transparency)
            
            if nargin < 5
                sigma_velocity = 5;  % Default smoothing sigma for velocity (cm/s)
            end
            
            % Get position bin configuration for this group
            pos_edges = obj.bin_config.(group).edges;
            n_pos_bins = obj.bin_config.(group).n_bins;
            
            % Get trials for this group
            trials_struct = obj.trial_groups.(group);
            n_trials = length(trials_struct);
            
            % Collect all velocities during motion to compute equal-occupancy bins
            all_velocities = [];
            for t = 1:n_trials
                trial = trials_struct(t).trial;
                motion_mask = trial.motion_mask();
                velocity = trial.velocity(motion_mask);
                all_velocities = [all_velocities; velocity(:)];
            end
            
            % Compute equal-occupancy velocity bins using percentiles
            prc_per_bin = 100 / n_vel_bins;
            prc_bounds = 0 : prc_per_bin : 100+eps;
            vel_bin_edges = prctile(all_velocities, prc_bounds);
            
            % Initialize 2D arrays for each trial
            rate_maps_per_trial = nan(n_trials, n_pos_bins, n_vel_bins);
            occupancy_per_trial = nan(n_trials, n_pos_bins, n_vel_bins);
            bin_sampled_per_trial = zeros(n_trials, n_pos_bins, n_vel_bins);  % Track which bins are sampled per trial
            
            % Find cluster index in the stored data
            cluster_idx = find(obj.cluster_ids == cluster.id);
            if isempty(cluster_idx)
                error('Cluster %d not found in analyzer cluster_ids', cluster.id);
            end
            
            for t = 1:n_trials
                trial = trials_struct(t).trial;
                
                % Get the already-computed smoothed spatial firing rate for this trial
                % rate_per_trial_smooth is stored as a cell array {n_clusters x 1}
                % where each cell contains a matrix (n_trials x n_bins)
                rate_per_trial_smooth_cluster = obj.firing_rates.(group).rate_per_trial_smooth{cluster_idx};
                spatial_rate_smooth = rate_per_trial_smooth_cluster(t, :);
                
                % Get position and velocity during motion
                motion_mask = trial.motion_mask();
                tvec = trial.probe_t(motion_mask);
                position = trial.position(motion_mask);
                velocity = trial.velocity(motion_mask);
                
                % Initialize 2D rate map for this trial
                rate_map_trial = nan(n_pos_bins, n_vel_bins);
                
                % For each position bin, assign the spatial firing rate to the appropriate velocity bin
                for p = 1:n_pos_bins
                    % Find times when animal was in this position bin
                    in_pos_bin = position >= pos_edges(p) & position < pos_edges(p+1);
                    
                    if any(in_pos_bin)
                        % Get velocities when in this position bin
                        vel_in_bin = velocity(in_pos_bin);
                        
                        % Assign to velocity bins
                        for v = 1:n_vel_bins
                            in_vel_bin = vel_in_bin >= vel_bin_edges(v) & vel_in_bin < vel_bin_edges(v+1);
                            
                            % Mark as sampled if this position×velocity combination occurred
                            if any(in_vel_bin)
                                bin_sampled_per_trial(t, p, v) = 1;
                                % Assign the spatial firing rate (already smoothed)
                                rate_map_trial(p, v) = spatial_rate_smooth(p);
                            end
                        end
                    end
                end
                
                rate_maps_per_trial(t, :, :) = rate_map_trial;
            end
            
            % Compute median firing rate across trials
            rate_map_2d = squeeze(nanmedian(rate_maps_per_trial, 1));
            
            % Count how many trials sampled each bin
            trial_count_map = squeeze(sum(bin_sampled_per_trial, 1));
        end
        
        function [rate_map_2d, accel_bin_edges, trial_count_map] = compute_2d_position_acceleration_tuning(obj, cluster, group, n_accel_bins, sigma_accel)
            % Compute 2D position × acceleration firing rate map with equal-occupancy acceleration bins
            %
            %   [rate_map_2d, accel_bin_edges, trial_count_map] = analyzer.compute_2d_position_acceleration_tuning(cluster, group, n_accel_bins, sigma_accel)
            %
            %   Inputs:
            %       cluster         - Cluster object
            %       group           - 'long' or 'short'
            %       n_accel_bins    - Number of acceleration bins (will have equal occupancy)
            %       sigma_accel     - Smoothing sigma for acceleration dimension (cm/s²)
            %
            %   Outputs:
            %       rate_map_2d     - 2D firing rate map (position × acceleration), median across trials
            %       accel_bin_edges - Acceleration bin edges (variable width for equal occupancy)
            %       trial_count_map - Number of trials that sampled each bin (for alpha transparency)
            
            if nargin < 5
                sigma_accel = 10;  % Default smoothing sigma for acceleration (cm/s²) - reduced to preserve acceleration structure
            end
            
            % Get position bin configuration for this group
            pos_edges = obj.bin_config.(group).edges;
            n_pos_bins = obj.bin_config.(group).n_bins;
            
            % Get trials for this group
            trials_struct = obj.trial_groups.(group);
            n_trials = length(trials_struct);
            
            % Collect all accelerations during motion to compute equal-occupancy bins
            all_accelerations = [];
            for t = 1:n_trials
                trial = trials_struct(t).trial;
                motion_mask = trial.motion_mask();
                acceleration = trial.acceleration(motion_mask);
                all_accelerations = [all_accelerations; acceleration(:)];
            end
            
            % Compute equal-occupancy acceleration bins using percentiles
            prc_per_bin = 100 / n_accel_bins;
            prc_bounds = 0 : prc_per_bin : 100+eps;
            accel_bin_edges = prctile(all_accelerations, prc_bounds);
            
            % Initialize 2D arrays for each trial
            rate_maps_per_trial = nan(n_trials, n_pos_bins, n_accel_bins);
            occupancy_per_trial = nan(n_trials, n_pos_bins, n_accel_bins);
            bin_sampled_per_trial = zeros(n_trials, n_pos_bins, n_accel_bins);  % Track which bins are sampled per trial
            
            % Find cluster index in the stored data
            cluster_idx = find(obj.cluster_ids == cluster.id);
            if isempty(cluster_idx)
                error('Cluster %d not found in analyzer cluster_ids', cluster.id);
            end
            
            for t = 1:n_trials
                trial = trials_struct(t).trial;
                
                % Get the already-computed smoothed spatial firing rate for this trial
                % rate_per_trial_smooth is stored as a cell array {n_clusters x 1}
                % where each cell contains a matrix (n_trials x n_bins)
                rate_per_trial_smooth_cluster = obj.firing_rates.(group).rate_per_trial_smooth{cluster_idx};
                spatial_rate_smooth = rate_per_trial_smooth_cluster(t, :);
                
                % Get position and acceleration during motion
                motion_mask = trial.motion_mask();
                tvec = trial.probe_t(motion_mask);
                position = trial.position(motion_mask);
                acceleration = trial.acceleration(motion_mask);
                
                % Initialize 2D rate map for this trial
                rate_map_trial = nan(n_pos_bins, n_accel_bins);
                
                % For each position bin, assign the spatial firing rate to the appropriate acceleration bin
                for p = 1:n_pos_bins
                    % Find times when animal was in this position bin
                    in_pos_bin = position >= pos_edges(p) & position < pos_edges(p+1);
                    
                    if any(in_pos_bin)
                        % Get accelerations when in this position bin
                        accel_in_bin = acceleration(in_pos_bin);
                        
                        % Assign to acceleration bins
                        for a = 1:n_accel_bins
                            in_accel_bin = accel_in_bin >= accel_bin_edges(a) & accel_in_bin < accel_bin_edges(a+1);
                            
                            % Mark as sampled if this position×acceleration combination occurred
                            if any(in_accel_bin)
                                bin_sampled_per_trial(t, p, a) = 1;
                                % Assign the spatial firing rate (already smoothed)
                                rate_map_trial(p, a) = spatial_rate_smooth(p);
                            end
                        end
                    end
                end
                
                rate_maps_per_trial(t, :, :) = rate_map_trial;
            end
            
            % Compute median firing rate across trials
            rate_map_2d = squeeze(nanmedian(rate_maps_per_trial, 1));
            
            % Count how many trials sampled each bin
            trial_count_map = squeeze(sum(bin_sampled_per_trial, 1));
        end
        
        function compute_ttg_analysis(obj)
            % Compute time-to-goal analysis for all clusters
            %
            %   analyzer.compute_ttg_analysis()
            %
            %   Computes TTG tuning and statistics for each cluster.
            %   Results are stored in obj.ttg_analysis.
            
            if isempty(obj.trial_groups)
                error('Trial groups not initialized. Run collect_trials() first.');
            end
            
            n_clusters = length(obj.clusters);
            fprintf('  Computing TTG analysis for %d clusters...\n', n_clusters);
            
            % Initialize storage
            obj.ttg_analysis = struct();
            
            % Process each cluster
            for c = 1:n_clusters
                cluster = obj.clusters(c);
                cluster_id = cluster.id;
                cluster_field = sprintf('cluster_%d', cluster_id);
                
                % Get global trial indices for this cluster (wrap in cell for indexing)
                all_global_trial_indices_by_group = struct();
                for g = 1:2
                    group = obj.group_names{g};
                    all_global_trial_indices_by_group.(group) = ...
                        {obj.firing_rates.(group).global_trial_indices{c}};
                end
                
                % Compute TTG data (use cluster_idx = 1 since we wrapped data in cell)
                [ttg_data, ttg_norm_bin_centers] = compute_time_to_goal_tuning(...
                    cluster, obj.trial_groups, obj.group_names, ...
                    all_global_trial_indices_by_group, 1);
                
                % Compute TTG statistics
                ttg_stats = [];
                if obj.stats_params.run_statistics
                    ttg_stats = compute_ttg_tuning_statistics(...
                        ttg_data, obj.group_names, obj.stats_params);
                end
                
                % Store results
                obj.ttg_analysis.(cluster_field) = struct(...
                    'ttg_data', ttg_data, ...
                    'ttg_norm_bin_centers', ttg_norm_bin_centers, ...
                    'ttg_stats', ttg_stats);
            end
            
            fprintf('  TTG analysis complete\n');
        end
        
        function compute_distribution_comparisons(obj)
            % Compute Kruskal-Wallis distribution comparisons for all clusters
            %
            %   analyzer.compute_distribution_comparisons()
            %
            %   Compares spatial distributions between long and short trials
            %   using Kruskal-Wallis test. Results stored in obj.distribution_comparisons.
            
            if isempty(obj.firing_rates)
                error('Firing rates not computed. Run analyze_all_clusters() first.');
            end
            
            if isempty(obj.ttg_analysis)
                error('TTG analysis not computed. Run compute_ttg_analysis() first.');
            end
            
            n_clusters = length(obj.clusters);
            fprintf('  Computing distribution comparisons for %d clusters...\n', n_clusters);
            
            % Initialize storage
            obj.distribution_comparisons = cell(n_clusters, 1);
            
            % Counter for successful comparisons
            n_spatially_tuned = 0;
            n_successful = 0;
            
            % Process each cluster
            for c = 1:n_clusters
                cluster = obj.clusters(c);
                cluster_id = cluster.id;
                cluster_field = sprintf('cluster_%d', cluster_id);
                
                % Check if spatially tuned
                is_spatially_tuned = false;
                if ~isempty(obj.tuning_stats) && isfield(obj.tuning_stats, cluster_field)
                    cluster_stats = obj.tuning_stats.(cluster_field);
                    if (isfield(cluster_stats, 'long') && cluster_stats.long.isSpatiallyTuned) || ...
                       (isfield(cluster_stats, 'short') && cluster_stats.short.isSpatiallyTuned)
                        is_spatially_tuned = true;
                        n_spatially_tuned = n_spatially_tuned + 1;
                    end
                end
                
                % Only compare if spatially tuned
                if is_spatially_tuned
                    try
                        % Get rate data
                        rate_per_trial_long = obj.firing_rates.long.rate_per_trial{c};
                        rate_per_trial_short = obj.firing_rates.short.rate_per_trial{c};
                        bin_centers_long = obj.bin_config.long.centers;
                        bin_centers_short = obj.bin_config.short.centers;
                        
                        % Get TTG rates
                        ttg_data = obj.ttg_analysis.(cluster_field).ttg_data;
                        ttg_rate_long = ttg_data.long.rate_per_trial_smooth;
                        ttg_rate_short = ttg_data.short.rate_per_trial_smooth;
                        
                        % Compute downsampled long data
                        [rate_long_ds_smooth, ~, ~, ~, bin_centers_long_ds, edges_long_ds] = ...
                            compute_downsampled_long_rates(obj.trial_groups, cluster, ...
                                obj.bin_size_cm, obj.gauss_sigma_cm);
                        
                        % Perform comparison
                        results = compare_spatial_distributions(...
                            rate_per_trial_long, rate_per_trial_short, ...
                            bin_centers_long, bin_centers_short, ...
                            ttg_rate_long, ttg_rate_short, ...
                            rate_long_ds_smooth, bin_centers_long_ds);
                        
                        obj.distribution_comparisons{c} = results;
                        n_successful = n_successful + 1;
                    catch ME
                        warning('Failed to compute distribution comparison for cluster %d: %s', ...
                            cluster_id, ME.message);
                        obj.distribution_comparisons{c} = [];
                    end
                else
                    obj.distribution_comparisons{c} = [];
                end
            end
            
            fprintf('  Distribution comparisons complete: %d spatially tuned, %d successful comparisons\n', ...
                n_spatially_tuned, n_successful);
        end
        
        function save_cache(obj, cache_filepath)
            % Save analysis results to cache file
            %
            %   analyzer.save_cache(cache_filepath)
            
            % Convert to old format for compatibility with plotting functions
            all_rate_smooth_by_group = struct();
            all_bin_centers_by_group = struct();
            all_avg_velocity_by_group = struct();
            all_Q1_rate_smooth_by_group = struct();
            all_Q2_rate_smooth_by_group = struct();
            all_Q3_rate_smooth_by_group = struct();
            all_occ_by_group = struct();
            all_rate_per_trial_by_group = struct();
            all_rate_per_trial_smooth_by_group = struct();
            all_spike_positions_by_group = struct();
            all_global_trial_indices_by_group = struct();
            
            for g = 1:2
                group = obj.group_names{g};
                all_rate_smooth_by_group.(group) = obj.firing_rates.(group).rate_smooth;
                all_bin_centers_by_group.(group) = obj.bin_config.(group).plot_centers;
                all_avg_velocity_by_group.(group) = obj.firing_rates.(group).avg_velocity;
                all_Q1_rate_smooth_by_group.(group) = obj.firing_rates.(group).Q1_smooth;
                all_Q2_rate_smooth_by_group.(group) = obj.firing_rates.(group).Q2_smooth;
                all_Q3_rate_smooth_by_group.(group) = obj.firing_rates.(group).Q3_smooth;
                all_occ_by_group.(group) = obj.firing_rates.(group).occupancy;
                all_rate_per_trial_by_group.(group) = obj.firing_rates.(group).rate_per_trial;
                all_rate_per_trial_smooth_by_group.(group) = obj.firing_rates.(group).rate_per_trial_smooth;
                all_spike_positions_by_group.(group) = obj.firing_rates.(group).spike_positions;
                all_global_trial_indices_by_group.(group) = obj.firing_rates.(group).global_trial_indices;
            end
            
            spatial_tuning_stats = obj.tuning_stats;
            cluster_ids = obj.cluster_ids;
            n_clusters = length(obj.clusters);
            group_names = obj.group_names;
            group_labels = obj.group_labels;
            bin_size_cm = obj.bin_size_cm;
            gauss_sigma_cm = obj.gauss_sigma_cm;
            
            % New cached data
            ttg_analysis = obj.ttg_analysis;
            distribution_comparisons = obj.distribution_comparisons;
            
            save(cache_filepath, 'all_rate_smooth_by_group', ...
                 'all_bin_centers_by_group', 'all_avg_velocity_by_group', ...
                 'all_Q1_rate_smooth_by_group', 'all_Q2_rate_smooth_by_group', ...
                 'all_Q3_rate_smooth_by_group', 'all_occ_by_group', ...
                 'all_rate_per_trial_by_group', 'all_rate_per_trial_smooth_by_group', ...
                 'all_spike_positions_by_group', 'all_global_trial_indices_by_group', ...
                 'spatial_tuning_stats', ...
                 'cluster_ids', 'n_clusters', 'group_names', 'group_labels', ...
                 'bin_size_cm', 'gauss_sigma_cm', ...
                 'ttg_analysis', 'distribution_comparisons', '-v7.3');
        end
        
        function load_cache(obj, cache_filepath)
            % Load analysis results from cache file
            %
            %   analyzer.load_cache(cache_filepath)
            
            loaded_vars = load(cache_filepath);
            
            % Load core variables
            for g = 1:2
                group = obj.group_names{g};
                obj.firing_rates.(group).rate_smooth = loaded_vars.all_rate_smooth_by_group.(group);
                obj.firing_rates.(group).Q1_smooth = loaded_vars.all_Q1_rate_smooth_by_group.(group);
                obj.firing_rates.(group).Q2_smooth = loaded_vars.all_Q2_rate_smooth_by_group.(group);
                obj.firing_rates.(group).Q3_smooth = loaded_vars.all_Q3_rate_smooth_by_group.(group);
                obj.firing_rates.(group).avg_velocity = loaded_vars.all_avg_velocity_by_group.(group);
                obj.firing_rates.(group).occupancy = loaded_vars.all_occ_by_group.(group);
                obj.firing_rates.(group).rate_per_trial = loaded_vars.all_rate_per_trial_by_group.(group);
                obj.firing_rates.(group).rate_per_trial_smooth = loaded_vars.all_rate_per_trial_smooth_by_group.(group);
                obj.firing_rates.(group).spike_positions = loaded_vars.all_spike_positions_by_group.(group);
                obj.firing_rates.(group).global_trial_indices = loaded_vars.all_global_trial_indices_by_group.(group);
            end
            
            obj.cluster_ids = loaded_vars.cluster_ids;
            obj.tuning_stats = loaded_vars.spatial_tuning_stats;
            obj.bin_size_cm = loaded_vars.bin_size_cm;
            obj.gauss_sigma_cm = loaded_vars.gauss_sigma_cm;
            
            % Load TTG analysis and distribution comparisons
            obj.ttg_analysis = loaded_vars.ttg_analysis;
            obj.distribution_comparisons = loaded_vars.distribution_comparisons;
            
            % Reconstruct bin configuration from loaded parameters
            obj.setup_bins();
        end
        
        function setup_parallel_pool(obj)
            % Initialize or verify parallel pool
            current_pool = gcp('nocreate');
            if isempty(current_pool)
                fprintf('  Starting parallel pool with %d workers...\n', obj.max_workers);
                parpool('local', obj.max_workers);
            elseif current_pool.NumWorkers ~= obj.max_workers
                fprintf('  Restarting parallel pool with %d workers (was %d)...\n', ...
                    obj.max_workers, current_pool.NumWorkers);
                delete(current_pool);
                parpool('local', obj.max_workers);
            else
                fprintf('  Using existing parallel pool with %d workers\n', current_pool.NumWorkers);
            end
        end
        
        function status = parallel_status(obj)
            % Get string describing parallel processing status
            if obj.use_parallel
                status = sprintf('parallelized, %d workers', obj.max_workers);
            else
                status = 'sequential';
            end
        end
        
        function print_summary(obj)
            % Print summary of analysis results
            fprintf('\n=== Spatial Tuning Analysis Summary ===\n');
            fprintf('Clusters analyzed: %d\n', length(obj.clusters));
            fprintf('Bin size: %.1f cm, Smoothing sigma: %.1f cm\n', ...
                obj.bin_size_cm, obj.gauss_sigma_cm);
            
            if ~isempty(obj.tuning_stats)
                n_tuned_long = 0;
                n_tuned_short = 0;
                n_tuned_both = 0;
                
                cluster_fields = fieldnames(obj.tuning_stats);
                for i = 1:length(cluster_fields)
                    stats = obj.tuning_stats.(cluster_fields{i});
                    tuned_long = stats.long.isSpatiallyTuned;
                    tuned_short = stats.short.isSpatiallyTuned;
                    
                    if tuned_long
                        n_tuned_long = n_tuned_long + 1;
                    end
                    if tuned_short
                        n_tuned_short = n_tuned_short + 1;
                    end
                    if tuned_long && tuned_short
                        n_tuned_both = n_tuned_both + 1;
                    end
                end
                
                fprintf('Spatially tuned in long trials: %d (%.1f%%)\n', ...
                    n_tuned_long, 100*n_tuned_long/length(cluster_fields));
                fprintf('Spatially tuned in short trials: %d (%.1f%%)\n', ...
                    n_tuned_short, 100*n_tuned_short/length(cluster_fields));
                fprintf('Spatially tuned in both: %d (%.1f%%)\n', ...
                    n_tuned_both, 100*n_tuned_both/length(cluster_fields));
            end
            fprintf('=======================================\n\n');
        end
    end
end
