% GLM cross-profile validation for speed generalization
%
% This script performs an out-of-distribution (OOD) validation where the GLM
% is trained on one speed profile and tested on the other. It mirrors the
% core data flow of glm_single_cluster_analysis.m while changing the CV split:
%   - Standard CV: trial-level condition-stratified K-fold
%   - Cross-profile CV: fold 1 = profile A, fold 2 = profile B

%% ====================================================================
%  Section 1: Configuration
%  ====================================================================
close all; clearvars;

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
    fprintf('Started parallel pool with %d workers\n', pool.NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', pool.NumWorkers);
end

log_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
log_filename = sprintf('glm_cross_profile_log_%s.txt', log_timestamp);
log_file_tmp = fullfile(tempdir, log_filename);
diary(log_file_tmp);
fprintf('Cross-profile analysis started: %s\n\n', datestr(now));

experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};
restricted              = get_restricted_flag();

time_bin_width          = 0.02;     % seconds
n_cv_folds              = 5;
n_speed_bases           = 5;
n_tf_bases              = 5;
n_onset_bases           = 6;
onset_range             = [0, 2.0];
min_trials_per_profile  = 2;
n_profiles              = 2;
n_profile_samples       = 80;       % for profile-shape extraction

% Quick-run selectors (set [] to disable each filter)
quick_probe_indices     = [];        % e.g., 1 for first probe, [1 3] for multiple
quick_max_trials        = [];       % max trials per selected probe
quick_max_clusters      = [];        % max VISp clusters per selected probe
quick_cluster_ids       = [];       % explicit cluster IDs override quick_max_clusters
quick_random_seed       = 1;        % reproducible trial subsampling

save_figs               = true;
figure_dir              = {'glm_cross_profile_validation'};

tuning_curves_dir       = 'D:\mvelez\formatted_data\csvs\tuning_curves';
tf_tuning_curves_dir    = 'D:\mvelez\formatted_data\csvs\tf_tuning_curves';
mc_sequence_path        = 'D:\mvelez\mateoData_mc\motion_cloud_sequence_250414.mat';
mc_folders_path         = 'D:\mvelez\mateoData_mc\image_folders.mat';

exclude_patterns = {'theta0p000_Btheta3p142_sf00p006_Bsf0p004_VX0p000_BV2p000'};

sf_values_map = containers.Map(...
    {'sf00p003', 'sf00p006', 'sf00p012'}, ...
    {0.003, 0.006, 0.012});
sf_keys = keys(sf_values_map);

or_values_map = containers.Map(...
    {'theta0p000', 'theta0p785', 'theta1p571', 'theta-0p785'}, ...
    {0, pi/4, pi/2, -pi/4});
or_keys = keys(or_values_map);

batch_patterns = { ...
    {'sf00p003_Bsf0p002_VX1p002', 'sf00p006_Bsf0p002_VX0p501', 'sf00p012_Bsf0p002_VX0p250'}, ...
    {'sf00p003_Bsf0p002_VX2p003', 'sf00p006_Bsf0p002_VX1p002', 'sf00p012_Bsf0p002_VX0p501'}, ...
    {'sf00p003_Bsf0p002_VX4p006', 'sf00p006_Bsf0p002_VX2p003', 'sf00p012_Bsf0p002_VX1p002'} };
batch_gains = [1/30, 2/30, 4/30];

fprintf('Configuration:\n');
fprintf('  Time-bin width: %.0f ms\n', time_bin_width*1000);
fprintf('  Standard CV folds: %d\n', n_cv_folds);
fprintf('  Profiles in cross-profile CV: %d\n', n_profiles);
fprintf('  Speed bases: %d | TF bases: %d | Onset bases: %d\n', n_speed_bases, n_tf_bases, n_onset_bases);
fprintf('  Restricted mode: %d\n\n', restricted);
fprintf('  Quick probe indices: %s\n', mat2str(quick_probe_indices));
fprintf('  Quick max trials/probe: %s\n', mat2str(quick_max_trials));
fprintf('  Quick max clusters/probe: %s\n\n', mat2str(quick_max_clusters));
fprintf('  Quick random seed: %d\n\n', quick_random_seed);

%% ====================================================================
%  Section 2: Metadata and controllers
%  ====================================================================
fprintf('--- Loading motion cloud metadata ---\n');

[mc_sequence, cloud_names] = glm_helpers.load_motion_cloud_metadata(mc_sequence_path, mc_folders_path);

if isempty(mc_sequence) || isempty(cloud_names)
    error('Could not load motion cloud metadata (presentation sequence or cloud names).');
end

fprintf('  Loaded %d sequence entries, %d cloud names\n', length(mc_sequence), length(cloud_names));

ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});

if ~isempty(quick_probe_indices)
    probe_idx = quick_probe_indices(:)';
    probe_idx = probe_idx(probe_idx >= 1 & probe_idx <= numel(probe_ids));
    probe_idx = unique(probe_idx, 'stable');
    probe_ids = probe_ids(probe_idx);
end

if isempty(probe_ids)
    error('No probes selected after applying quick_probe_indices.');
end

ctl.setup_figures(figure_dir, save_figs);

diary off;
log_file = fullfile(ctl.figs.curr_dir, log_filename);
movefile(log_file_tmp, log_file);
diary(log_file);
fprintf('Log file: %s\n', log_file);

%% ====================================================================
%  Section 3: Trial stimulus lookup and speed/TF ranges
%  ====================================================================
fprintf('\n--- Building trial stimulus lookup ---\n');

trial_stim = glm_helpers.build_trial_stim_lookup(mc_sequence, cloud_names, exclude_patterns, ...
    sf_keys, sf_values_map, or_keys, or_values_map, batch_patterns, batch_gains);

if isempty(probe_ids)
    error('No probes found for experiment group(s): %s', strjoin(experiment_groups, ', '));
end

first_spd_file = fullfile(tuning_curves_dir, [probe_ids{1} '.mat']);
first_tf_file  = fullfile(tf_tuning_curves_dir, [probe_ids{1} '.mat']);
if ~exist(first_spd_file, 'file') || ~exist(first_tf_file, 'file')
    error('Could not read tuning tables for the first probe (%s).', probe_ids{1});
end

D_tmp = load(first_spd_file);
spd_bin_edges = D_tmp.tuning_curves{1}(1).bin_edges(:)';
D_tmp = load(first_tf_file);
tf_bin_edges  = D_tmp.tuning_curves{1}(1).bin_edges(:)';

speed_range = [spd_bin_edges(1), spd_bin_edges(end)];
tf_range    = [tf_bin_edges(1), tf_bin_edges(end)];
fprintf('  Speed range: [%.2f, %.2f] cm/s\n', speed_range(1), speed_range(2));
fprintf('  TF range: [%.4f, %.4f] Hz\n', tf_range(1), tf_range(2));

%% ====================================================================
%  Section 4: Build time-bin master table with profile IDs
%  ====================================================================
fprintf('\n--- Building time-bin data with profile assignments [PARALLEL] ---\n');

all_probe_tables = cell(length(probe_ids), 1);

for probe_i = 1:length(probe_ids)
    pid = probe_ids{probe_i};
    fprintf('\nProbe %d/%d: %s\n', probe_i, length(probe_ids), pid);

    data = ctl.load_formatted_data(pid);
    visp_clusters = data.VISp_clusters();
    if isempty(visp_clusters)
        fprintf('  No VISp clusters, skipping\n');
        continue;
    end
    visp_ids = arrayfun(@(x) x.id, visp_clusters);

    if restricted
        keep_restr = true(size(visp_ids));
        for ci = 1:length(visp_ids)
            [~, ~, dir_T] = data.is_stationary_vs_motion_significant(visp_ids(ci), 'T_Vstatic');
            [~, ~, dir_V] = data.is_stationary_vs_motion_significant(visp_ids(ci), 'V');
            if dir_T == 0 || dir_V == 0
                keep_restr(ci) = false;
            end
        end
        visp_ids = visp_ids(keep_restr);
    end
    if isempty(visp_ids)
        fprintf('  No VISp clusters after restriction, skipping\n');
        continue;
    end

    if ~isempty(quick_cluster_ids)
        visp_ids = intersect(visp_ids, quick_cluster_ids, 'stable');
    elseif ~isempty(quick_max_clusters)
        visp_ids = visp_ids(1:min(length(visp_ids), quick_max_clusters));
    end

    if isempty(visp_ids)
        fprintf('  No VISp clusters after quick-run cluster filter, skipping\n');
        continue;
    end

    all_trial_ids_map = containers.Map('KeyType', 'int32', 'ValueType', 'char');
    for cond_i = 1:length(trial_group_labels)
        condition = trial_group_labels{cond_i};
        trials_cond = data.get_trials_with_trial_group_label(condition);
        for ti_k = 1:length(trials_cond)
            tid_k = trials_cond{ti_k}.trial_id;
            if tid_k >= 1 && tid_k <= length(trial_stim) && ~trial_stim(tid_k).excluded
                all_trial_ids_map(int32(tid_k)) = condition;
            end
        end
    end

    all_trial_ids = cell2mat(keys(all_trial_ids_map));
    all_trial_conditions = values(all_trial_ids_map, num2cell(all_trial_ids));

    if ~isempty(quick_max_trials) && numel(all_trial_ids) > quick_max_trials
        rng(quick_random_seed + probe_i, 'twister');
        trial_sel = randperm(numel(all_trial_ids), quick_max_trials);
        all_trial_ids = all_trial_ids(trial_sel);
        all_trial_conditions = all_trial_conditions(trial_sel);
    end

    n_trials_probe = length(all_trial_ids);
    if n_trials_probe == 0
        fprintf('  No valid trials after exclusion, skipping\n');
        continue;
    end
    fprintf('  Trials to process: %d\n', n_trials_probe);

    trial_objs = cell(n_trials_probe, 1);
    for ti_pre = 1:n_trials_probe
        trial_id = all_trial_ids(ti_pre);
        try
            trial_obj = data.get_trials_with_trial_ids(trial_id);
            if iscell(trial_obj)
                trial_obj = trial_obj{1};
            end
            trial_objs{ti_pre} = trial_obj.to_aligned;
        catch
            trial_objs{ti_pre} = [];
        end
    end

    % Extract per-trial speed templates and cluster into profile IDs
    trial_templates = nan(n_trials_probe, n_profile_samples);
    valid_template = false(n_trials_probe, 1);
    for ti = 1:n_trials_probe
        aligned_obj = trial_objs{ti};
        if isempty(aligned_obj)
            continue;
        end
        [ok, template] = glm_helpers.get_trial_speed_template(aligned_obj, n_profile_samples);
        if ok
            trial_templates(ti, :) = template;
            valid_template(ti) = true;
        end
    end

    trial_profile_ids = zeros(n_trials_probe, 1);
    if sum(valid_template) >= n_profiles
        X_prof = trial_templates(valid_template, :);
        X_prof = X_prof(:, any(~isnan(X_prof), 1));

        if size(X_prof, 2) < 2
            fprintf('  Not enough profile-feature dimensions, skipping probe\n');
            continue;
        end

        if any(isnan(X_prof(:)))
            col_means = mean(X_prof, 1, 'omitnan');
            for ci_col = 1:size(X_prof, 2)
                bad = isnan(X_prof(:, ci_col));
                X_prof(bad, ci_col) = col_means(ci_col);
            end
        end

        try
            [cluster_idx, centroids] = kmeans(X_prof, n_profiles, 'Replicates', 10, 'MaxIter', 500, 'Display', 'off');
        catch
            avg_speed = mean(X_prof, 2, 'omitnan');
            thr = median(avg_speed, 'omitnan');
            cluster_idx = ones(size(avg_speed));
            cluster_idx(avg_speed > thr) = 2;
            centroids = [mean(X_prof(cluster_idx == 1, :), 1, 'omitnan'); mean(X_prof(cluster_idx == 2, :), 1, 'omitnan')];
        end

        centroid_speed = mean(centroids, 2, 'omitnan');
        [~, ord] = sort(centroid_speed, 'ascend');
        relabel = zeros(n_profiles, 1);
        for ii = 1:n_profiles
            relabel(ord(ii)) = ii;
        end
        cluster_idx = arrayfun(@(x) relabel(x), cluster_idx);

        trial_profile_ids(valid_template) = cluster_idx;
    end

    n_p1 = sum(trial_profile_ids == 1);
    n_p2 = sum(trial_profile_ids == 2);
    fprintf('  Assigned profile IDs: P1=%d, P2=%d, unassigned=%d\n', n_p1, n_p2, sum(trial_profile_ids == 0));

    if n_p1 < min_trials_per_profile || n_p2 < min_trials_per_profile
        fprintf('  Insufficient profile coverage for this probe, skipping\n');
        continue;
    end

    n_visp = length(visp_ids);
    spike_times_all = cell(n_visp, 1);
    valid_cids_mask = false(n_visp, 1);
    for ci_cache = 1:n_visp
        cid_cache = visp_ids(ci_cache);
        try
            cluster_obj = data.get_cluster_with_id(cid_cache);
            spike_times_all{ci_cache} = cluster_obj.spike_times;
            valid_cids_mask(ci_cache) = true;
        catch
        end
    end
    valid_cids = visp_ids(valid_cids_mask);
    spike_times_valid = spike_times_all(valid_cids_mask);
    n_valid_clusters = length(valid_cids);
    if n_valid_clusters == 0
        fprintf('  No valid clusters with spike_times, skipping\n');
        continue;
    end
    fprintf('  Spike times cached for %d clusters\n', n_valid_clusters);

    trial_sf_vals = zeros(n_trials_probe, 1);
    trial_gain_vals = zeros(n_trials_probe, 1);
    trial_or_vals = zeros(n_trials_probe, 1);
    for ti_pre = 1:n_trials_probe
        stim = trial_stim(all_trial_ids(ti_pre));
        trial_sf_vals(ti_pre) = stim.sf;
        trial_gain_vals(ti_pre) = stim.batch_gain;
        trial_or_vals(ti_pre) = stim.or;
    end

    trial_results = cell(n_trials_probe, 1);

    parfor ti_main = 1:n_trials_probe
        trial_id = all_trial_ids(ti_main);
        profile_id = trial_profile_ids(ti_main);
        if profile_id < 1
            continue;
        end

        condition = all_trial_conditions{ti_main};
        sf_val = trial_sf_vals(ti_main);
        gain_val = trial_gain_vals(ti_main);
        or_val = trial_or_vals(ti_main);

        aligned_obj = trial_objs{ti_main};
        if isempty(aligned_obj)
            continue;
        end

        tr_probe_t = aligned_obj.probe_t;
        tr_vel     = aligned_obj.velocity;
        tr_mmask   = aligned_obj.motion_mask;
        tr_smask   = aligned_obj.stationary_mask;

        motion_idx = find(tr_mmask);
        if isempty(motion_idx)
            continue;
        end
        motion_start_idx = motion_idx(1);
        motion_end_idx   = motion_idx(end);

        t_motion_start = tr_probe_t(motion_start_idx);
        t_motion_end   = tr_probe_t(motion_end_idx);
        motion_dur = t_motion_end - t_motion_start;
        if motion_dur < time_bin_width
            continue;
        end

        bin_edges_t = t_motion_start : time_bin_width : t_motion_end;
        n_tbins = length(bin_edges_t) - 1;
        if n_tbins < 1
            continue;
        end

        sample_bin = discretize(tr_probe_t, bin_edges_t);
        sample_bin(isnan(sample_bin)) = 0;

        valid_mask   = sample_bin > 0;
        bins_vec     = sample_bin(valid_mask);
        mmask_vec    = tr_mmask(valid_mask);
        absvel_vec   = abs(tr_vel(valid_mask)) .* mmask_vec;

        n_samp_per_bin    = accumarray(bins_vec(:), 1, [n_tbins 1]);
        n_motion_per_bin  = accumarray(bins_vec(:), mmask_vec(:), [n_tbins 1]);
        sum_speed_per_bin = accumarray(bins_vec(:), absvel_vec(:), [n_tbins 1]);

        good_bins = find(n_samp_per_bin > 0 & n_motion_per_bin >= n_samp_per_bin * 0.5);
        if isempty(good_bins)
            continue;
        end

        mean_speed_per_bin = sum_speed_per_bin(good_bins) ./ n_motion_per_bin(good_bins);
        bin_centres = (bin_edges_t(1:end-1) + bin_edges_t(2:end)) / 2;
        time_in_trial_vec = bin_centres(good_bins)' - t_motion_start;

        n_good = length(good_bins);
        switch condition
            case 'T_Vstatic'
                speed_vec = mean_speed_per_bin;
                tf_vec    = zeros(n_good, 1);
                sf_vec    = NaN(n_good, 1);
                or_vec    = NaN(n_good, 1);
                gain_vec  = NaN(n_good, 1);
            case 'V'
                speed_vec = zeros(n_good, 1);
                tf_vec    = gain_val * mean_speed_per_bin;
                sf_vec    = repmat(sf_val, n_good, 1);
                or_vec    = repmat(or_val, n_good, 1);
                gain_vec  = repmat(gain_val, n_good, 1);
            case 'VT'
                speed_vec = mean_speed_per_bin;
                tf_vec    = gain_val * mean_speed_per_bin;
                sf_vec    = repmat(sf_val, n_good, 1);
                or_vec    = repmat(or_val, n_good, 1);
                gain_vec  = repmat(gain_val, n_good, 1);
            otherwise
                speed_vec = mean_speed_per_bin;
                tf_vec    = zeros(n_good, 1);
                sf_vec    = NaN(n_good, 1);
                or_vec    = NaN(n_good, 1);
                gain_vec  = NaN(n_good, 1);
        end

        motion_rows = cell(n_valid_clusters, 1);
        for ci = 1:n_valid_clusters
            st = spike_times_valid{ci};
            spike_counts_all = histcounts(st, bin_edges_t)';
            spike_counts_good = spike_counts_all(good_bins);

            motion_rows{ci} = [repmat(valid_cids(ci), n_good, 1), ...
                               repmat(trial_id, n_good, 1), ...
                               repmat(profile_id, n_good, 1), ...
                               speed_vec, tf_vec, sf_vec, or_vec, gain_vec, ...
                               spike_counts_good, time_in_trial_vec, time_in_trial_vec];
        end
        motion_data = vertcat(motion_rows{:});
        motion_conds = repmat({condition}, size(motion_data, 1), 1);

        stat_data = [];
        stat_conds = {};
        stationary_before_motion = tr_smask & (1:length(tr_smask))' < motion_start_idx;
        stationary_idx = find(stationary_before_motion);

        if ~isempty(stationary_idx)
            stat_end_idx = stationary_idx(end);
            gaps = find(diff(stationary_idx) > 1);
            if ~isempty(gaps)
                stat_start_idx = stationary_idx(gaps(end) + 1);
            else
                stat_start_idx = stationary_idx(1);
            end

            t_stat_start = tr_probe_t(stat_start_idx);
            t_stat_end   = tr_probe_t(stat_end_idx);
            stat_dur = t_stat_end - t_stat_start;

            if stat_dur >= time_bin_width
                bin_edges_stat = t_stat_start : time_bin_width : t_stat_end;
                n_stat_bins = length(bin_edges_stat) - 1;

                if n_stat_bins >= 1
                    n_stat = n_stat_bins;
                    bin_centres_stat = (bin_edges_stat(1:end-1) + bin_edges_stat(2:end)) / 2;
                    time_since_onset_stat = bin_centres_stat' - t_motion_start;

                    stat_rows = cell(n_valid_clusters, 1);
                    for ci = 1:n_valid_clusters
                        st = spike_times_valid{ci};
                        spike_counts_stat = histcounts(st, bin_edges_stat)';

                        stat_rows{ci} = [repmat(valid_cids(ci), n_stat, 1), ...
                                         repmat(trial_id, n_stat, 1), ...
                                         repmat(profile_id, n_stat, 1), ...
                                         zeros(n_stat, 1), zeros(n_stat, 1), ...
                                         NaN(n_stat, 1), NaN(n_stat, 1), NaN(n_stat, 1), ...
                                         spike_counts_stat, zeros(n_stat, 1), time_since_onset_stat];
                    end
                    stat_data = vertcat(stat_rows{:});
                    stat_conds = repmat({'stationary'}, size(stat_data, 1), 1);
                end
            end
        end

        all_data = [motion_data; stat_data];
        all_conds = [motion_conds; stat_conds];
        trial_results{ti_main} = struct('data', all_data, 'conds', {all_conds}, 'pid', pid);
    end

    probe_data = [];
    probe_conds = {};
    probe_pids = {};
    for ti_main = 1:n_trials_probe
        if ~isempty(trial_results{ti_main}) && ~isempty(trial_results{ti_main}.data)
            probe_data = [probe_data; trial_results{ti_main}.data]; %#ok<AGROW>
            probe_conds = [probe_conds; trial_results{ti_main}.conds]; %#ok<AGROW>
            probe_pids = [probe_pids; repmat({pid}, size(trial_results{ti_main}.data, 1), 1)]; %#ok<AGROW>
        end
    end

    if ~isempty(probe_data)
        all_probe_tables{probe_i} = table(...
            string(probe_pids), ...
            probe_data(:,1), probe_data(:,2), probe_data(:,3), string(probe_conds), ...
            probe_data(:,4), probe_data(:,5), probe_data(:,6), probe_data(:,7), probe_data(:,8), ...
            probe_data(:,9), probe_data(:,10), probe_data(:,11), ...
            'VariableNames', {'probe_id', 'cluster_id', 'trial_id', 'profile_id', 'condition', ...
                'speed', 'tf', 'sf', 'orientation', 'batch_gain', ...
                'spike_count', 'time_in_trial', 'time_since_onset'});
        fprintf('  Probe rows: %d\n', height(all_probe_tables{probe_i}));
    else
        fprintf('  Probe produced no rows\n');
    end
end

valid_tables = all_probe_tables(~cellfun(@isempty, all_probe_tables));
if isempty(valid_tables)
    error('No probe generated usable time-bin data with profile assignments.');
end

T_master_time = vertcat(valid_tables{:});
fprintf('\nTime-bin table built: %d rows x %d columns\n', height(T_master_time), width(T_master_time));

%% ====================================================================
%  Section 5: Cluster-level GLM metrics (Null vs M0_Speed)
%  ====================================================================
fprintf('\n--- Fitting per-cluster GLMs (Null vs M0_Speed) ---\n');

unique_clusters = unique(T_master_time(:, {'probe_id', 'cluster_id'}), 'rows');
n_clusters = height(unique_clusters);
fprintf('Clusters to evaluate: %d\n', n_clusters);

results = table('Size', [n_clusters, 14], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double','double','double','double','double','double','logical'}, ...
    'VariableNames', {'probe_id','cluster_id','n_rows','n_trials','n_trials_p1','n_trials_p2', ...
                      'null_bps_std','speed_bps_std','delta_bps_std', ...
                      'null_bps_cross','speed_bps_cross','delta_bps_cross', ...
                      'delta_drop','usable'});

for ci = 1:n_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);

    idx = T_master_time.probe_id == pid & T_master_time.cluster_id == cid & T_master_time.profile_id > 0;
    T_cluster = T_master_time(idx, :);

    results.probe_id(ci) = pid;
    results.cluster_id(ci) = cid;
    results.n_rows(ci) = height(T_cluster);
    results.usable(ci) = false;

    if height(T_cluster) < 20
        continue;
    end

    tr_unique = unique(T_cluster(:, {'trial_id', 'profile_id'}), 'rows');
    n_trials_total = height(tr_unique);
    n_trials_p1 = sum(tr_unique.profile_id == 1);
    n_trials_p2 = sum(tr_unique.profile_id == 2);

    results.n_trials(ci) = n_trials_total;
    results.n_trials_p1(ci) = n_trials_p1;
    results.n_trials_p2(ci) = n_trials_p2;

    if n_trials_p1 < min_trials_per_profile || n_trials_p2 < min_trials_per_profile
        continue;
    end

    y = double(T_cluster.spike_count);
    offset = log(time_bin_width) * ones(height(T_cluster), 1);

    speed_v = T_cluster.speed;
    tf_v = T_cluster.tf;
    sf_v = T_cluster.sf; sf_v(isnan(sf_v)) = 0;
    or_v = T_cluster.orientation; or_v(isnan(or_v)) = 0;
    t_since = T_cluster.time_since_onset;

    B_speed = glm_helpers.make_raised_cosine_basis(speed_v, n_speed_bases, speed_range(1), speed_range(2));
    B_tf    = glm_helpers.make_raised_cosine_basis(tf_v, n_tf_bases, tf_range(1), tf_range(2));
    B_onset = glm_helpers.make_onset_kernel_basis(t_since, n_onset_bases, onset_range(2));

    [X_null, ~]  = glm_helpers.assemble_design_matrix(B_speed, B_tf, B_onset, sf_v, or_v, 'Null');
    [X_speed, ~] = glm_helpers.assemble_design_matrix(B_speed, B_tf, B_onset, sf_v, or_v, 'M0_Speed');

    % Standard CV: trial-level, condition-stratified
    unique_trials = unique(T_cluster(:, {'trial_id', 'condition'}), 'rows');
    trial_fold = zeros(height(unique_trials), 1);
    conds = unique(unique_trials.condition);
    for cond_i = 1:length(conds)
        c_idx = find(unique_trials.condition == conds(cond_i));
        n_c = length(c_idx);
        perm = randperm(n_c);
        trial_fold(c_idx(perm)) = mod((1:n_c)' - 1, n_cv_folds) + 1;
    end
    [~, ~, row_to_trial] = unique(T_cluster(:, {'trial_id', 'condition'}), 'rows');
    fold_std = trial_fold(row_to_trial);

    % Cross-profile CV: strict profile split
    fold_cross = double(T_cluster.profile_id);

    [~, null_bps_std] = glm_helpers.cross_validate_glm(X_null, y, offset, fold_std, 0);
    [~, spd_bps_std]  = glm_helpers.cross_validate_glm(X_speed, y, offset, fold_std, 0);

    [~, null_bps_cross] = glm_helpers.cross_validate_glm(X_null, y, offset, fold_cross, 0);
    [~, spd_bps_cross]  = glm_helpers.cross_validate_glm(X_speed, y, offset, fold_cross, 0);

    delta_std = spd_bps_std - null_bps_std;
    delta_cross = spd_bps_cross - null_bps_cross;

    results.null_bps_std(ci) = null_bps_std;
    results.speed_bps_std(ci) = spd_bps_std;
    results.delta_bps_std(ci) = delta_std;
    results.null_bps_cross(ci) = null_bps_cross;
    results.speed_bps_cross(ci) = spd_bps_cross;
    results.delta_bps_cross(ci) = delta_cross;
    results.delta_drop(ci) = delta_std - delta_cross;
    results.usable(ci) = true;
end

valid_idx = results.usable & isfinite(results.delta_bps_std) & isfinite(results.delta_bps_cross);
n_valid = sum(valid_idx);
fprintf('Usable clusters: %d / %d\n', n_valid, n_clusters);

if n_valid == 0
    error('No clusters passed usability checks for cross-profile validation.');
end

standard_delta_bps = results.delta_bps_std(valid_idx);
cross_profile_delta_bps = results.delta_bps_cross(valid_idx);

%% ====================================================================
%  Section 6: Save outputs and figures
%  ====================================================================
fprintf('\n--- Saving results ---\n');

out_mat = fullfile(ctl.figs.curr_dir, 'glm_cross_profile_results.mat');
out_csv = fullfile(ctl.figs.curr_dir, 'glm_cross_profile_results.csv');
save(out_mat, 'results', 'standard_delta_bps', 'cross_profile_delta_bps', 'time_bin_width', ...
    'n_speed_bases', 'n_tf_bases', 'n_onset_bases', 'n_cv_folds', 'n_profiles', 'n_profile_samples');
writetable(results, out_csv);

fig = glm_plotting.plot_cross_profile_cv_comparison(standard_delta_bps, cross_profile_delta_bps);
if save_figs
    glm_plotting.save_figure(fig, fullfile(ctl.figs.curr_dir, 'cross_profile_generalization.png'));
end

fprintf('Saved:\n');
fprintf('  %s\n', out_mat);
fprintf('  %s\n', out_csv);
fprintf('  %s\n', fullfile(ctl.figs.curr_dir, 'cross_profile_generalization.png'));

fprintf('\nCross-profile validation complete.\n');
diary off;


function restricted = get_restricted_flag()
%GET_RESTRICTED_FLAG Runtime helper to keep restricted mode configurable.
    restricted = false;
end
