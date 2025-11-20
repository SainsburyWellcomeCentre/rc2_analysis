% Plot baseline vs. response firing rates for motion cloud experiments
% Similar to single_cluster_svm_unity_plot.m but with 2-column plots instead of unity plots
% and filtering for specific motion cloud names

experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'motionClouds', 'passive_same_luminance_mc', 'stationary_vs_motion_orientations', 'single_cluster'};


% Parameters: Motion cloud trial filters

name_prefix1 = 'theta-0p785';
name_prefix2 = 'theta0p000';
exclude_substrings_prefix2 = {'VX0p000'}; % e.g. exclude zero-velocity variants
name_prefix3 = 'theta0p785';
name_prefix4 = 'theta1p571';

% Portable substring matcher (older MATLAB lacks contains/startsWith)
contains_any = @(str, subs) any(cellfun(@(s) ...
    (exist('contains','builtin') && contains(str, s)) || (~exist('contains','builtin') && ~isempty(strfind(str, s))), subs));

% Load motion cloud sequence and folder names
mc_sequence = [];
cloud_names = {};
proto_seq_path = fullfile('D:\mvelez\mateoData_mc', 'motion_cloud_sequence_250414.mat');
if exist(proto_seq_path,'file')
    P = load(proto_seq_path);
    if isfield(P,'presentation_sequence')
        mc_sequence = P.presentation_sequence;
    else
        fns = fieldnames(P);
        for i3=1:numel(fns)
            v = P.(fns{i3});
            if isnumeric(v) && (isvector(v) || ismatrix(v))
                mc_sequence = v; break;
            end
        end
    end
end

% Load cloud names from image_folders.mat
folders_path = fullfile('D:\mvelez\mateoData_mc', 'image_folders.mat');
if exist(folders_path,'file')
    S = load(folders_path);
    fns = fieldnames(S);
    for i3=1:numel(fns)
        v = S.(fns{i3});
        if iscell(v)
            cloud_names = v; break;
        elseif isstring(v)
            cloud_names = cellstr(v(:)); break;
        elseif ischar(v)
            cloud_names = cellstr(v); break;
        elseif isstruct(v)
            if isfield(v, 'name')
                try
                    cloud_names = {v.name}; break;
                catch
                end
            end
        end
    end
end
if ~isempty(cloud_names), cloud_names = cloud_names(:)'; end

%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

x_all1                  = cell(1, length(probe_ids));
y_all1                  = cell(1, length(probe_ids));
p_val1                  = cell(1, length(probe_ids));
direction1              = cell(1, length(probe_ids));
x_all2                  = cell(1, length(probe_ids));
y_all2                  = cell(1, length(probe_ids));
p_val2                  = cell(1, length(probe_ids));
direction2              = cell(1, length(probe_ids));
x_all3                  = cell(1, length(probe_ids));
y_all3                  = cell(1, length(probe_ids));
p_val3                  = cell(1, length(probe_ids));
direction3              = cell(1, length(probe_ids));
x_all4                  = cell(1, length(probe_ids));
y_all4                  = cell(1, length(probe_ids));
p_val4                  = cell(1, length(probe_ids));
direction4              = cell(1, length(probe_ids));
cluster_ids             = cell(1, length(probe_ids));
region_str              = cell(1, length(probe_ids));
% Results and metadata for non-parametric pipeline
wilcoxon_p_bsl_vs_resp = cell(1, length(probe_ids));
kruskal_p_orientation  = cell(1, length(probe_ids));
% interaction is not computed in the non-parametric pipeline
n_paired_trials        = cell(1, length(probe_ids));
n_sig_orients_tukey    = cell(1, length(probe_ids));
% Tukey post-hoc pair strings per probe/cluster/trial_group
tukey_pairs            = cell(1, length(probe_ids));
% Store tukey pairwise p-values as a 6-element vector per probe/cluster/trial_group
% ordering: [1v2, 1v3, 1v4, 2v3, 2v4, 3v4]
tukey_pairwise_pvals   = cell(1, length(probe_ids));

for ii = 1 : length(probe_ids)
    
    data                = ctl.load_formatted_data(probe_ids{ii});
    clusters            = data.VISp_clusters;

    % get list of cluster IDs and the regions
    cluster_ids{ii}     = data.VISp_cluster_ids;
    region_str{ii}      = {clusters(:).region_str};

    for jj = 1 : length(trial_group_labels)
        
        % skip if the trial group label is not in the experiment
        if ~data.check_trial_group(trial_group_labels{jj})
            continue
        end
        
        % Precompute kept trial IDs per protocol using available API
        trials = data.get_trials_with_trial_group_label(trial_group_labels{jj});
        trial_ids_list = [];
        if ~isempty(trials)
            if iscell(trials)
                trial_ids_list = cellfun(@(t) t.trial_id, trials);
            else
                trial_ids_list = trials.trial_id;
            end
            trial_ids_list = trial_ids_list(:)';
        end

        for kk = 1 : length(cluster_ids{ii})
            
            % Get baseline/response with their trial_ids and align
            [baseline_fr, base_ids] = data.stationary_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            [response_fr, resp_ids] = data.motion_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            [common_ids, ia, ib] = intersect(base_ids(:), resp_ids(:));
            baseline_fr = baseline_fr(ia);
            response_fr = response_fr(ib);

            % Build kept masks for both prefixes
            keep1 = false(size(common_ids));
            keep2 = false(size(common_ids));
            keep3 = false(size(common_ids));
            keep4 = false(size(common_ids));
            for ti = 1:numel(common_ids)
                tid = common_ids(ti);
                if ~isempty(mc_sequence) && ~isempty(cloud_names) && tid>=1 && tid<=length(mc_sequence)
                    mc_id = mc_sequence(tid);
                    if mc_id>=1 && mc_id<=length(cloud_names)
                        cname = cloud_names{mc_id};
                        if (exist('startsWith','builtin') && startsWith(cname, name_prefix1)) || ...
                           (~exist('startsWith','builtin') && strncmp(cname, name_prefix1, length(name_prefix1)))
                            keep1(ti) = true;
                        end
                        if ( (exist('startsWith','builtin') && startsWith(cname, name_prefix2)) || ...
                             (~exist('startsWith','builtin') && strncmp(cname, name_prefix2, length(name_prefix2))) )
                            if ~contains_any(cname, exclude_substrings_prefix2)
                                keep2(ti) = true;
                            end
                        end
                        if (exist('startsWith','builtin') && startsWith(cname, name_prefix3)) || ...
                           (~exist('startsWith','builtin') && strncmp(cname, name_prefix3, length(name_prefix3)))
                            keep3(ti) = true;
                        end
                        if (exist('startsWith','builtin') && startsWith(cname, name_prefix4)) || ...
                           (~exist('startsWith','builtin') && strncmp(cname, name_prefix4, length(name_prefix4)))
                            keep4(ti) = true;
                        end
                    end
                end
            end

            x_all1{ii}{kk}{jj} = baseline_fr(keep1);
            y_all1{ii}{kk}{jj} = response_fr(keep1);
            x_all2{ii}{kk}{jj} = baseline_fr(keep2);
            y_all2{ii}{kk}{jj} = response_fr(keep2);
            x_all3{ii}{kk}{jj} = baseline_fr(keep3);
            y_all3{ii}{kk}{jj} = response_fr(keep3);
            x_all4{ii}{kk}{jj} = baseline_fr(keep4);
            y_all4{ii}{kk}{jj} = response_fr(keep4);
            
            % -- New statistics: repeated-measures (within) Condition x Orientation (between)
            % Build orientation labels for the filtered trials
            orientation_idx = nan(size(common_ids));
            % Assign orientation indices in the order of prefixes (1..4). If a trial
            % matches multiple prefixes, the first match (1..4) is used.
            orientation_idx(keep1) = 1;
            orientation_idx(isnan(orientation_idx) & keep2) = 2;
            orientation_idx(isnan(orientation_idx) & keep3) = 3;
            orientation_idx(isnan(orientation_idx) & keep4) = 4;

            sel_mask = ~isnan(orientation_idx);
            baseline_sel = baseline_fr(sel_mask);
            response_sel = response_fr(sel_mask);
            orient_sel = orientation_idx(sel_mask);

            % initialize ANOVA outputs for this probe/cluster/trial_group
            pCondition = NaN; pOrientation = NaN; pInteraction = NaN; n_total = numel(baseline_sel);

            if n_total >= 2
                % (no diagnostic prints)

                % Non-parametric pipeline per request:
                % 1) Paired Wilcoxon signed-rank across all Baseline vs Response (responsiveness)
                % 2) Kruskal-Wallis on Response only to test Orientation effect
                % 3) Tukey-Kramer post-hoc (multcompare) to find which orientations differ
                try
                    % Paired Wilcoxon signed-rank (Condition main effect)
                    if exist('signrank','file') == 2
                        try
                            pCondition = signrank(baseline_sel, response_sel);
                        catch
                            pCondition = NaN;
                        end
                    else
                        pCondition = NaN;
                    end

                    % Kruskal-Wallis for Orientation effect (responses only)
                    pOrientation = NaN;
                    pInteraction = NaN; % interaction not tested in this non-parametric pipeline
                    tukey_n_sig_orients = 0;
                    try
                        uniq_orients = unique(orient_sel);
                        if numel(uniq_orients) >= 2
                            % suppress figure with 'off'
                            [pOrientation, ~, stats] = kruskalwallis(response_sel, orient_sel, 'off');
                            % Tukey-Kramer post-hoc (multcompare) to identify differing orientations
                            try
                                c = multcompare(stats, 'CType', 'tukey-kramer', 'Display', 'off');
                                % c columns: [group1, group2, lower, diff, upper, pvalue]
                                % Populate a 6-element vector of pairwise p-values for canonical pairs
                                % ordering: [1v2, 1v3, 1v4, 2v3, 2v4, 3v4]
                                pvals6 = nan(1,6);
                                if ~isempty(c)
                                    pvals_pairs = c(:,6);
                                    % Map each returned comparison to the canonical index
                                    for r = 1:size(c,1)
                                        g1 = c(r,1); g2 = c(r,2);
                                        if g1 > g2, tmp = g1; g1 = g2; g2 = tmp; end
                                        idx = [];
                                        if g1==1 && g2==2, idx = 1; end
                                        if g1==1 && g2==3, idx = 2; end
                                        if g1==1 && g2==4, idx = 3; end
                                        if g1==2 && g2==3, idx = 4; end
                                        if g1==2 && g2==4, idx = 5; end
                                        if g1==3 && g2==4, idx = 6; end
                                        if ~isempty(idx)
                                            pvals6(idx) = c(r,6);
                                        end
                                    end
                                    sig_pairs = pvals_pairs < 0.05;
                                    if any(sig_pairs)
                                        sig_groups = unique([c(sig_pairs,1); c(sig_pairs,2)]);
                                        tukey_n_sig_orients = numel(sig_groups);
                                        % build pair strings like '1v2;1v3'
                                        pair_idx = find(sig_pairs);
                                        pair_strs = arrayfun(@(r) sprintf('%dv%d', c(r,1), c(r,2)), pair_idx, 'UniformOutput', false);
                                        tukey_sig_pairs = strjoin(pair_strs, ';');
                                        % store in top-level container for strings
                                        if numel(tukey_pairs) < ii || isempty(tukey_pairs{ii})
                                            tukey_pairs{ii} = cell(1, length(cluster_ids{ii}));
                                        end
                                        if numel(tukey_pairs{ii}) < kk || isempty(tukey_pairs{ii}{kk})
                                            tukey_pairs{ii}{kk} = cell(1, numel(trial_group_labels));
                                        end
                                        tukey_pairs{ii}{kk}{jj} = tukey_sig_pairs;
                                    else
                                        tukey_n_sig_orients = 0;
                                        tukey_sig_pairs = '';
                                    end
                                else
                                    tukey_n_sig_orients = 0;
                                    tukey_sig_pairs = '';
                                end
                                % store numeric pairwise p-values in top-level container
                                if numel(tukey_pairwise_pvals) < ii || isempty(tukey_pairwise_pvals{ii})
                                    tukey_pairwise_pvals{ii} = cell(1, length(cluster_ids{ii}));
                                end
                                if numel(tukey_pairwise_pvals{ii}) < kk || isempty(tukey_pairwise_pvals{ii}{kk})
                                    tukey_pairwise_pvals{ii}{kk} = cell(1, numel(trial_group_labels));
                                end
                                tukey_pairwise_pvals{ii}{kk}{jj} = pvals6;
                            catch
                                tukey_n_sig_orients = 0;
                                tukey_sig_pairs = '';
                            end
                        else
                            pOrientation = NaN;
                            tukey_n_sig_orients = 0;
                        end
                    catch
                        pOrientation = NaN;
                        tukey_n_sig_orients = 0;
                    end
                catch
                    % leave NaNs on unexpected errors
                    pCondition = NaN; pOrientation = NaN; pInteraction = NaN; tukey_n_sig_orients = 0;
                end
            end

            % Store summary metadata for non-parametric pipeline
            wilcoxon_p_bsl_vs_resp{ii}{kk}(jj) = pCondition;
            kruskal_p_orientation{ii}{kk}(jj) = pOrientation;
            n_paired_trials{ii}{kk}(jj) = n_total;

            % Now perform per-orientation paired tests and Bonferroni-correct
            p_unc = NaN(1,4);
            dir_vals = zeros(1,4);
            for pidx = 1:4
                xb = [];
                yr = [];
                switch pidx
                    case 1
                        xb = x_all1{ii}{kk}{jj}(:);
                        yr = y_all1{ii}{kk}{jj}(:);
                    case 2
                        xb = x_all2{ii}{kk}{jj}(:);
                        yr = y_all2{ii}{kk}{jj}(:);
                    case 3
                        xb = x_all3{ii}{kk}{jj}(:);
                        yr = y_all3{ii}{kk}{jj}(:);
                    case 4
                        xb = x_all4{ii}{kk}{jj}(:);
                        yr = y_all4{ii}{kk}{jj}(:);
                end
                good = isfinite(xb) & isfinite(yr);
                xb = xb(good); yr = yr(good);
                if numel(xb) >= 2
                    % prefer signrank when available
                    if exist('signrank','file') == 2
                        try
                            p_unc(pidx) = signrank(xb, yr);
                        catch
                            try
                                [~, p_unc(pidx)] = ttest(xb, yr);
                            catch
                                p_unc(pidx) = NaN;
                            end
                        end
                    else
                        try
                            [~, p_unc(pidx)] = ttest(xb, yr);
                        catch
                            p_unc(pidx) = NaN;
                        end
                    end
                else
                    p_unc(pidx) = NaN;
                end
            end

            % Bonferroni correction across the 4 orientations
            p_corr = min(1, p_unc * 4);

            % Determine directions only for orientations that survive Bonferroni
            sig_mask = p_corr < 0.05;
            n_sig_orients_tukey{ii}{kk}(jj) = sum(sig_mask);
            for pidx = 1:4
                xb = [];
                yr = [];
                switch pidx
                    case 1
                        xb = x_all1{ii}{kk}{jj}(:);
                        yr = y_all1{ii}{kk}{jj}(:);
                    case 2
                        xb = x_all2{ii}{kk}{jj}(:);
                        yr = y_all2{ii}{kk}{jj}(:);
                    case 3
                        xb = x_all3{ii}{kk}{jj}(:);
                        yr = y_all3{ii}{kk}{jj}(:);
                    case 4
                        xb = x_all4{ii}{kk}{jj}(:);
                        yr = y_all4{ii}{kk}{jj}(:);
                end
                good = isfinite(xb) & isfinite(yr);
                xb = xb(good); yr = yr(good);
                pv = p_corr(pidx);
                % direction only meaningful if orientation survived Bonferroni
                dirn = 0;
                if ~isnan(pv) && pv < 0.05 && numel(xb) >= 1
                    d = median(yr - xb);
                    if d > 0
                        dirn = 1;
                    elseif d < 0
                        dirn = -1;
                    else
                        dirn = 0;
                    end
                end
                % Store corrected p-value and direction per orientation
                switch pidx
                    case 1
                        p_val1{ii}{kk}(jj) = pv;
                        direction1{ii}{kk}(jj) = dirn;
                    case 2
                        p_val2{ii}{kk}(jj) = pv;
                        direction2{ii}{kk}(jj) = dirn;
                    case 3
                        p_val3{ii}{kk}(jj) = pv;
                        direction3{ii}{kk}(jj) = dirn;
                    case 4
                        p_val4{ii}{kk}(jj) = pv;
                        direction4{ii}{kk}(jj) = dirn;
                end
            end
            % Prefer Tukey-Kramer post-hoc count (from kruskalwallis + multcompare)
            if exist('tukey_n_sig_orients','var') && ~isempty(tukey_n_sig_orients)
                n_sig_orients_tukey{ii}{kk}(jj) = tukey_n_sig_orients;
            end
        end
    end
end

% Print summary
for ii = 1:length(probe_ids)
    data = ctl.load_formatted_data(probe_ids{ii});
    kept_counts1 = zeros(1, numel(trial_group_labels));
    kept_counts2 = zeros(1, numel(trial_group_labels));
    kept_counts3 = zeros(1, numel(trial_group_labels));
    kept_counts4 = zeros(1, numel(trial_group_labels));
    for jj = 1:numel(trial_group_labels)
        if ~data.check_trial_group(trial_group_labels{jj})
            kept_counts1(jj) = 0; kept_counts2(jj) = 0; continue
        end
        trials = data.get_trials_with_trial_group_label(trial_group_labels{jj});
        if isempty(trials)
            kept_counts1(jj) = 0; kept_counts2(jj) = 0; kept_counts3(jj) = 0; kept_counts4(jj) = 0; continue
        end
        if iscell(trials)
            trial_ids_list = cellfun(@(t) t.trial_id, trials);
        else
            trial_ids_list = trials.trial_id;
        end
        keep_mask1 = false(size(trial_ids_list));
        keep_mask2 = false(size(trial_ids_list));
        keep_mask3 = false(size(trial_ids_list));
        keep_mask4 = false(size(trial_ids_list));
        for ti = 1:numel(trial_ids_list)
            tid = trial_ids_list(ti);
            if ~isempty(mc_sequence) && ~isempty(cloud_names) && tid>=1 && tid<=length(mc_sequence)
                mc_id = mc_sequence(tid);
                if mc_id>=1 && mc_id<=length(cloud_names)
                    cname = cloud_names{mc_id};
                    if (exist('startsWith','builtin') && startsWith(cname, name_prefix1)) || ...
                       (~exist('startsWith','builtin') && strncmp(cname, name_prefix1, length(name_prefix1)))
                        keep_mask1(ti) = true;
                    end
                    if ( (exist('startsWith','builtin') && startsWith(cname, name_prefix2)) || ...
                         (~exist('startsWith','builtin') && strncmp(cname, name_prefix2, length(name_prefix2))) )
                        if ~contains_any(cname, exclude_substrings_prefix2)
                            keep_mask2(ti) = true;
                        end
                    end
                    if (exist('startsWith','builtin') && startsWith(cname, name_prefix3)) || ...
                       (~exist('startsWith','builtin') && strncmp(cname, name_prefix3, length(name_prefix3)))
                        keep_mask3(ti) = true;
                    end
                    if (exist('startsWith','builtin') && startsWith(cname, name_prefix4)) || ...
                       (~exist('startsWith','builtin') && strncmp(cname, name_prefix4, length(name_prefix4)))
                        keep_mask4(ti) = true;
                    end
                end
            end
        end
        kept_counts1(jj) = sum(keep_mask1);
        kept_counts2(jj) = sum(keep_mask2);
        kept_counts3(jj) = sum(keep_mask3);
        kept_counts4(jj) = sum(keep_mask4);
    end
    fprintf('Probe %s: kept trials p1=%s [%d %d %d]; p2=%s [%d %d %d]; p3=%s [%d %d %d]; p4=%s [%d %d %d]\n', ...
        probe_ids{ii}, name_prefix1, kept_counts1(1), kept_counts1(2), kept_counts1(3), ...
        name_prefix2, kept_counts2(1), kept_counts2(2), kept_counts2(3), ...
        name_prefix3, kept_counts3(1), kept_counts3(2), kept_counts3(3), ...
        name_prefix4, kept_counts4(1), kept_counts4(2), kept_counts4(3));
end

%%
ctl.setup_figures(figure_dir, save_figs);

% Prepare CSV summary data
probe_id_list = {};
cluster_id_list = [];
trial_group_list = {};
prefix_list = {};
trial_number_list = [];
baseline_fr_list = [];
response_fr_list = [];
difference_list = [];
paired_p_bonferroni_list = [];
direction_list = [];
% Summary columns for the non-parametric pipeline
wilcoxon_p_bsl_vs_resp_list = [];
kruskal_p_orientation_list = [];
% Tukey pairwise p-values per CSV row (canonical six pairs: 1v2,1v3,1v4,2v3,2v4,3v4)
tukey_p_1v2_list = [];
tukey_p_1v3_list = [];
tukey_p_1v4_list = [];
tukey_p_2v3_list = [];
tukey_p_2v4_list = [];
tukey_p_3v4_list = [];

for ii = 1 : length(probe_ids)
    for kk = 1 : length(cluster_ids{ii})
        for jj = 1 : length(trial_group_labels)
            % Get data for this combination (match the order used during ANOVA computation)
            xb1 = x_all1{ii}{kk}{jj};
            yr1 = y_all1{ii}{kk}{jj};
            xb2 = x_all2{ii}{kk}{jj};
            yr2 = y_all2{ii}{kk}{jj};
            xb3 = x_all3{ii}{kk}{jj};
            yr3 = y_all3{ii}{kk}{jj};
            xb4 = x_all4{ii}{kk}{jj};
            yr4 = y_all4{ii}{kk}{jj};

            % Add data for each prefix
            prefixes = {name_prefix1, name_prefix2, name_prefix3, name_prefix4};
            baseline_data = {xb1, xb2, xb3, xb4};
            response_data = {yr1, yr2, yr3, yr4};

            for p = 1:4
                baseline = baseline_data{p};
                response = response_data{p};

                if ~isempty(baseline) && ~isempty(response)
                    n_trials = min(length(baseline), length(response));
                    for t = 1:n_trials
                        probe_id_list{end+1} = probe_ids{ii};
                        cluster_id_list(end+1) = cluster_ids{ii}(kk);
                        trial_group_list{end+1} = trial_group_labels{jj};
                        prefix_list{end+1} = prefixes{p};
                        trial_number_list(end+1) = t;
                        baseline_fr_list(end+1) = baseline(t);
                        response_fr_list(end+1) = response(t);
                        difference_list(end+1) = response(t) - baseline(t);
                        % pick the correct p-value and direction for this prefix
                        switch p
                            case 1
                                pv = p_val1{ii}{kk}(jj);
                                dn = direction1{ii}{kk}(jj);
                            case 2
                                pv = p_val2{ii}{kk}(jj);
                                dn = direction2{ii}{kk}(jj);
                            case 3
                                pv = p_val3{ii}{kk}(jj);
                                dn = direction3{ii}{kk}(jj);
                            case 4
                                pv = p_val4{ii}{kk}(jj);
                                dn = direction4{ii}{kk}(jj);
                        end
                        paired_p_bonferroni_list(end+1) = pv;
                        direction_list(end+1) = dn;
                        % Summary-level values (same for all prefixes within this cluster/trial-group)
                        wilcoxon_p_bsl_vs_resp_list(end+1) = wilcoxon_p_bsl_vs_resp{ii}{kk}(jj);
                        kruskal_p_orientation_list(end+1) = kruskal_p_orientation{ii}{kk}(jj);
                        % Tukey pairwise p-values (six canonical pairs)
                        pvals6 = nan(1,6);
                        if exist('tukey_pairwise_pvals','var') && numel(tukey_pairwise_pvals) >= ii && ~isempty(tukey_pairwise_pvals{ii}) && numel(tukey_pairwise_pvals{ii}) >= kk && ~isempty(tukey_pairwise_pvals{ii}{kk}) && numel(tukey_pairwise_pvals{ii}{kk}) >= jj
                            pvals6 = tukey_pairwise_pvals{ii}{kk}{jj};
                        end
                        tukey_p_1v2_list(end+1) = pvals6(1);
                        tukey_p_1v3_list(end+1) = pvals6(2);
                        tukey_p_1v4_list(end+1) = pvals6(3);
                        tukey_p_2v3_list(end+1) = pvals6(4);
                        tukey_p_2v4_list(end+1) = pvals6(5);
                        tukey_p_3v4_list(end+1) = pvals6(6);
                    end
                end
            end
        end
    end
end

% Convert to table and save CSV
    if ~isempty(probe_id_list)
    csv_table = table(probe_id_list', cluster_id_list', trial_group_list', prefix_list', ...
                     trial_number_list', baseline_fr_list', response_fr_list', difference_list', ...
                     wilcoxon_p_bsl_vs_resp_list', kruskal_p_orientation_list', paired_p_bonferroni_list', direction_list', ...
                     tukey_p_1v2_list', tukey_p_1v3_list', tukey_p_1v4_list', tukey_p_2v3_list', tukey_p_2v4_list', tukey_p_3v4_list', ...
                     'VariableNames', {'probe_id','cluster_id','trial_group','prefix','trial_number','baseline_fr','response_fr','difference','wilcoxon_p_bsl_vs_resp','kruskal_p_orientation','paired_p_bonferroni','response_direction','tukey_p_1v2','tukey_p_1v3','tukey_p_1v4','tukey_p_2v3','tukey_p_2v4','tukey_p_3v4'});
    
    csv_filename = fullfile(ctl.figs.curr_dir, 'motion_cloud_orientations_analysis_summary.csv');
    writetable(csv_table, csv_filename);
    fprintf('Saved analysis summary to: %s\n', csv_filename);
end

for ii = 1 : length(probe_ids)
    for jj = 1 : length(cluster_ids{ii})
        fprintf('Probe %s: plotting cluster_id=%d\n', probe_ids{ii}, cluster_ids{ii}(jj));
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        axs                     = gobjects(0);
        panel_mins              = [];
        panel_maxs              = [];

        for kk = 1 : length(trial_group_labels)
            
            % Left column: baseline vs response plots
            pos         = plot_array.get_position(kk);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            axs(end+1)  = h_ax;
            
            % Get data for this trial group (prefix1)
            xb1 = x_all1{ii}{jj}{kk};
            yr1 = y_all1{ii}{jj}{kk};
            % and for prefix2
            xb2 = x_all2{ii}{jj}{kk};
            yr2 = y_all2{ii}{jj}{kk};
            % and for prefix3/4
            xb3p = x_all3{ii}{jj}{kk};
            yr3p = y_all3{ii}{jj}{kk};
            xb4p = x_all4{ii}{jj}{kk};
            yr4p = y_all4{ii}{jj}{kk};
            
            if ~isempty(xb1) && ~isempty(yr1)
                xb = xb1(:); yr = yr1(:);
                scatter(h_ax, ones(numel(xb),1), xb, 20, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax, 'on');
                scatter(h_ax, 2*ones(numel(yr),1), yr, 20, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
                n_pairs = min(numel(xb), numel(yr));
                for i = 1:n_pairs
                    plot(h_ax, [1 2], [xb(i) yr(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.3]);
                end
                % Medians
                med_b = median(xb);
                med_r = median(yr);
                scatter(h_ax, 1, med_b, 60, [0 0 0], 'filled');
                scatter(h_ax, 2, med_r, 60, [0 0 0], 'filled');
                if p_val1{ii}{jj}(kk) < 0.05
                    if direction1{ii}{jj}(kk) > 0
                        scatter(h_ax, 2*ones(numel(yr),1), yr, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
                    else
                        scatter(h_ax, 2*ones(numel(yr),1), yr, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
                    end
                else
                    scatter(h_ax, 2*ones(numel(yr),1), yr, 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.8);
                end
            end

            if ~isempty(xb2) && ~isempty(yr2)
                xb = xb2(:); yr = yr2(:);
                scatter(h_ax, 3*ones(numel(xb),1), xb, 20, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax, 'on');
                scatter(h_ax, 4*ones(numel(yr),1), yr, 20, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
            n_pairs = min(numel(xb), numel(yr));
                for i = 1:n_pairs
                    plot(h_ax, [3 4], [xb(i) yr(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.3]);
                end
                % Medians
                med_b = median(xb);
                med_r = median(yr);
                scatter(h_ax, 3, med_b, 60, [0 0 0], 'filled');
                scatter(h_ax, 4, med_r, 60, [0 0 0], 'filled');
                if p_val2{ii}{jj}(kk) < 0.05
                    if direction2{ii}{jj}(kk) > 0
                        scatter(h_ax, 4*ones(numel(yr),1), yr, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
                    else
                        scatter(h_ax, 4*ones(numel(yr),1), yr, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
                    end
                else
                    scatter(h_ax, 4*ones(numel(yr),1), yr, 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.8);
                end
            end

            if ~isempty(xb3p) && ~isempty(yr3p)
                xb = xb3p(:); yr = yr3p(:);
                scatter(h_ax, 5*ones(numel(xb),1), xb, 20, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax, 'on');
                scatter(h_ax, 6*ones(numel(yr),1), yr, 20, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
                n_pairs = min(numel(xb), numel(yr));
                for i = 1:n_pairs
                    plot(h_ax, [5 6], [xb(i) yr(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.3]);
                end
                % Medians
                med_b = median(xb);
                med_r = median(yr);
                scatter(h_ax, 5, med_b, 60, [0 0 0], 'filled');
                scatter(h_ax, 6, med_r, 60, [0 0 0], 'filled');
                if p_val3{ii}{jj}(kk) < 0.05
                    if direction3{ii}{jj}(kk) > 0
                        scatter(h_ax, 6*ones(numel(yr),1), yr, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
                    else
                        scatter(h_ax, 6*ones(numel(yr),1), yr, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
                    end
                else
                    scatter(h_ax, 6*ones(numel(yr),1), yr, 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.8);
                end
            end

            if ~isempty(xb4p) && ~isempty(yr4p)
                xb = xb4p(:); yr = yr4p(:);
                scatter(h_ax, 7*ones(numel(xb),1), xb, 20, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax, 'on');
                scatter(h_ax, 8*ones(numel(yr),1), yr, 20, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
                n_pairs = min(numel(xb), numel(yr));
                for i = 1:n_pairs
                    plot(h_ax, [7 8], [xb(i) yr(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.3]);
                end
                % Medians
                med_b = median(xb);
                med_r = median(yr);
                scatter(h_ax, 7, med_b, 60, [0 0 0], 'filled');
                scatter(h_ax, 8, med_r, 60, [0 0 0], 'filled');
                if p_val4{ii}{jj}(kk) < 0.05
                    if direction4{ii}{jj}(kk) > 0
                        scatter(h_ax, 8*ones(numel(yr),1), yr, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
                    else
                        scatter(h_ax, 8*ones(numel(yr),1), yr, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
                    end
                else
                    scatter(h_ax, 8*ones(numel(yr),1), yr, 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.8);
                end
            end

            xlim(h_ax, [0.5 8.5]);
            xlabel(h_ax, 'Condition');
            ylabel(h_ax, 'Firing Rate (Hz)');
            title(h_ax, trial_group_labels{kk});
            
            % Set x-axis labels
            set(h_ax, 'XTick', [1 2 3 4 5 6 7 8], 'XTickLabel', {'B1','R1','B2','R2','B3','R3','B4','R4'});
            
            % Store min/max for axis synchronization
            local_vals = [];
            if exist('xb1','var') && exist('yr1','var') && ~isempty(xb1) && ~isempty(yr1)
                local_vals = [local_vals; xb1(:); yr1(:)];
            end
            if exist('xb2','var') && exist('yr2','var') && ~isempty(xb2) && ~isempty(yr2)
                local_vals = [local_vals; xb2(:); yr2(:)];
            end
            if exist('xb3p','var') && exist('yr3p','var') && ~isempty(xb3p) && ~isempty(yr3p)
                local_vals = [local_vals; xb3p(:); yr3p(:)];
            end
            if exist('xb4p','var') && exist('yr4p','var') && ~isempty(xb4p) && ~isempty(yr4p)
                local_vals = [local_vals; xb4p(:); yr4p(:)];
            end
            if ~isempty(local_vals)
                panel_mins(end+1) = min(local_vals);
                panel_maxs(end+1) = max(local_vals);
            else
                panel_mins(end+1) = 0;
                panel_maxs(end+1) = 1;
            end
            
            % Right column: response-baseline difference plots
            pos2        = plot_array.get_position(kk + 3); % Right column
            h_ax2       = axes('units', 'centimeters', 'position', pos2);
            axs(end+1)  = h_ax2;
            
            % Plot differences for each prefix
            if ~isempty(xb1) && ~isempty(yr1)
                diff1 = yr1(:) - xb1(:);
                scatter(h_ax2, 1*ones(numel(diff1),1), diff1, 50, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax2, 'on');
                % Median difference
                med_diff1 = median(diff1);
                scatter(h_ax2, 1, med_diff1, 70, [0 0 0], 'filled');
            end
            
            if ~isempty(xb2) && ~isempty(yr2)
                diff2 = yr2(:) - xb2(:);
                scatter(h_ax2, 2*ones(numel(diff2),1), diff2, 50, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax2, 'on');
                % Median difference
                med_diff2 = median(diff2);
                scatter(h_ax2, 2, med_diff2, 70, [0 0 0], 'filled');
            end
            
            if ~isempty(xb3p) && ~isempty(yr3p)
                diff3 = yr3p(:) - xb3p(:);
                scatter(h_ax2, 3*ones(numel(diff3),1), diff3, 50, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax2, 'on');
                % Median difference
                med_diff3 = median(diff3);
                scatter(h_ax2, 3, med_diff3, 70, [0 0 0], 'filled');
            end
            
            if ~isempty(xb4p) && ~isempty(yr4p)
                diff4 = yr4p(:) - xb4p(:);
                scatter(h_ax2, 4*ones(numel(diff4),1), diff4, 50, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
                hold(h_ax2, 'on');
                % Median difference
                med_diff4 = median(diff4);
                scatter(h_ax2, 4, med_diff4, 70, [0 0 0], 'filled');
            end
            
            xlim(h_ax2, [0.5 4.5]);
            xlabel(h_ax2, 'Motion Cloud Group');
            ylabel(h_ax2, 'Response - Baseline (Hz)');
            title(h_ax2, [trial_group_labels{kk} ' (Differences)']);
            
            % Set x-axis labels for differences
            set(h_ax2, 'XTick', [1 2 3 4], 'XTickLabel', {'P1','P2','P3','P4'});
            
            % Store min/max for difference plot
            diff_vals = [];
            if exist('diff1','var') && ~isempty(diff1)
                diff_vals = [diff_vals; diff1(:)];
            end
            if exist('diff2','var') && ~isempty(diff2)
                diff_vals = [diff_vals; diff2(:)];
            end
            if exist('diff3','var') && ~isempty(diff3)
                diff_vals = [diff_vals; diff3(:)];
            end
            if exist('diff4','var') && ~isempty(diff4)
                diff_vals = [diff_vals; diff4(:)];
            end
            if ~isempty(diff_vals)
                panel_mins(end+1) = min(diff_vals);
                panel_maxs(end+1) = max(diff_vals);
            else
                panel_mins(end+1) = 0;
                panel_maxs(end+1) = 1;
            end
        end

        % Sync y-axes across all panels
        if ~isempty(panel_mins)
            y_min = min(panel_mins);
            y_max = max(panel_maxs);
            for kk = 1:length(axs)
                ylim(axs(kk), [y_min, y_max]);
            end
        end

         % give the page a title
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', probe_ids{ii}, cluster_ids{ii}(jj), region_str{ii}{jj}));
        
        ctl.figs.save_fig_to_join(); %
    end

    ctl.figs.join_figs(sprintf('%s.pdf', probe_ids{ii}), overwrite);
    ctl.figs.clear_figs();
end