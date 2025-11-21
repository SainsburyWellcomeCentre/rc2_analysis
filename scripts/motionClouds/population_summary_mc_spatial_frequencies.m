% population_summary_mc_spatial_frequencies.m
% Summarise population classification for spatial frequencies (V and VT)
% Produces percentage bar plots and CSV summaries in the summary folder.

input_csv = fullfile('C:','Users','lee','Documents','mvelez','figures','motionClouds', ...
    'passive_same_luminance_mc','stationary_vs_motion_spatial_frequencies','single_cluster', ...
    'motion_cloud_spatial_frequencies_analysis_summary.csv');
if ~exist(input_csv,'file')
    error('Input CSV not found: %s', input_csv);
end

out_dir = fullfile('C:','Users','lee','Documents','mvelez','figures','motionClouds', ...
    'passive_same_luminance_mc','stationary_vs_motion_spatial_frequencies','summary');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

T = readtable(input_csv);
% normalize prefix list
if iscell(T.prefix), prefixes_all = unique(T.prefix,'stable'); else prefixes_all = unique(cellstr(T.prefix),'stable'); end
prefixes_all = cellstr(prefixes_all(:)');

% unique clusters
[cluster_keys, ~] = unique(T(:,{'probe_id','cluster_id'}),'rows','stable');
cluster_keys = cluster_keys(:,:);

nClusters = height(cluster_keys);

% build clusters struct like other scripts
trial_groups = {'V','VT'};
clusters = struct();
for i = 1:nClusters
    pid = cluster_keys.probe_id{i};
    cid = cluster_keys.cluster_id(i);
    clusters(i).probe_id = pid;
    clusters(i).cluster_id = cid;
    for g = 1:length(trial_groups)
        grp = trial_groups{g};
        clusters(i).(grp).medians = nan(1,numel(prefixes_all));
        clusters(i).(grp).direction = nan(1,numel(prefixes_all));
        clusters(i).(grp).wilcoxon = NaN;
        clusters(i).(grp).kruskal = NaN;
        for p = 1:numel(prefixes_all)
            pref = prefixes_all{p};
            sel = strcmp(T.probe_id, pid) & (T.cluster_id == cid) & strcmp(T.trial_group, grp) & strcmp(T.prefix, pref);
            if any(sel)
                clusters(i).(grp).medians(p) = median(T.difference(sel),'omitnan');
                rd = unique(T.response_direction(sel)); rd = rd(~isnan(rd));
                if isempty(rd), clusters(i).(grp).direction(p) = 0; else clusters(i).(grp).direction(p) = rd(1); end
            else
                clusters(i).(grp).medians(p) = NaN;
                clusters(i).(grp).direction(p) = NaN;
            end
        end
        sel_all = strcmp(T.probe_id, pid) & (T.cluster_id == cid) & strcmp(T.trial_group, grp);
        if any(sel_all)
            if ismember('wilcoxon_p_bsl_vs_resp', T.Properties.VariableNames)
                wvals = unique(T.wilcoxon_p_bsl_vs_resp(sel_all)); wvals = wvals(~isnan(wvals));
                if ~isempty(wvals), clusters(i).(grp).wilcoxon = wvals(1); end
            end
            if ismember('kruskal_p_orientation', T.Properties.VariableNames)
                kvals = unique(T.kruskal_p_orientation(sel_all)); kvals = kvals(~isnan(kvals));
                if ~isempty(kvals), clusters(i).(grp).kruskal = kvals(1); end
            end
        end
    end
end

% Helper: tuned bound = num prefixes - 1
tuned_upper = numel(prefixes_all) - 1;

% run for both groups
groups = {'V','VT'};
for gi = 1:2
    grp = groups{gi};
    [counts, tuned_idx, pref_counts] = classify_group(clusters, nClusters, prefixes_all, grp, tuned_upper);
    % build category table and percentages
    cats = {'tuned','responsive_not_tuned','not_responsive','other'}';
    cnts = [counts.tuned; counts.resp_not_tuned; counts.not_responsive; counts.other];
    pct = 100 * cnts / nClusters;
    Tcats = table(cats, cnts, pct, 'VariableNames', {'category','count','pct'});
    out_csv = fullfile(out_dir, sprintf('population_spatial_summary_%s_categories.csv', grp));
    writetable(Tcats, out_csv);
    % plot categories as percentage bar
    fig = figure('Visible','off');
    bar(categorical(cats), pct);
    ylabel('Percentage of clusters');
    title(sprintf('Population categories (%s) - spatial frequencies', grp));
    ylim([0 100]);
    pdfname = fullfile(out_dir, sprintf('population_spatial_summary_%s_categories.pdf', grp));
    set(fig,'PaperPositionMode','auto'); print(fig, pdfname, '-dpdf'); close(fig);

    % preferred distribution for tuned clusters
    pref_cnts = pref_counts(:);
    if sum(pref_cnts) > 0
        pref_pct = 100 * pref_cnts / sum(pref_cnts);
        Tpref = table(prefixes_all(:), pref_cnts, pref_pct, 'VariableNames', {'prefix','count','pct_of_tuned'});
        out_csv2 = fullfile(out_dir, sprintf('population_spatial_summary_%s_pref_distribution.csv', grp));
        writetable(Tpref, out_csv2);
        fig2 = figure('Visible','off');
        bar(categorical(prefixes_all), pref_pct);
        ylabel('Percent of tuned clusters');
        xlabel('Prefix');
        title(sprintf('Preferred prefix distribution (%s) - tuned clusters (spatial)', grp));
        ylim([0 100]);
        pdfname2 = fullfile(out_dir, sprintf('population_spatial_summary_%s_pref_distribution.pdf', grp));
        set(fig2,'PaperPositionMode','auto'); print(fig2, pdfname2, '-dpdf'); close(fig2);
    else
        fprintf('No tuned clusters for %s, skipping preferred distribution output.\n', grp);
    end
    fprintf('Saved summaries for %s -> %s\n', grp, out_dir);
end

fprintf('Population spatial summaries complete. Files written to: %s\n', out_dir);

% local function definitions (must appear at end of script)
function [counts, tuned_idx, pref_counts] = classify_group(clusters, nClusters, prefixes_all, grp, tuned_upper)
    is_tuned = false(1,nClusters);
    is_resp_not_tuned = false(1,nClusters);
    is_not_responsive = false(1,nClusters);
    tuned_pref_idx = nan(1,nClusters);
    for i=1:nClusters
        d = clusters(i).(grp).direction;
        w = clusters(i).(grp).wilcoxon;
        k = clusters(i).(grp).kruskal;
        med_overall = median(clusters(i).(grp).medians,'omitnan');
        nOnes = sum(~isnan(d) & d == 1);
        allOnes = all(d == 1);
        allZeros = all(d == 0);
        % Part1: tuned: 1..(prefixes-1) responsive OR all ones & kruskal<0.05
        if (nOnes >= 1 && nOnes <= tuned_upper) || (allOnes && ~isnan(k) && k < 0.05)
            is_tuned(i) = true;
            med = clusters(i).(grp).medians;
            [~,pref_ind] = max(med);
            if ~isnan(med(pref_ind))
                tuned_pref_idx(i) = pref_ind;
            end
        elseif (allZeros && ~isnan(w) && w < 0.05 && ~isnan(med_overall) && med_overall > 0) || (allOnes && ~isnan(k) && k > 0.05)
            is_resp_not_tuned(i) = true;
        elseif allZeros && ~isnan(w) && w > 0.05 && ~isnan(k) && k > 0.05
            is_not_responsive(i) = true;
        end
    end
    counts.tuned = sum(is_tuned);
    counts.resp_not_tuned = sum(is_resp_not_tuned);
    counts.not_responsive = sum(is_not_responsive);
    counts.other = nClusters - (counts.tuned + counts.resp_not_tuned + counts.not_responsive);
    tuned_idx = find(is_tuned);
    pref_counts = zeros(1,numel(prefixes_all));
    for ii = tuned_idx
        idx = tuned_pref_idx(ii);
        if ~isnan(idx)
            pref_counts(idx) = pref_counts(idx) + 1;
        end
    end
end
