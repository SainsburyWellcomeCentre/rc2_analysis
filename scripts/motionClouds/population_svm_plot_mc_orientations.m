% population_svm_plot_mc_orientations.m
% Read single-cluster orientations CSV and produce population plots/CSVs
% Parts:
% 1) clusters responsive in 'V' with orientation tuning (at least one response_direction==1)
% 2) clusters responsive in 'V' but no orientation tuning (all response_direction==0 and wilcoxon<0.05)
% 3) clusters not responsive in 'V' but responsive & tuned in 'VT'
% 4) clusters not responsive in 'V' but responsive in 'VT' without tuning
%
% Saves PDFs and CSV summaries to the same figures folder as the input CSV.

input_csv = fullfile('C:','Users','lee','Documents','mvelez','figures','motionClouds','passive_same_luminance_mc','stationary_vs_motion_orientations','single_cluster','motion_cloud_orientations_analysis_summary.csv');
if ~exist(input_csv,'file')
    error('Input CSV not found: %s', input_csv);
end

T = readtable(input_csv);
% normalize some column types
if iscell(T.prefix), prefixes_all = unique(T.prefix,'stable'); else prefixes_all = unique(cellstr(T.prefix),'stable'); end
prefixes_all = cellstr(prefixes_all(:)');
trial_groups = {'V','VT'};
lightV = [0.6 0.85 1.0]; darkV = [0 0.2 0.6];
lightVT = [1.0 0.7 0.7]; darkVT = [0.6 0 0];
[cluster_keys, ia] = unique(T(:,{'probe_id','cluster_id'}),'rows','stable');
cluster_keys = cluster_keys(:,:);

nClusters = height(cluster_keys);

% precompute per-cluster, per-group, per-prefix medians and directions and wilcoxon
% store in struct array keyed by (i)
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
                % response_direction should be identical for all rows of that grouping; take mode/non-NaN
                rd = unique(T.response_direction(sel));
                rd = rd(~isnan(rd));
                if isempty(rd), clusters(i).(grp).direction(p) = 0; else clusters(i).(grp).direction(p) = rd(1); end
            else
                clusters(i).(grp).medians(p) = NaN;
                clusters(i).(grp).direction(p) = NaN;
            end
        end
        % wilcoxon per cluster/group (take first non-NaN if present)
        sel_all = strcmp(T.probe_id, pid) & (T.cluster_id == cid) & strcmp(T.trial_group, grp);
        if any(sel_all)
            if ismember('wilcoxon_p_bsl_vs_resp', T.Properties.VariableNames)
                wvals = unique(T.wilcoxon_p_bsl_vs_resp(sel_all));
                wvals = wvals(~isnan(wvals));
                if ~isempty(wvals), clusters(i).(grp).wilcoxon = wvals(1); end
            end
            if ismember('kruskal_p_orientation', T.Properties.VariableNames)
                kvals = unique(T.kruskal_p_orientation(sel_all));
                kvals = kvals(~isnan(kvals));
                if ~isempty(kvals), clusters(i).(grp).kruskal = kvals(1); end
            end
        end
    end
end

% small helper: rotate order so that preferred index ends up at position 2
rotate_to_pos2 = @(arr,k) circshift(arr, 2-k);

out_dir = fullfile('C:','Users','lee','Documents','mvelez','figures','motionClouds','passive_same_luminance_mc','stationary_vs_motion_orientations','population');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

% Part 1: clusters with in 'V' at least one orientation with response_direction == 1
part1_idx = false(1,nClusters);
for i=1:nClusters
    dV = clusters(i).V.direction;
    wV = clusters(i).V.wilcoxon;
    kV = clusters(i).V.kruskal;
    % count how many prefixes have response_direction == 1 (ignore NaNs)
    nOnesV = sum(~isnan(dV) & dV == 1);
    allOnesV = all(dV == 1);
    % Part1: one-to-three response_direction==1 OR all response_direction==1 AND kruskal_p_orientation<0.05
    if (nOnesV >= 1 && nOnesV <= (numel(prefixes_all)-1)) || (allOnesV && ~isnan(kV) && kV < 0.05)
        part1_idx(i) = true;
    end
end
part1_clusters = find(part1_idx);

    if ~isempty(part1_clusters)
    fig = figure('Visible','off'); hold on;
    summary_rows = {}; % use cell array to collect structs to avoid struct-field mismatch issues
    % arrays to collect rotated medians for across-cluster median
    medV_all = nan(numel(part1_clusters), numel(prefixes_all));
    medVT_all = nan(numel(part1_clusters), numel(prefixes_all));
    row_ctr = 0;
    lightV = [0.6 0.85 1.0]; darkV = [0 0.2 0.6];
    lightVT = [1.0 0.7 0.7]; darkVT = [0.6 0 0];
    for ii = 1:numel(part1_clusters)
        i = part1_clusters(ii);
        medV = clusters(i).V.medians;
        % pick orientation with highest median (ignore NaNs)
        [~,pref_ind] = max(medV); % if all NaN, pref_ind=1 but will be NaN
        if isnan(medV(pref_ind))
            continue
        end
        order_idx = rotate_to_pos2(1:numel(prefixes_all), pref_ind);
        medV_ord = medV(order_idx);
        medVT_ord = clusters(i).VT.medians(order_idx);
        x = 1:numel(prefixes_all);
        % collect per-cluster medians (we'll plot Q1/Q3 across clusters instead of individual traces)
        row_ctr = row_ctr + 1;
        medV_all(row_ctr, :) = medV_ord;
        medVT_all(row_ctr, :) = medVT_ord;
        % summary row
        row.probe_id = clusters(i).probe_id;
        row.cluster_id = clusters(i).cluster_id;
        row.pref = prefixes_all{pref_ind};
        % store medians as columns
        for p = 1:numel(prefixes_all)
            row.(sprintf('V_med_%d',p)) = medV_ord(p);
            row.(sprintf('VT_med_%d',p)) = medVT_ord(p);
        end
        summary_rows{end+1} = row;
    end
    % trim unused preallocated rows if any
    if row_ctr < size(medV_all,1)
        medV_all = medV_all(1:row_ctr,:);
        medVT_all = medVT_all(1:row_ctr,:);
    end
    % compute across-cluster medians (omit NaNs) and Q1/Q3
    if ~isempty(medV_all)
        nx = numel(prefixes_all);
        q1V = nan(1,nx); q3V = nan(1,nx); q1VT = nan(1,nx); q3VT = nan(1,nx);
        for xi = 1:nx
            v = medV_all(:,xi); v = v(~isnan(v)); if ~isempty(v), qq = quantile(v,[0.25 0.75]); q1V(xi)=qq(1); q3V(xi)=qq(2); end
            vt = medVT_all(:,xi); vt = vt(~isnan(vt)); if ~isempty(vt), qq = quantile(vt,[0.25 0.75]); q1VT(xi)=qq(1); q3VT(xi)=qq(2); end
        end
        x = 1:nx;
        if any(~isnan(q1V))
            fill([x fliplr(x)], [q1V fliplr(q3V)], lightV, 'FaceAlpha',0.25, 'EdgeColor','none');
        end
        if any(~isnan(q1VT))
            fill([x fliplr(x)], [q1VT fliplr(q3VT)], lightVT, 'FaceAlpha',0.25, 'EdgeColor','none');
        end
        medianV = median(medV_all, 1, 'omitnan');
        medianVT = median(medVT_all, 1, 'omitnan');
        % plot medians in dark colors and include in legend
        hV = plot(x, medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
        hVT = plot(x, medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
    end
    title(sprintf('Part1: clusters responsive in V with tuning (n=%d out of %d)', row_ctr, nClusters));
    xticks(1:numel(prefixes_all));
    % Label only the preferred-orientation position (x=2) as V_pref_ori
    labs = repmat({''}, 1, numel(prefixes_all));
    labs{2} = 'V_pref_ori';
    xticklabels(labs);
    xlabel('Orientation (rotated so pref at x=2)'); ylabel('Median(Response - Baseline) (Hz)');
    if exist('hV','var') && exist('hVT','var')
        legend([hV,hVT], {'V (median)','VT (median)'});
    else
        legend({'V (median)','VT (median)'});
    end
    % save figure
    pdfname = fullfile(out_dir, 'population_orientations_part1_V_tuned.pdf');
    set(fig,'PaperPositionMode','auto'); print(fig, pdfname, '-dpdf'); close(fig);
    % build table from summary_rows
    if ~isempty(summary_rows)
        S = struct2table([summary_rows{:}]);
        writetable(S, fullfile(out_dir,'population_orientations_part1_summary.csv'));
    end
else
    fprintf('Part1: no clusters found\n');
end

% Part 2: clusters responsive in V, but no orientation tuning
part2_idx = false(1,nClusters);
for i=1:nClusters
    dV = clusters(i).V.direction;
    wV = clusters(i).V.wilcoxon;
    kV = clusters(i).V.kruskal;
    medV_overall = median(clusters(i).V.medians, 'omitnan');
    allZerosV = all(dV == 0);
    allOnesV = all(dV == 1);
    % Part2: (all V directions == 0 AND wilcoxon<0.05) OR (all V directions == 1 AND kruskal>0.05)
    if (allZerosV && ~isnan(wV) && wV < 0.05 && ~isnan(medV_overall) && medV_overall > 0) || (allOnesV && ~isnan(kV) && kV > 0.05)
        part2_idx(i) = true;
    end
end
part2_clusters = find(part2_idx);

    if ~isempty(part2_clusters)
    fig = figure('Visible','off'); hold on;
    summary_rows = {}; % collect as cell array of structs
    medV_all = nan(numel(part2_clusters), numel(prefixes_all));
    medVT_all = nan(numel(part2_clusters), numel(prefixes_all));
    row_ctr = 0;
    lightV = [0.6 0.85 1.0]; darkV = [0 0.2 0.6];
    lightVT = [1.0 0.7 0.7]; darkVT = [0.6 0 0];
    for ii = 1:numel(part2_clusters)
        i = part2_clusters(ii);
        medV = clusters(i).V.medians;
        medVT = clusters(i).VT.medians;
        x = 1:numel(prefixes_all);
        % collect per-cluster medians (we'll plot Q1/Q3 across clusters instead of individual traces)
        row_ctr = row_ctr + 1;
        medV_all(row_ctr,:) = medV;
        medVT_all(row_ctr,:) = medVT;
        % summary
        row.probe_id = clusters(i).probe_id;
        row.cluster_id = clusters(i).cluster_id;
        for p = 1:numel(prefixes_all)
            row.(sprintf('V_med_%d',p)) = medV(p);
            row.(sprintf('VT_med_%d',p)) = medVT(p);
        end
        summary_rows{end+1} = row;
    end
    if row_ctr < size(medV_all,1)
        medV_all = medV_all(1:row_ctr,:);
        medVT_all = medVT_all(1:row_ctr,:);
    end
    % compute medians and Q1/Q3 across clusters
    nx = numel(prefixes_all);
    q1V = nan(1,nx); q3V = nan(1,nx); q1VT = nan(1,nx); q3VT = nan(1,nx);
    for xi = 1:nx
        v = medV_all(:,xi); v = v(~isnan(v)); if ~isempty(v), qq = quantile(v,[0.25 0.75]); q1V(xi)=qq(1); q3V(xi)=qq(2); end
        vt = medVT_all(:,xi); vt = vt(~isnan(vt)); if ~isempty(vt), qq = quantile(vt,[0.25 0.75]); q1VT(xi)=qq(1); q3VT(xi)=qq(2); end
    end
    x = 1:nx;
    if any(~isnan(q1V)), fill([x fliplr(x)], [q1V fliplr(q3V)], lightV, 'FaceAlpha',0.25, 'EdgeColor','none'); end
    if any(~isnan(q1VT)), fill([x fliplr(x)], [q1VT fliplr(q3VT)], lightVT, 'FaceAlpha',0.25, 'EdgeColor','none'); end
    medianV = median(medV_all,1,'omitnan');
    medianVT = median(medVT_all,1,'omitnan');
    hV = plot(x, medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
    hVT = plot(x, medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
    title(sprintf('Part2: V responsive but no orientation tuning (n=%d out of %d)', row_ctr, nClusters));
    xticks(1:numel(prefixes_all)); xticklabels(prefixes_all);
    xlabel('Orientation'); ylabel('Median(Response - Baseline) (Hz)'); legend([hV,hVT],{'V (median)','VT (median)'});
    pdfname = fullfile(out_dir, 'population_orientations_part2_V_responsive_no_tuning.pdf');
    set(fig,'PaperPositionMode','auto'); print(fig, pdfname, '-dpdf'); close(fig);
    if ~isempty(summary_rows)
        S = struct2table([summary_rows{:}]);
        writetable(S, fullfile(out_dir,'population_orientations_part2_summary.csv'));
    end
else
    fprintf('Part2: no clusters found\n');
end

% Part 3: not responsive in V, but responsive & tuned in VT
part3_idx = false(1,nClusters);
for i=1:nClusters
    dV = clusters(i).V.direction;
    wV = clusters(i).V.wilcoxon;
    dVT = clusters(i).VT.direction;
    wVT = clusters(i).VT.wilcoxon;
    kVT = clusters(i).VT.kruskal;
    allZerosV = all(dV == 0);
    % VT tuned count
    nOnesVT = sum(~isnan(dVT) & dVT == 1);
    allOnesVT = all(dVT == 1);
    % Part3: V all zeros & wilcoxon>0.05 AND (VT has 1-3 ones OR all VT ones & kruskal<0.05)
    if allZerosV && (~isnan(wV) && wV > 0.05)
        if (nOnesVT >= 1 && nOnesVT <= (numel(prefixes_all)-1)) || (allOnesVT && ~isnan(kVT) && kVT < 0.05)
            part3_idx(i) = true;
        end
    end
end
part3_clusters = find(part3_idx);

    if ~isempty(part3_clusters)
    fig = figure('Visible','off'); hold on;
    summary_rows = {}; % collect as cell array of structs
    medVT_all = nan(numel(part3_clusters), numel(prefixes_all));
    medV_all = nan(numel(part3_clusters), numel(prefixes_all));
    row_ctr = 0;
    lightV = [0.6 0.85 1.0]; darkV = [0 0.2 0.6];
    lightVT = [1.0 0.7 0.7]; darkVT = [0.6 0 0];
    for ii = 1:numel(part3_clusters)
        i = part3_clusters(ii);
        medVT = clusters(i).VT.medians;
        [~,pref_ind] = max(medVT);
        if isnan(medVT(pref_ind)), continue; end
        order_idx = rotate_to_pos2(1:numel(prefixes_all), pref_ind);
        medVT_ord = medVT(order_idx);
        medV_ord = clusters(i).V.medians(order_idx);
        x = 1:numel(prefixes_all);
        % collect per-cluster medians (we'll plot Q1/Q3 across clusters instead of individual traces)
        row_ctr = row_ctr + 1;
        medVT_all(row_ctr,:) = medVT_ord;
        medV_all(row_ctr,:) = medV_ord;
        row.probe_id = clusters(i).probe_id;
        row.cluster_id = clusters(i).cluster_id;
        row.pref = prefixes_all{pref_ind};
        for p = 1:numel(prefixes_all)
            row.(sprintf('VT_med_%d',p)) = medVT_ord(p);
            row.(sprintf('V_med_%d',p)) = medV_ord(p);
        end
        summary_rows{end+1} = row;
    end
    if row_ctr < size(medVT_all,1)
        medVT_all = medVT_all(1:row_ctr,:);
        medV_all = medV_all(1:row_ctr,:);
    end
    % compute medians and Q1/Q3 across clusters
    nx = numel(prefixes_all);
    q1V = nan(1,nx); q3V = nan(1,nx); q1VT = nan(1,nx); q3VT = nan(1,nx);
    for xi = 1:nx
        vt = medVT_all(:,xi); vt = vt(~isnan(vt)); if ~isempty(vt), qq = quantile(vt,[0.25 0.75]); q1VT(xi)=qq(1); q3VT(xi)=qq(2); end
        v = medV_all(:,xi); v = v(~isnan(v)); if ~isempty(v), qq = quantile(v,[0.25 0.75]); q1V(xi)=qq(1); q3V(xi)=qq(2); end
    end
    x = 1:nx;
    if any(~isnan(q1V)), fill([x fliplr(x)], [q1V fliplr(q3V)], lightV, 'FaceAlpha',0.25, 'EdgeColor','none'); end
    if any(~isnan(q1VT)), fill([x fliplr(x)], [q1VT fliplr(q3VT)], lightVT, 'FaceAlpha',0.25, 'EdgeColor','none'); end
    medianVT = median(medVT_all,1,'omitnan');
    medianV = median(medV_all,1,'omitnan');
    hVT = plot(x, medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
    hV = plot(x, medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
    title(sprintf('Part3: VT tuned but V not responsive (n=%d out of %d)', row_ctr, nClusters));
    xticks(1:numel(prefixes_all));
    % Label only the preferred-orientation position (x=2) as VT_pref_ori
    labs = repmat({''}, 1, numel(prefixes_all));
    labs{2} = 'VT_pref_ori';
    xticklabels(labs);
    xlabel('Orientation (rotated so VT pref at x=2)'); ylabel('Median(Response - Baseline) (Hz)'); legend([hVT,hV],{'VT (median)','V (median)'});
    pdfname = fullfile(out_dir, 'population_orientations_part3_V_not_responsive_VT_tuned.pdf');
    set(fig,'PaperPositionMode','auto'); print(fig, pdfname, '-dpdf'); close(fig);
    if ~isempty(summary_rows)
        S = struct2table([summary_rows{:}]);
        writetable(S, fullfile(out_dir,'population_orientations_part3_summary.csv'));
    end
else
    fprintf('Part3: no clusters found\n');
end

% Part 4: V not responsive, VT responsive but not orientation tuned
part4_idx = false(1,nClusters);
for i=1:nClusters
    dV = clusters(i).V.direction;
    wV = clusters(i).V.wilcoxon;
    dVT = clusters(i).VT.direction;
    wVT = clusters(i).VT.wilcoxon;
    kVT = clusters(i).VT.kruskal;
    medVT_overall = median(clusters(i).VT.medians, 'omitnan');
    allZerosV = all(dV == 0);
    allZerosVT = all(dVT == 0);
    allOnesVT = all(dVT == 1);
    % Part4: (V all zeros & wV>0.05) AND ((VT all zeros & wVT<0.05) OR (all VT ones & kruskal>0.05))
    if allZerosV && (~isnan(wV) && wV > 0.05)
        if (allZerosVT && ~isnan(wVT) && wVT < 0.05 && ~isnan(medVT_overall) && medVT_overall > 0) || (allOnesVT && ~isnan(kVT) && kVT > 0.05)
            part4_idx(i) = true;
        end
    end
end
part4_clusters = find(part4_idx);

    if ~isempty(part4_clusters)
    fig = figure('Visible','off'); hold on;
    summary_rows = {}; % collect as cell array of structs
    medVT_all = nan(numel(part4_clusters), numel(prefixes_all));
    medV_all = nan(numel(part4_clusters), numel(prefixes_all));
    row_ctr = 0;
    lightV = [0.6 0.85 1.0]; darkV = [0 0.2 0.6];
    lightVT = [1.0 0.7 0.7]; darkVT = [0.6 0 0];
    for ii = 1:numel(part4_clusters)
        i = part4_clusters(ii);
        medVT = clusters(i).VT.medians;
        medV = clusters(i).V.medians;
        x = 1:numel(prefixes_all);
        % collect per-cluster medians (we'll plot Q1/Q3 across clusters instead of individual traces)
        row_ctr = row_ctr + 1;
        medVT_all(row_ctr,:) = medVT;
        medV_all(row_ctr,:) = medV;
        row.probe_id = clusters(i).probe_id;
        row.cluster_id = clusters(i).cluster_id;
        for p = 1:numel(prefixes_all)
            row.(sprintf('VT_med_%d',p)) = medVT(p);
            row.(sprintf('V_med_%d',p)) = medV(p);
        end
        summary_rows{end+1} = row;
    end
    if row_ctr < size(medVT_all,1)
        medVT_all = medVT_all(1:row_ctr,:);
        medV_all = medV_all(1:row_ctr,:);
    end
    % compute medians and Q1/Q3 across clusters
    nx = numel(prefixes_all);
    q1V = nan(1,nx); q3V = nan(1,nx); q1VT = nan(1,nx); q3VT = nan(1,nx);
    for xi = 1:nx
        vt = medVT_all(:,xi); vt = vt(~isnan(vt)); if ~isempty(vt), qq = quantile(vt,[0.25 0.75]); q1VT(xi)=qq(1); q3VT(xi)=qq(2); end
        v = medV_all(:,xi); v = v(~isnan(v)); if ~isempty(v), qq = quantile(v,[0.25 0.75]); q1V(xi)=qq(1); q3V(xi)=qq(2); end
    end
    x = 1:nx;
    if any(~isnan(q1V)), fill([x fliplr(x)], [q1V fliplr(q3V)], lightV, 'FaceAlpha',0.25, 'EdgeColor','none'); end
    if any(~isnan(q1VT)), fill([x fliplr(x)], [q1VT fliplr(q3VT)], lightVT, 'FaceAlpha',0.25, 'EdgeColor','none'); end
    medianVT = median(medVT_all,1,'omitnan');
    medianV = median(medV_all,1,'omitnan');
    hVT = plot(x, medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
    hV = plot(x, medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
    title(sprintf('Part4: VT responsive no tuning, V not responsive (n=%d out of %d)', row_ctr, nClusters));
    xticks(1:numel(prefixes_all)); xticklabels(prefixes_all);
    xlabel('Orientation'); ylabel('Median(Response - Baseline) (Hz)'); legend([hVT,hV],{'VT (median)','V (median)'});
    pdfname = fullfile(out_dir, 'population_orientations_part4_VT_responsive_no_tuning.pdf');
    set(fig,'PaperPositionMode','auto'); print(fig, pdfname, '-dpdf'); close(fig);
    if ~isempty(summary_rows)
        S = struct2table([summary_rows{:}]);
        writetable(S, fullfile(out_dir,'population_orientations_part4_summary.csv'));
    end
else
    fprintf('Part4: no clusters found\n');
end

fprintf('Population plots and CSVs written to: %s\n', out_dir);

% Diagnostic: report clusters that were included in more than one Part
cluster_parts = cell(nClusters,1);
for i=1:nClusters
    if exist('part1_idx','var') && part1_idx(i), cluster_parts{i} = [cluster_parts{i}, 1]; end
    if exist('part2_idx','var') && part2_idx(i), cluster_parts{i} = [cluster_parts{i}, 2]; end
    if exist('part3_idx','var') && part3_idx(i), cluster_parts{i} = [cluster_parts{i}, 3]; end
    if exist('part4_idx','var') && part4_idx(i), cluster_parts{i} = [cluster_parts{i}, 4]; end
end
multi_idx = find(cellfun(@numel, cluster_parts) > 1);
if ~isempty(multi_idx)
    fprintf('Diagnostic: %d clusters assigned to multiple parts\n', numel(multi_idx));
    diag_rows = struct(); diag_rows = repmat(diag_rows, numel(multi_idx), 1);
    for k = 1:numel(multi_idx)
        ii = multi_idx(k);
        parts_str = sprintf('%d,', cluster_parts{ii}); parts_str = parts_str(1:end-1);
        fprintf('  probe %s cluster %d -> parts: %s\n', clusters(ii).probe_id, clusters(ii).cluster_id, parts_str);
        diag_rows(k).probe_id = clusters(ii).probe_id;
        diag_rows(k).cluster_id = clusters(ii).cluster_id;
        diag_rows(k).parts = parts_str;
    end
    try
        Tdiag = struct2table(diag_rows);
        writetable(Tdiag, fullfile(out_dir,'cluster_parts_diagnostic_orientations.csv'));
    catch
        warning('Could not write diagnostic CSV for orientations');
    end
else
    fprintf('Diagnostic: no clusters assigned to multiple parts\n');
end
