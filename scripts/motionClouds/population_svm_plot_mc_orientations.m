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

% helper: unique clusters (probe_id + cluster_id)
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
        if any(sel_all) && ismember('wilcoxon_p_bsl_vs_resp', T.Properties.VariableNames)
            wvals = unique(T.wilcoxon_p_bsl_vs_resp(sel_all));
            wvals = wvals(~isnan(wvals));
            if ~isempty(wvals), clusters(i).(grp).wilcoxon = wvals(1); end
        end
    end
end

% small helper: rotate order so that preferred index ends up at position 2
rotate_to_pos2 = @(arr,k) circshift(arr, 2-k);

out_dir = fileparts(input_csv);

% Part 1: clusters with in 'V' at least one orientation with response_direction == 1
part1_idx = false(1,nClusters);
for i=1:nClusters
    % check V directions
    dV = clusters(i).V.direction;
    if any(dV == 1)
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
        % plot individual clusters in light colors, hide from legend
        lightV = [0.6 0.85 1.0]; darkV = [0 0.2 0.6];
        lightVT = [1.0 0.7 0.7]; darkVT = [0.6 0 0];
        plot(x, medV_ord, '-o', 'Color', lightV, 'MarkerFaceColor', lightV, 'LineWidth',0.8, 'HandleVisibility','off');
        plot(x, medVT_ord, '-o', 'Color', lightVT, 'MarkerFaceColor', lightVT, 'LineWidth',0.8, 'HandleVisibility','off');
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
    % compute across-cluster medians (omit NaNs)
    if ~isempty(medV_all)
        medianV = median(medV_all, 1, 'omitnan');
        medianVT = median(medVT_all, 1, 'omitnan');
        % plot medians in dark colors and include in legend
        hV = plot(1:numel(prefixes_all), medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
        hVT = plot(1:numel(prefixes_all), medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
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
    % require V non-tuned (all directions 0), significant responsiveness (wilcoxon<0.05)
    % and positive median difference across V prefixes
    medV_overall = median(clusters(i).V.medians, 'omitnan');
    if all(dV == 0) && ~isnan(wV) && wV < 0.05 && ~isnan(medV_overall) && medV_overall > 0
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
        plot(x, medV, '-o', 'Color', lightV, 'MarkerFaceColor', lightV, 'LineWidth',0.8, 'HandleVisibility','off');
        plot(x, medVT, '-o', 'Color', lightVT, 'MarkerFaceColor', lightVT, 'LineWidth',0.8, 'HandleVisibility','off');
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
    medianV = median(medV_all,1,'omitnan');
    medianVT = median(medVT_all,1,'omitnan');
    hV = plot(1:numel(prefixes_all), medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
    hVT = plot(1:numel(prefixes_all), medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
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
    if all(dV == 0) && (~isnan(wV) && wV > 0.05) && any(dVT == 1)
        part3_idx(i) = true;
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
        plot(x, medVT_ord, '-o', 'Color', lightVT, 'MarkerFaceColor', lightVT, 'LineWidth',0.8, 'HandleVisibility','off');
        plot(x, medV_ord, '-o', 'Color', lightV, 'MarkerFaceColor', lightV, 'LineWidth',0.8, 'HandleVisibility','off');
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
    medianVT = median(medVT_all,1,'omitnan');
    medianV = median(medV_all,1,'omitnan');
    hVT = plot(1:numel(prefixes_all), medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
    hV = plot(1:numel(prefixes_all), medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
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
    % require V not responsive, VT responsive without tuning, and positive median in VT
    medVT_overall = median(clusters(i).VT.medians, 'omitnan');
    if all(dV == 0) && (~isnan(wV) && wV > 0.05) && all(dVT == 0) && (~isnan(wVT) && wVT < 0.05) && ~isnan(medVT_overall) && medVT_overall > 0
        part4_idx(i) = true;
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
        plot(x, medVT, '-o', 'Color', lightVT, 'MarkerFaceColor', lightVT, 'LineWidth',0.8, 'HandleVisibility','off');
        plot(x, medV, '-o', 'Color', lightV, 'MarkerFaceColor', lightV, 'LineWidth',0.8, 'HandleVisibility','off');
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
    medianVT = median(medVT_all,1,'omitnan');
    medianV = median(medV_all,1,'omitnan');
    hVT = plot(1:numel(prefixes_all), medianVT, '-o', 'Color', darkVT, 'MarkerFaceColor', darkVT, 'LineWidth',1.5);
    hV = plot(1:numel(prefixes_all), medianV, '-o', 'Color', darkV, 'MarkerFaceColor', darkV, 'LineWidth',1.5);
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
