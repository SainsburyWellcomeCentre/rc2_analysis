% Plot baseline vs. response firing rates for motion cloud experiments
% Similar to single_cluster_svm_unity_plot.m but with 2-column plots instead of unity plots
% and filtering for specific motion cloud names

experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'motionClouds', 'passive_same_luminance_mc', 'stationary_vs_motion_spatial_frequencies', 'single_cluster'};

% Parameters: Motion cloud trial filters

name_prefix1 = 'sf00p003';
name_prefix2 = 'sf00p006';
exclude_substrings_prefix2 = {'VX0p000'}; % e.g. exclude zero-velocity variants
name_prefix3 = 'sf00p012';

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
cluster_ids             = cell(1, length(probe_ids));
region_str              = cell(1, length(probe_ids));

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
            for ti = 1:numel(common_ids)
                tid = common_ids(ti);
                if ~isempty(mc_sequence) && ~isempty(cloud_names) && tid>=1 && tid<=length(mc_sequence)
                    mc_id = mc_sequence(tid);
                    if mc_id>=1 && mc_id<=length(cloud_names)
                        cname = cloud_names{mc_id};
                        if (exist('contains','builtin') && contains(cname, name_prefix1)) || ...
                           (~exist('contains','builtin') && ~isempty(strfind(cname, name_prefix1)))
                            keep1(ti) = true;
                        end
                        if ( (exist('contains','builtin') && contains(cname, name_prefix2)) || ...
                             (~exist('contains','builtin') && ~isempty(strfind(cname, name_prefix2))) )
                            if ~contains_any(cname, exclude_substrings_prefix2)
                                keep2(ti) = true;
                            end
                        end
                        if (exist('contains','builtin') && contains(cname, name_prefix3)) || ...
                           (~exist('contains','builtin') && ~isempty(strfind(cname, name_prefix3)))
                            keep3(ti) = true;
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
            
            % Calculate significance on filtered pairs (prefix1)
            xb2 = x_all1{ii}{kk}{jj}(:);
            yr2 = y_all1{ii}{kk}{jj}(:);
            good = isfinite(xb2) & isfinite(yr2);
            xb2 = xb2(good); yr2 = yr2(good);
            pv = NaN; dirn = 0;
            if numel(xb2) >= 2
                if exist('signrank','file') == 2
                    try
                        pv = signrank(xb2, yr2);
                    catch
                        try
                            [~, pv] = ttest(xb2, yr2);
                        catch
                            pv = NaN;
                        end
                    end
                else
                    try
                        [~, pv] = ttest(xb2, yr2);
                    catch
                        pv = NaN;
                    end
                end
                d = median(yr2 - xb2);
                if d > 0
                    dirn = 1;
                elseif d < 0
                    dirn = -1;
                else
                    dirn = 0;
                end
            end
            p_val1{ii}{kk}(jj) = pv;
            direction1{ii}{kk}(jj) = dirn;

            % Calculate significance on filtered pairs (prefix2)
            xb3 = x_all2{ii}{kk}{jj}(:);
            yr3 = y_all2{ii}{kk}{jj}(:);
            good = isfinite(xb3) & isfinite(yr3);
            xb3 = xb3(good); yr3 = yr3(good);
            pv = NaN; dirn = 0;
            if numel(xb3) >= 2
                if exist('signrank','file') == 2
                    try
                        pv = signrank(xb3, yr3);
                    catch
                        try
                            [~, pv] = ttest(xb3, yr3);
                        catch
                            pv = NaN;
                        end
                    end
                else
                    try
                        [~, pv] = ttest(xb3, yr3);
                    catch
                        pv = NaN;
                    end
                end
                d = median(yr3 - xb3);
                if d > 0
                    dirn = 1;
                elseif d < 0
                    dirn = -1;
                else
                    dirn = 0;
                end
            end
            p_val2{ii}{kk}(jj) = pv;
            direction2{ii}{kk}(jj) = dirn;

            % Calculate significance on filtered pairs (prefix3)
            xb4 = x_all3{ii}{kk}{jj}(:);
            yr4 = y_all3{ii}{kk}{jj}(:);
            good = isfinite(xb4) & isfinite(yr4);
            xb4 = xb4(good); yr4 = yr4(good);
            pv = NaN; dirn = 0;
            if numel(xb4) >= 2
                if exist('signrank','file') == 2
                    try
                        pv = signrank(xb4, yr4);
                    catch
                        try
                            [~, pv] = ttest(xb4, yr4);
                        catch
                            pv = NaN;
                        end
                    end
                else
                    try
                        [~, pv] = ttest(xb4, yr4);
                    catch
                        pv = NaN;
                    end
                end
                d = median(yr4 - xb4);
                if d > 0
                    dirn = 1;
                elseif d < 0
                    dirn = -1;
                else
                    dirn = 0;
                end
            end
            p_val3{ii}{kk}(jj) = pv;
            direction3{ii}{kk}(jj) = dirn;
        end
    end
end

% Print summary
for ii = 1:length(probe_ids)
    data = ctl.load_formatted_data(probe_ids{ii});
    kept_counts1 = zeros(1, numel(trial_group_labels));
    kept_counts2 = zeros(1, numel(trial_group_labels));
    kept_counts3 = zeros(1, numel(trial_group_labels));
    for jj = 1:numel(trial_group_labels)
        if ~data.check_trial_group(trial_group_labels{jj})
            kept_counts1(jj) = 0; kept_counts2(jj) = 0; kept_counts3(jj) = 0; continue
        end
        trials = data.get_trials_with_trial_group_label(trial_group_labels{jj});
        if isempty(trials)
            kept_counts1(jj) = 0; kept_counts2(jj) = 0; kept_counts3(jj) = 0; continue
        end
        if iscell(trials)
            trial_ids_list = cellfun(@(t) t.trial_id, trials);
        else
            trial_ids_list = trials.trial_id;
        end
        keep_mask1 = false(size(trial_ids_list));
        keep_mask2 = false(size(trial_ids_list));
        keep_mask3 = false(size(trial_ids_list));
        for ti = 1:numel(trial_ids_list)
            tid = trial_ids_list(ti);
            if ~isempty(mc_sequence) && ~isempty(cloud_names) && tid>=1 && tid<=length(mc_sequence)
                mc_id = mc_sequence(tid);
                if mc_id>=1 && mc_id<=length(cloud_names)
                    cname = cloud_names{mc_id};
                    if (exist('contains','builtin') && contains(cname, name_prefix1)) || ...
                       (~exist('contains','builtin') && ~isempty(strfind(cname, name_prefix1)))
                        keep_mask1(ti) = true;
                    end
                    if ( (exist('contains','builtin') && contains(cname, name_prefix2)) || ...
                         (~exist('contains','builtin') && ~isempty(strfind(cname, name_prefix2))) )
                        if ~contains_any(cname, exclude_substrings_prefix2)
                            keep_mask2(ti) = true;
                        end
                    end
                    if (exist('contains','builtin') && contains(cname, name_prefix3)) || ...
                       (~exist('contains','builtin') && ~isempty(strfind(cname, name_prefix3)))
                        keep_mask3(ti) = true;
                    end
                end
            end
        end
        kept_counts1(jj) = sum(keep_mask1);
        kept_counts2(jj) = sum(keep_mask2);
        kept_counts3(jj) = sum(keep_mask3);
    end
    fprintf('Probe %s: kept trials p1=%s [%d %d %d]; p2=%s [%d %d %d]; p3=%s [%d %d %d]\n', ...
        probe_ids{ii}, name_prefix1, kept_counts1(1), kept_counts1(2), kept_counts1(3), ...
        name_prefix2, kept_counts2(1), kept_counts2(2), kept_counts2(3), ...
        name_prefix3, kept_counts3(1), kept_counts3(2), kept_counts3(3));
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
p_value_list = [];
direction_list = [];

for ii = 1 : length(probe_ids)
    for jj = 1 : length(cluster_ids{ii})
        for kk = 1 : length(trial_group_labels)
            % Get data for this combination
            xb1 = x_all1{ii}{jj}{kk};
            yr1 = y_all1{ii}{jj}{kk};
            xb2 = x_all2{ii}{jj}{kk};
            yr2 = y_all2{ii}{jj}{kk};
            xb3 = x_all3{ii}{jj}{kk};
            yr3 = y_all3{ii}{jj}{kk};
            
            % Add data for each prefix
            prefixes = {name_prefix1, name_prefix2, name_prefix3};
            baseline_data = {xb1, xb2, xb3};
            response_data = {yr1, yr2, yr3};
            
            for p = 1:3
                baseline = baseline_data{p};
                response = response_data{p};
                
                if ~isempty(baseline) && ~isempty(response)
                    n_trials = min(length(baseline), length(response));
                    for t = 1:n_trials
                        probe_id_list{end+1} = probe_ids{ii};
                        cluster_id_list(end+1) = cluster_ids{ii}(jj);
                        trial_group_list{end+1} = trial_group_labels{kk};
                        prefix_list{end+1} = prefixes{p};
                        trial_number_list(end+1) = t;
                        baseline_fr_list(end+1) = baseline(t);
                        response_fr_list(end+1) = response(t);
                        difference_list(end+1) = response(t) - baseline(t);
                        if p == 1
                            p_value_list(end+1) = p_val1{ii}{jj}(kk);
                            direction_list(end+1) = direction1{ii}{jj}(kk);
                        elseif p == 2
                            p_value_list(end+1) = p_val2{ii}{jj}(kk);
                            direction_list(end+1) = direction2{ii}{jj}(kk);
                        elseif p == 3
                            p_value_list(end+1) = p_val3{ii}{jj}(kk);
                            direction_list(end+1) = direction3{ii}{jj}(kk);
                        end
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
                     p_value_list', direction_list');
    
    csv_filename = fullfile(ctl.figs.curr_dir, 'motion_cloud_spatial_frequencies_analysis_summary.csv');
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
            % and for prefix3
            xb3p = x_all3{ii}{jj}{kk};
            yr3p = y_all3{ii}{jj}{kk};
            
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

            xlim(h_ax, [0.5 6.5]);
            xlabel(h_ax, 'Condition');
            ylabel(h_ax, 'Firing Rate (Hz)');
            title(h_ax, trial_group_labels{kk});
            
            % Set x-axis labels
            set(h_ax, 'XTick', [1 2 3 4 5 6], 'XTickLabel', {'B1','R1','B2','R2','B3','R3'});
            
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
            
            xlim(h_ax2, [0.5 3.5]);
            xlabel(h_ax2, 'Motion Cloud Group');
            ylabel(h_ax2, 'Response - Baseline (Hz)');
            title(h_ax2, [trial_group_labels{kk} ' (Differences)']);
            
            % Set x-axis labels for differences
            set(h_ax2, 'XTick', [1 2 3], 'XTickLabel', {'P1','P2','P3'});
            
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
