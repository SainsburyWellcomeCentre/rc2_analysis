% Plot baseline vs. response firing rates for motion cloud experiments
% Similar to single_cluster_svm_unity_plot.m but with 2-column plots instead of unity plots
% and filtering for specific motion cloud names

experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'motionClouds', 'passive_same_luminance_mc', 'stationary_vs_motion', 'single_cluster'};

% Motion cloud filtering
name_prefix = 'theta-0p785';

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

x_all                   = cell(1, length(probe_ids));
y_all                   = cell(1, length(probe_ids));
p_val                   = cell(1, length(probe_ids));
direction               = cell(1, length(probe_ids));
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
        kept_ids = [];
        if ~isempty(trials)
            if iscell(trials)
                trial_ids_list = cellfun(@(t) t.trial_id, trials);
            else
                trial_ids_list = trials.trial_id;
            end
            trial_ids_list = trial_ids_list(:)';
            keep_mask = false(size(trial_ids_list));
            for ti = 1:numel(trial_ids_list)
                tid = trial_ids_list(ti);
                if ~isempty(mc_sequence) && ~isempty(cloud_names) && tid>=1 && tid<=length(mc_sequence)
                    mc_id = mc_sequence(tid);
                    if mc_id>=1 && mc_id<=length(cloud_names)
                        cname = cloud_names{mc_id};
                        if (exist('startsWith','builtin') && startsWith(cname, name_prefix)) || ...
                           (~exist('startsWith','builtin') && strncmp(cname, name_prefix, length(name_prefix)))
                            keep_mask(ti) = true;
                        end
                    end
                end
            end
            kept_ids = trial_ids_list(keep_mask);
        end

        for kk = 1 : length(cluster_ids{ii})
            
            % Get baseline/response with their trial_ids and align
            [baseline_fr, base_ids] = data.stationary_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            [response_fr, resp_ids] = data.motion_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            [common_ids, ia, ib] = intersect(base_ids(:), resp_ids(:));
            baseline_fr = baseline_fr(ia);
            response_fr = response_fr(ib);

            if ~isempty(kept_ids)
                keep_common = ismember(common_ids, kept_ids);
                x_all{ii}{kk}{jj} = baseline_fr(keep_common);
                y_all{ii}{kk}{jj} = response_fr(keep_common);
            else
                x_all{ii}{kk}{jj} = [];
                y_all{ii}{kk}{jj} = [];
            end
            
            % Calculate significance on filtered pairs only
            xb2 = x_all{ii}{kk}{jj}(:);
            yr2 = y_all{ii}{kk}{jj}(:);
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
            p_val{ii}{kk}(jj) = pv;
            direction{ii}{kk}(jj) = dirn;
        end
    end
end

% Print summary
for ii = 1:length(probe_ids)
    data = ctl.load_formatted_data(probe_ids{ii});
    kept_counts = zeros(1, numel(trial_group_labels));
    for jj = 1:numel(trial_group_labels)
        if ~data.check_trial_group(trial_group_labels{jj})
            kept_counts(jj) = 0; continue
        end
        trials = data.get_trials_with_trial_group_label(trial_group_labels{jj});
        if isempty(trials)
            kept_counts(jj) = 0; continue
        end
        if iscell(trials)
            trial_ids_list = cellfun(@(t) t.trial_id, trials);
        else
            trial_ids_list = trials.trial_id;
        end
        keep_mask = false(size(trial_ids_list));
        for ti = 1:numel(trial_ids_list)
            tid = trial_ids_list(ti);
            if ~isempty(mc_sequence) && ~isempty(cloud_names) && tid>=1 && tid<=length(mc_sequence)
                mc_id = mc_sequence(tid);
                if mc_id>=1 && mc_id<=length(cloud_names)
                    cname = cloud_names{mc_id};
                    if (exist('startsWith','builtin') && startsWith(cname, name_prefix)) || ...
                       (~exist('startsWith','builtin') && strncmp(cname, name_prefix, length(name_prefix)))
                        keep_mask(ti) = true;
                    end
                end
            end
        end
        kept_counts(jj) = sum(keep_mask);
    end
    fprintf('Probe %s: kept trials per protocol_id [VT(1) V(2) T_Vstatic(3)] = [%d %d %d]\n', ...
        probe_ids{ii}, kept_counts(1), kept_counts(2), kept_counts(3));
end

%%
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    for jj = 1 : length(cluster_ids{ii})
        fprintf('Probe %s: plotting cluster_id=%d\n', probe_ids{ii}, cluster_ids{ii}(jj));
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        axs                     = gobjects(0);
        panel_mins              = [];
        panel_maxs              = [];

        for kk = 1 : length(trial_group_labels)
            
            pos         = plot_array.get_position(kk);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            axs(end+1)  = h_ax;
            
            % Get data for this trial group
            x_data = x_all{ii}{jj}{kk};
            y_data = y_all{ii}{jj}{kk};
            
            if isempty(x_data) || isempty(y_data)
                continue;
            end
            
            % Plot baseline (x=1) and response (x=2) as separate columns
            xb = x_data(:);
            yr = y_data(:);
            
            % Plot individual points
            scatter(h_ax, ones(numel(xb),1), xb, 50, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
            hold(h_ax, 'on');
            scatter(h_ax, 2*ones(numel(yr),1), yr, 50, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
            
            % Add connecting lines between paired baseline/response
            n_pairs = min(numel(xb), numel(yr));
            for i = 1:n_pairs
                plot(h_ax, [1 2], [xb(i) yr(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.3]);
            end
            
            % Add significance coloring for response points
            if p_val{ii}{jj}(kk) < 0.05
                if direction{ii}{jj}(kk) > 0
                    scatter(h_ax, 2*ones(numel(yr),1), yr, 80, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
                else
                    scatter(h_ax, 2*ones(numel(yr),1), yr, 80, 'b', 'filled', 'MarkerFaceAlpha', 0.8);
                end
            else
                scatter(h_ax, 2*ones(numel(yr),1), yr, 80, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.8);
            end
            
            xlim(h_ax, [0.5 2.5]);
            xlabel(h_ax, 'Condition');
            ylabel(h_ax, 'Firing Rate (Hz)');
            title(h_ax, trial_group_labels{kk});
            
            % Set x-axis labels
            set(h_ax, 'XTick', [1 2], 'XTickLabel', {'Baseline', 'Response'});
            
            % Store min/max for axis synchronization
            panel_mins(end+1) = min([xb; yr]);
            panel_maxs(end+1) = max([xb; yr]);
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