% population_best_motion_cloud.m
% For each individual motion cloud stimulus, compute the population-level
% median (response - baseline) across all VISp clusters and all probes,
% using V trials only.
% Ranks all clouds by this metric and saves a bar chart + CSV.

experiment_groups = {'passive_same_luminance_mc'};
trial_group_label = 'V';

% Load motion cloud sequence and cloud names (goggles setup)
mc_sequence  = [];
cloud_names  = {};
proto_seq_path = fullfile('D:\mvelez\mateoData_mc', 'motion_clouds_goggles_sequence_260420.mat');
if exist(proto_seq_path, 'file')
    P = load(proto_seq_path);
    if isfield(P, 'presentation_sequence')
        mc_sequence = P.presentation_sequence;
    else
        fns = fieldnames(P);
        for i = 1:numel(fns)
            v = P.(fns{i});
            if isnumeric(v) && (isvector(v) || ismatrix(v))
                mc_sequence = v; break;
            end
        end
    end
end

folders_path = fullfile('D:\mvelez\mateoData_mc', 'image_folders.mat');
if exist(folders_path, 'file')
    S = load(folders_path);
    fns = fieldnames(S);
    for i = 1:numel(fns)
        v = S.(fns{i});
        if iscell(v)
            cloud_names = v; break;
        elseif isstring(v)
            cloud_names = cellstr(v(:)); break;
        elseif ischar(v)
            cloud_names = cellstr(v); break;
        elseif isstruct(v) && isfield(v, 'name')
            try cloud_names = {v.name}; break; catch, end
        end
    end
end
if ~isempty(cloud_names), cloud_names = cloud_names(:)'; end

if isempty(mc_sequence) || isempty(cloud_names)
    error('Could not load mc_sequence or cloud_names. Check paths.');
end

% Accumulate per-cloud differences across all probes and clusters
% containers.Map: cloud_name -> vector of (response - baseline) values
cloud_diffs = containers.Map('KeyType', 'char', 'ValueType', 'any');

ctl       = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});

for ii = 1:length(probe_ids)
    data        = ctl.load_formatted_data(probe_ids{ii});
    cluster_ids = data.VISp_cluster_ids;

    if ~data.check_trial_group(trial_group_label)
        fprintf('Probe %s: no %s trials, skipping.\n', probe_ids{ii}, trial_group_label);
        continue;
    end

    for kk = 1:length(cluster_ids)
        [baseline_fr, base_ids] = data.stationary_fr_for_trial_group(cluster_ids(kk), trial_group_label);
        [response_fr, resp_ids] = data.motion_fr_for_trial_group(cluster_ids(kk), trial_group_label);

        if isempty(baseline_fr) || isempty(response_fr), continue; end

        [common_ids, ia, ib] = intersect(base_ids(:), resp_ids(:));
        baseline_fr = baseline_fr(ia);
        response_fr = response_fr(ib);

        for ti = 1:numel(common_ids)
            tid = common_ids(ti);
            if tid < 1 || tid > length(mc_sequence), continue; end
            mc_id = mc_sequence(tid);
            if mc_id < 1 || mc_id > length(cloud_names), continue; end
            cname    = cloud_names{mc_id};
            diff_val = response_fr(ti) - baseline_fr(ti);
            if isnan(diff_val), continue; end
            if isKey(cloud_diffs, cname)
                cloud_diffs(cname) = [cloud_diffs(cname), diff_val];
            else
                cloud_diffs(cname) = diff_val;
            end
        end
    end
    fprintf('Probe %s: done.\n', probe_ids{ii});
end

% Compute median per cloud
all_clouds = keys(cloud_diffs);
n_clouds   = numel(all_clouds);
med_diffs  = nan(1, n_clouds);
n_obs      = zeros(1, n_clouds);
for ci = 1:n_clouds
    vals          = cloud_diffs(all_clouds{ci});
    med_diffs(ci) = median(vals, 'omitnan');
    n_obs(ci)     = numel(vals);
end

% Sort descending
[med_diffs_sorted, sort_idx] = sort(med_diffs, 'descend');
clouds_sorted = all_clouds(sort_idx);
n_obs_sorted  = n_obs(sort_idx);

% Print top 10 to console
fprintf('\n--- Top 10 motion clouds (V trials, population median response-baseline) ---\n');
for ci = 1:min(10, n_clouds)
    fprintf('%2d. %-80s  median=%+.3f Hz  (n=%d)\n', ...
        ci, clouds_sorted{ci}, med_diffs_sorted(ci), n_obs_sorted(ci));
end
fprintf('\nBest cloud: %s\n  median = %.3f Hz\n', clouds_sorted{1}, med_diffs_sorted(1));

% Output directory
out_dir = fullfile('C:', 'Users', 'lee', 'Documents', 'mvelez', 'figures', ...
    'motionClouds_mouseGoggles', 'passive_same_luminance_mc', ...
    'best_motion_clouds');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Save CSV (all clouds, sorted best to worst)
T = table(clouds_sorted(:), med_diffs_sorted(:), n_obs_sorted(:), ...
    'VariableNames', {'cloud_name', 'median_diff_Hz', 'n_observations'});
csv_out = fullfile(out_dir, 'population_best_motion_cloud_V.csv');
writetable(T, csv_out);
fprintf('CSV saved to: %s\n', csv_out);

% Horizontal bar chart (best cloud at top)
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 30 max(15, n_clouds*0.5)]);
barh(fliplr(med_diffs_sorted));
yticks(1:n_clouds);
yticklabels(fliplr(clouds_sorted));
xlabel('Median (Response - Baseline) Hz');
title(sprintf('Population response per motion cloud (%s trials)', trial_group_label));
xline(0, 'k--');
set(fig, 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape');
pdf_out = fullfile(out_dir, 'population_best_motion_cloud_V.pdf');
print(fig, pdf_out, '-dpdf');
close(fig);
fprintf('Plot saved to: %s\n', pdf_out);
