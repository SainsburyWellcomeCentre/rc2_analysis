function data = read_ks_data(data_dir, type)


npy_types = {'amplitudes', 'channel_map', 'channel_positions', ...
    'mean_waveforms', 'overlap_matrix', 'overlap_summary', ...
    'pc_feature_ind', 'pc_features', 'similar_templates', ...
    'spike_clusters', 'spike_templates', 'spike_times', ...
    'template_feature_ind', 'template_features', 'templates', ...
    'templates_ind', 'templates_unw', 'whitening_mat', ...
    'whitening_mat_inv'};

data = [];


% get the filename for the requested type
if any(strcmp(npy_types, type))
    fname = fullfile(data_dir, sprintf('%s.npy', type));
elseif strcmp('metrics', type)
    fname = fullfile(data_dir, 'csv', 'metrics.csv');
elseif strcmp('cluster_groups', type)
    fname = fullfile(data_dir, 'cluster_groups.csv');
elseif strcmp('ks_label', type)
    fname = fullfile(data_dir, 'cluster_KSLabel.tsv');
elseif strcmp('params', type)
    fname = fullfile(data_dir, 'params.py');
elseif strcmp('trigger', type)
    fname = fullfile(data_dir, 'trigger.mat');
elseif strcmp('rez', type)
    fname = fullfile(data_dir, 'rez.mat');
elseif strcmp('rez2', type)
    fname = fullfile(data_dir, 'rez2.mat');
elseif strcmp('selected_clusters', type)
    fname = fullfile(data_dir, 'selected_clusters.txt');
end


% exit early if filename doesn't exist
if ~exist(fname, 'file')
    warning('Filename: %s doesn''t exist', fname);
    return
end


% load the data for the requested type
if any(strcmp(npy_types, type))
    % load the npy data
    data = readNPY(fname);
elseif strcmp('metrics', type)
    % load the csv data
    data = read_metrics_csv(fname);
elseif any(strcmp({'cluster_groups', 'ks_label'}, type))
    % load the csv data (same for ks_label)
    data = read_cluster_groups(fname);
elseif strcmp('params', type)
    data = read_params_py(fname);
elseif strcmp('trigger', type)
    data = load(fname);
    data = data.trigger;
elseif strcmp('rez', type)
    data = load(fname);
    data = data.rez;
elseif strcmp('rez2', type)
    data = load(fname);
    data = data.rez;
elseif strcmp('selected_clusters', type)
    data = read_cluster_id_list(fname);
end
