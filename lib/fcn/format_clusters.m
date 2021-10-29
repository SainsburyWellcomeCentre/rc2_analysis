function clusters = format_clusters(ks_dir)
%%FORMAT_CLUSTERS
%   CLUSTERS = FORMAT_CLUSTERS(KS_DIR)
%
%   Takes the files in the output of kilosort2/ecephys_spike_sorting or
%   related fork, and produces a sequence of "cluster" objects with all the
%   associated information.
%
%   Input:  ks_dir - directory of kilosort2 output
%   Output: clusters - structure array of clusters with several properties

% read the kilosort data
params            = read_ks_data(ks_dir, 'params');
spike_clusters    = read_ks_data(ks_dir, 'spike_clusters');
spike_templates   = read_ks_data(ks_dir, 'spike_templates');
spike_times       = read_ks_data(ks_dir, 'spike_times');
amplitudes        = read_ks_data(ks_dir, 'amplitudes');
templates         = read_ks_data(ks_dir, 'templates');
chan_map          = read_ks_data(ks_dir, 'channel_map');
chan_pos          = read_ks_data(ks_dir, 'channel_positions');
qm_table          = read_ks_data(ks_dir, 'metrics');
cluster_groups    = read_ks_data(ks_dir, 'cluster_groups');
ks_label          = read_ks_data(ks_dir, 'ks_label');

%
probe_track       = read_track_csv(ks_dir);

% correct probe track
if ~isempty(probe_track)
    track_offset      = read_track_offset(ks_dir);
    probe_track.pos   = probe_track.pos + track_offset;
end


% get fs from params.py
fs = params.sample_rate;

% which quality metrics to unpack
qm = {'firing_rate', 'presence_ratio', 'isi_viol', 'amplitude_cutoff', ...
    'isolation_distance', 'l_ratio', 'd_prime', 'nn_hit_rate', ...
    'nn_miss_rate', 'silhouette_score', 'max_drift', 'cumulative_drift', ...
    'peak_channel', 'snr', 'duration', 'halfwidth', 'PT_ratio', ...
    'repolarization_slope', 'recovery_slope', 'amplitude', 'spread', ...
    'velocity_above', 'velocity_below'};

% all clusters with spikes
cluster_ids = sort(unique(spike_clusters));

% iterate over clusters
for i = length(cluster_ids) : -1 : 1
    
    % location in the total spike vector of this cluster
    spike_mask = spike_clusters == cluster_ids(i);
    
    % cluster ID
    clusters(i).id = cluster_ids(i);
    
    % times of spikes in probe
    clusters(i).spike_times = double(spike_times(spike_mask)) / fs;
    
    % sample points of spikes on probe
    clusters(i).spike_sample_point = spike_times(spike_mask);
    
    % amplitude of each spike
    clusters(i).amplitudes = amplitudes(spike_mask);
    
    % most numerous template before merging and splitting
    %   required to get template and KS properties which are not modified
    %   by phy... more for merged units (as split units will all have the
    %   same cluster_id beforehand.
    ori_cluster_id = mode(spike_templates(spike_mask));
    
    % get template of this cluster
    %   TODO: should we unwhiten the templates?
    template = permute(templates(ori_cluster_id+1, :, :), [2, 3, 1]); % +1 because id is 0 indexed
    
    % get amplitude of waveform on each channel
    amp = max(template, [], 1) - min(template, [], 1);
    
    % get channel of max amplitude
    [~, max_chan_good] = max(amp);
    
    % index correctly
    max_chan = chan_map(max_chan_good);
    y_dist = chan_pos(max_chan_good, 2);
    
    % store best channel and distance from the probe tip
    clusters(i).best_channel = max_chan;
    clusters(i).distance_from_probe_tip = 200 + y_dist; % TODO: load from settings
    
    % find position in the probe track file
    if ~isempty(probe_track)
        
        [~, idx] = min(abs(clusters(i).distance_from_probe_tip - probe_track.pos));
        
        clusters(i).region_id = probe_track.region_id(idx);
        clusters(i).region_str = probe_track.region_str{idx};
        
        % inverse position
        clusters(i).depth = probe_track.pos(1) - clusters(i).distance_from_probe_tip;
    else
        clusters(i).region_id = nan;
        clusters(i).region_str = "";
        clusters(i).depth = nan;
    end
    
    % cluster groups info
    %   this is changed by Phy so we can use the current cluster
    if ~isempty(cluster_groups)
        cluster_groups_idx = str2double(cluster_groups(:, 1)) == cluster_ids(i);
        clusters(i).class = cluster_groups(cluster_groups_idx, 2);
    else
        clusters(i).class = "";
    end
    
    % KS2 metrics
    if ~isempty(ks_label)
        % use the original template as this is NOT changed after phy
        %   takes the original cluster with most spikes contributing merged
        %   cluster
        %   if split, the original cluster will be the same label for all new
        %   clusters.
        ks_label_idx = str2double(ks_label(:, 1)) == ori_cluster_id;
        clusters(i).ks2_class = ks_label(ks_label_idx, 2);
    else
        clusters(i).ks2_class = nan;
    end
    
    % if quality metrics were read successfully
    for j = 1 : length(qm)    
        if ~isempty(qm_table)
            % metric index
            metric_idx = qm_table.cluster_id == cluster_ids(i);
            clusters(i).(qm{j}) = qm_table.(qm{j})(metric_idx);
        else
            clusters(i).(qm{j}) = nan;
        end
    end
    
    % additional metrics
    median_amp = median(clusters(i).amplitudes);
    min_amp = min(clusters(i).amplitudes);
    clusters(i).amplitude_ratio = median_amp/min_amp;
    
end
