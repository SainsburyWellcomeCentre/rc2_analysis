%% ANALYSE FIRST PASSIVE_SAME_LUMINANCE AND VISUAL_FLOW
% Visual flow is then required to normalize the mismatch dataset
% Calculate fold change across all conditions

ctl  = RC2Analysis();

probe_ids_groups = {...
    ctl.get_probe_ids('visual_flow'), ...                       % visual_flow
    {'CAA-1121765_rec1_rec2', 'CAA-1121766_rec1'},...           % passive_same_luminance
};

% Order of the cathegories respects the order that is created by
% groupsummary (alphabetical)
trial_types_groups = {
    {'RV', 'RVT', 'V', 'VT'}, ...                               % visual_flow - edited trial types
    {'T_Vstatic', 'V', 'VT'}, ...                               % passive_same_luminance
};


for group_i = 1 : length(trial_types_groups)
    
    probe_ids = probe_ids_groups(group_i);
    trial_types = trial_types_groups(group_i);
    
    group_i
    probe_ids{1}
    
    all_tables = cell(1,length(probe_ids{1}));
    cluster_count = 0;
    for probe_i = 1 : length(probe_ids{1})
        data   = ctl.load_formatted_data(probe_ids{1}{probe_i});
        clusters  = data.VISp_clusters();
        
        clusters_id = arrayfun(@(x) x.id, clusters)';
        T_data = data.svm_table;
        T_subset = T_data(ismember(T_data.cluster_id, clusters_id),:);

        all_tables{probe_i} = T_subset;
        cluster_count = cluster_count + length(clusters);
    end
    cluster_count
    T = vertcat(all_tables{:});
    
    % substitute strings, VT and V
    if group_i == 1
        T.trial_group_label(strcmp(T.trial_group_label,'V_RV')) = {'V'};
        T.trial_group_label(strcmp(T.trial_group_label,'V_RVT')) = {'V'};
        T.trial_group_label(strcmp(T.trial_group_label,'VT_RV')) = {'VT'};
        T.trial_group_label(strcmp(T.trial_group_label,'VT_RVT')) = {'VT'};
    end

    median_responses = groupsummary(T, {'probe_id', 'cluster_id', 'trial_group_label'}, 'median', {'stationary_fr', 'motion_fr'});
    mean_stationary = groupsummary(median_responses, 'trial_group_label', 'mean', 'median_stationary_fr');
    
    % Use the mean of the stationary to normalise both the stationary and
    % motion response 
    len = size(median_responses);
    median_responses.normalised_stationary = zeros(len(1), 1);
    median_responses.normalised_motion = zeros(len(1), 1);
    for type_i = 1 : length(trial_types{1})
        idxs = find(strcmp(median_responses{:, 'trial_group_label'}, trial_types{1}(type_i)));
        median_responses{idxs, 'normalised_stationary'} = median_responses{idxs, 'median_stationary_fr'}./mean_stationary{type_i, 'mean_median_stationary_fr'};
        median_responses{idxs, 'normalised_motion'} = median_responses{idxs, 'median_motion_fr'}./mean_stationary{type_i, 'mean_median_stationary_fr'};
    end
    
    % verify mean stationary is 1
    grpstats(median_responses, "trial_group_label", "mean", "DataVars", "normalised_stationary");
    
    % print fold changes, mean and sem
    stats = grpstats(median_responses, "trial_group_label", ["mean", "sem"], "DataVars", "normalised_motion")
    
    if group_i == 1
        four_way_dataset_stats = stats;
    else
        % Normalise arithmetic sum dataset to the visual_flow (four way
        % protocol) dataset
        median_responses.T_norm2_four_way_dataset_stats = zeros(len(1), 1);

        type_idx_4w = find(strcmp(four_way_dataset_stats{:, 'trial_group_label'}, 'VT'));
        type_idx_stats = find(strcmp(stats{:, 'trial_group_label'}, 'VT'));        
        
        idx_T = find(strcmp(median_responses{:, 'trial_group_label'}, 'T_Vstatic'));
        idx_VT = find(strcmp(median_responses{:, 'trial_group_label'}, 'VT'));
        median_responses{idx_T, 'T_norm2_four_way_dataset_stats'} = median_responses{idx_T, 'normalised_motion'}.*four_way_dataset_stats{type_idx_4w, 'mean_normalised_motion'}./stats{type_idx_stats, 'mean_normalised_motion'};

        % print fold changes, mean and sem
        passive_same_luminance_norm2_four_way_dataset_stats = grpstats(median_responses, "trial_group_label", ["mean", "sem"], "DataVars", "T_norm2_four_way_dataset_stats")
    end
end


%% NOW ANALYSE THE MISMATCH DATASET

% Calculate for mismatch separately
ctl  = RC2Analysis();

probe_ids = ctl.get_probe_ids('mismatch_nov20', 'mismatch_jul21');
trial_types = {'RVT_gain_up', 'RV_gain_up'};


mm                  = MismatchAnalysis();
mm.method           = 'ranksum'; % in order to use the median...

column_names = {'probe_id', 'cluster_id', 'trial_group_label', 'avg_pre_slip', 'avg_slip'};
variable_types = {'string', 'int32', 'string', 'double', 'double'};
slip_data = table('Size', [0, numel(column_names)], 'VariableNames', column_names, 'VariableTypes', variable_types);

for ii = 1 : length(probe_ids)
    data   = ctl.load_formatted_data(probe_ids{ii});
    
    if isempty(data)
        this_data   = ctl.load_formatted_data(probe_ids{ii});
    else
        this_data   = get_data_for_probe_id(data, probe_ids{ii});
    end
    
    clusters        = this_data.VISp_clusters();
    
    for trial_i = 1 : length(trial_types)
        trials          = this_data.get_trials_with_trial_group_label(trial_types(trial_i));

        % remove trials with small mismatch period
        trials          = remove_invalid_mm_trials(trials);

        for jj = 1 : length(clusters)
            slip_data(end+1, :) = {probe_ids{ii}, clusters(jj).id, trial_types(trial_i), ...
                                             mm.get_avg_baseline_fr(clusters(jj), trials), ...
                                             mm.get_avg_response_fr(clusters(jj), trials)};

        end
    end
end

all_tables = cell(1,length(probe_ids{1}));
cluster_count = 0;
for probe_i = 1 : length(probe_ids)
    data   = ctl.load_formatted_data(probe_ids{probe_i});
    clusters  = data.VISp_clusters();

    clusters_id = arrayfun(@(x) x.id, clusters)';
    T_data = data.svm_table;
    T_subset = T_data(ismember(T_data.cluster_id, clusters_id),:);

    all_tables{probe_i} = T_subset;
    cluster_count = cluster_count + length(clusters);
end
cluster_count


T = vertcat(all_tables{:});
responses = groupsummary(T, {'probe_id', 'cluster_id', 'trial_group_label'}, 'median', {'stationary_fr'});
responses.probe_id = cellfun(@string, responses.probe_id);
responses.trial_group_label = cellfun(@string, responses.trial_group_label);

% keep only RV and RVT trials
idxs_RVT = strcmp(responses{:, 'trial_group_label'}, trial_types(1));
idxs_RV = strcmp(responses{:, 'trial_group_label'}, trial_types(2));
responses = responses(idxs_RVT | idxs_RV, :);

% include slip averages
responses = join(responses, slip_data);

% Use the mean of the stationary to normalise both the stationary and
% motion responses 
mean_stationary = groupsummary(responses, 'trial_group_label', 'mean', 'median_stationary_fr');

len = size(responses);
responses.normalised_stationary = zeros(len(1), 1);
responses.normalised_avg_pre_slip = zeros(len(1), 1);
responses.normalised_avg_slip = zeros(len(1), 1);
for type_i = 1 : length(trial_types)
    idxs = find(strcmp(responses{:, 'trial_group_label'}, trial_types(type_i)));
    responses{idxs, 'normalised_stationary'} = responses{idxs, 'median_stationary_fr'}./mean_stationary{type_i, 'mean_median_stationary_fr'};
    responses{idxs, 'normalised_avg_pre_slip'} = responses{idxs, 'avg_pre_slip'}./mean_stationary{type_i, 'mean_median_stationary_fr'};
    responses{idxs, 'normalised_avg_slip'} = responses{idxs, 'avg_slip'}./mean_stationary{type_i, 'mean_median_stationary_fr'};
end

% verify mean stationary is 1
grpstats(responses, "trial_group_label", "mean", "DataVars", "normalised_stationary")

% print fold changes, mean and sem
slip_stats = grpstats(responses, "trial_group_label", ["mean", "sem"], "DataVars", ["normalised_avg_pre_slip", "normalised_avg_slip"])



% now normalise these values based on the visual_flowdataset
responses.avg_pre_slip_norm2_four_way_dataset_stats = zeros(len(1), 1);
responses.avg_slip_norm2_four_way_dataset_stats = zeros(len(1), 1);
matched_type = {'RVT', 'RV'};
for type_i = 1 : length(trial_types)
    type_idx = find(strcmp(four_way_dataset_stats{:, 'trial_group_label'}, matched_type(type_i)));
    idxs = find(strcmp(responses{:, 'trial_group_label'}, trial_types(type_i)));
    responses{idxs, 'avg_pre_slip_norm2_four_way_dataset_stats'} = responses{idxs, 'normalised_avg_pre_slip'}.*four_way_dataset_stats{type_idx, 'mean_normalised_motion'}./slip_stats{type_i, 'mean_normalised_avg_pre_slip'};
    responses{idxs, 'avg_slip_norm2_four_way_dataset_stats'} = responses{idxs, 'normalised_avg_slip'}.*four_way_dataset_stats{type_idx, 'mean_normalised_motion'}./slip_stats{type_i, 'mean_normalised_avg_pre_slip'};
end


% print fold changes, mean and sem
slip_stats_norm2_four_way_dataset_stats = grpstats(responses, "trial_group_label", ["mean", "sem"], "DataVars", ["avg_pre_slip_norm2_four_way_dataset_stats", "avg_slip_norm2_four_way_dataset_stats"])







