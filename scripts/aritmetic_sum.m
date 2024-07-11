close all;

experiment_groups           = 'passive_same_luminance';
trial_types                 = {'T_Vstatic', 'V', 'VT'};
restricted                  = true;

ctl                         = RC2Analysis();
probe_ids                   = {'CAA-1121765_rec1_rec2', 'CAA-1121766_rec1'};

fname = 'Y:\mvelez\mateoData_rc2\passive_protocol_sequence.mat';
passive_protocol_sequence = load(fname);

all_tables = cell(1,length(probe_ids));
cluster_count = 0;

for probe_i = 1 : length(probe_ids)
    data   = ctl.load_formatted_data(probe_ids{probe_i});
    clusters  = data.VISp_clusters();
    
    clusters_id = arrayfun(@(x) x.id, clusters)';
    
    T_data = data.svm_table;
    
    
    if restricted
        direction_T = [];
        direction_V = [];
        for jj = 1 : length(clusters)
            [~, ~, direction_T(jj), ~, ~,] = ...
               data.is_stationary_vs_motion_significant(clusters(jj).id, trial_types{1});
            [~, ~, direction_V(jj), ~, ~,] = ...
               data.is_stationary_vs_motion_significant(clusters(jj).id, trial_types{2});
        end
        clusters_id = clusters_id(~(direction_T == 0) & ~(direction_V == 0));
    end
    
    T_subset=T_data(ismember(T_data.cluster_id, clusters_id),:);
    
    all_tables{probe_i} = T_subset;
    cluster_count = cluster_count + length(clusters_id);
end

cluster_count

T = vertcat(all_tables{:});

median_responses = groupsummary(T, {'probe_id', 'cluster_id', 'trial_group_label'}, 'median', {'stationary_fr', 'motion_fr'});
mean_stationary = groupsummary(median_responses, 'trial_group_label', 'mean', 'median_stationary_fr');
    
len = size(median_responses);
median_responses.normalised_stationary = zeros(len(1), 1);
median_responses.normalised_motion = zeros(len(1), 1);
for trial_i = 1 : length(trial_types)
    idxs = find(strcmp(median_responses{:, 'trial_group_label'}, trial_types(trial_i)));
    median_responses{idxs, 'normalised_stationary'} = median_responses{idxs, 'median_stationary_fr'}./mean_stationary{trial_i, 'mean_median_stationary_fr'};
    median_responses{idxs, 'normalised_motion'} = median_responses{idxs, 'median_motion_fr'}./mean_stationary{trial_i, 'mean_median_stationary_fr'};
end

median_responses.median_stationary_fr = [];
median_responses.median_motion_fr = [];
median_responses.normalised_stationary = [];

T_pivot = unstack(median_responses, 'normalised_motion', 'trial_group_label');
lreg = fitlm(T_pivot,'VT~V+T_Vstatic')
plot(lreg, 'Marker', 'o')


xlim ([-1 12]);
ylim([-1 12]);
title('Arithmetic sum')
xlabel('fold change of VS+T  +  VF')
ylabel('fold change of VF+T')
set(gca,'box','off')



