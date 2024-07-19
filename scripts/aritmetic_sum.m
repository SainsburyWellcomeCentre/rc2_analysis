% Fit a linear model on the passive same luminance dataset,
% on the equation : VT = b_0 + b_1 * T_Vstatic + b_2 * V

% How are fold changes calculated?
% M is the 3D matrix of motion responses, with shape trial_types x clusters x trials
% S is the 3D matrix of stationary responses, with shape n_trial_types x n_clusters x n_trials
%
% Example: mean(M, trials) means we are taking the mean of M across trials, obtaining a 
% trial_types x clusters 2D matrix
%
% fold change per cluster: FC(M, S, cluster, trial_type) = median(M, trials) / mean(median(S, trials))
% fold change per trial_type: FC(M, S, trial_type) = mean(median(M, trials) / mean(median(S, trials)), clusters)

close all;

% initialization
experiment_groups           = 'passive_same_luminance';
trial_types                 = {'T_Vstatic', 'V', 'VT'};
restricted                  = true; % If restricted true, select only for clusters that are positively or negatively modulated

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);

% load sequence of stimuli
fname = 'Y:\mvelez\mateoData_rc2\passive_protocol_sequence.mat';
passive_protocol_sequence = load(fname);

% initialize variables
all_tables = cell(1,length(probe_ids));
cluster_count = 0;

for probe_i = 1 : length(probe_ids)
    % load data
    data   = ctl.load_formatted_data(probe_ids{probe_i});
    clusters  = data.VISp_clusters();
    
    clusters_id = arrayfun(@(x) x.id, clusters)';
    
    % get stationary vs motion table
    T_data = data.svm_table;
    
    % filter clusters: retain only those of witch there is a significant difference 
    % between stationary and motion.
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

    % apply filter for VISp and in case for restricted
    T_subset=T_data(ismember(T_data.cluster_id, clusters_id),:);
    
    % join tables
    all_tables{probe_i} = T_subset;
    cluster_count = cluster_count + length(clusters_id);
end

cluster_count

% concatenate the tables
T = vertcat(all_tables{:});

% Calculations for the fold change
% get the median stationary and motion response for each combination of animal / cluster / trial type
% the median is calculated across trials
median_responses = groupsummary(T, {'probe_id', 'cluster_id', 'trial_group_label'}, 'median', {'stationary_fr', 'motion_fr'});

% get the mean of the median stationary for each trial group label
mean_stationary = groupsummary(median_responses, 'trial_group_label', 'mean', 'median_stationary_fr');

len = size(median_responses);
median_responses.normalised_stationary = zeros(len(1), 1);
median_responses.normalised_motion = zeros(len(1), 1);
for trial_i = 1 : length(trial_types)
    % find indices in table of a given trial type
    idxs = find(strcmp(median_responses{:, 'trial_group_label'}, trial_types(trial_i)));

    % normalize to the mean stationary : get the fold change for each cluster
    median_responses{idxs, 'normalised_stationary'} = median_responses{idxs, 'median_stationary_fr'}./mean_stationary{trial_i, 'mean_median_stationary_fr'};
    median_responses{idxs, 'normalised_motion'} = median_responses{idxs, 'median_motion_fr'}./mean_stationary{trial_i, 'mean_median_stationary_fr'};
end

median_responses.median_stationary_fr = [];
median_responses.median_motion_fr = [];
median_responses.normalised_stationary = [];

% fit the clusters fold changes with a linear model
T_pivot = unstack(median_responses, 'normalised_motion', 'trial_group_label');
lreg = fitlm(T_pivot,'VT~V+T_Vstatic')
plot(lreg, 'Marker', 'o')


xlim ([-1 12]);
ylim([-1 12]);
title('Arithmetic sum')
xlabel('fold change of VS+T  +  VF')
ylabel('fold change of VF+T')
set(gca,'box','off')



