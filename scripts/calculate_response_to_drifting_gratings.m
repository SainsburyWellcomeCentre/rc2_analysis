close all
clear all

sequence = load('C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sf_tf\sf_tf_rc2_20200706.mat');

ctl = RC2Analysis();
probe_ids  = ctl.get_probe_ids('visual_flow', ...
                               'mismatch_nov20', ...
                               'mismatch_jul21', ...
                               'mismatch_darkness_oct21', ...
                               'darkness');
mean_clusters_responses = [];
mean_clusters_baselines = [];
median_clusters_responses = [];
median_clusters_baselines = [];
direction = [];

columnNames = {'sf', 'tf', 'dir'};
preferred_sf_tf_dir = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'double'}, 'VariableNames', columnNames);

for probe_i = 1 : length(probe_ids)
    
    this_probe = probe_ids(probe_i)
    if strcmp(this_probe, 'CAA-1112872_rec1_rec1b_rec2_rec3') || strcmp(this_probe, 'CAA-1112874_rec1_rec2_rec3') ||  ...
       strcmp(this_probe, 'CAA-1114977_rec1_rec2_rec3') || strcmp(this_probe, 'CAA-1114978_rec1_rec2') || ...
        strcmp(this_probe, 'CAA-1114979_rec1_rec2') || strcmp(this_probe, 'CAA-1114980_rec1_rec2') || ...
        strcmp(this_probe, 'CAA-1115689_rec1_rec2') || strcmp(this_probe, 'CAA-1115691_rec1_rec2') || ...
        strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3') || strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3') || ...
        strcmp(this_probe, 'CAA-1110264_rec1_rec2')
        % skip this probes as the photodiode signal is corrupted or no sf_tf
        % session present
        continue
    end

    data = ctl.load_formatted_data(this_probe{1});
    drifting_gratings_session = data.sessions{1, 3};
    
    animal_id = split(this_probe, "_");
    if strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3') || strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3')
        animal_id = join([animal_id(1), animal_id(2), animal_id(3)], '_');
    else
        animal_id = animal_id(1);
    end
    starts = pd_times(animal_id, drifting_gratings_session, "sf_tf");
    
    
    clusters = data.VISp_clusters();
    
    for clust_i = 1 : length(clusters)
        responses = [];
        baselines = [];
       
        % Identify total spikes in the time interval of interest
        
        stim_sampled_points = 25000; % duration of one stimulus is assumed to be 2.5s
        
        % Take the first four portions of baseline
        
        for bs_i = 4 : -1 : 1
            start_ = starts(1) - stim_sampled_points * bs_i;
            end_ = starts(1) - stim_sampled_points * (bs_i - 1);
            baselines(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(start_ : end_)));
        end
        
        deltas = [];
        for i = 1 : length(starts) - 1
            delta = starts(i + 1) - starts(i);
            if delta < 50000
                deltas(end+1) = delta;
                responses(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i) : starts(i+1))));
            else
                mean_duration = cast(mean(deltas), "uint32");
                responses(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i): (starts(i) + mean_duration))));
                
                for bs_i = 1 : 4
                    start_ = starts(i) + stim_sampled_points * bs_i;
                    end_ = starts(i) + stim_sampled_points * (bs_i + 1);
                    baselines(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(start_ : end_)));
                end
            end
        end
        
        responses(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i + 1) : (starts(i + 1) + mean_duration))));

        mean_clusters_responses(end+1) = nanmean(responses);
        mean_clusters_baselines(end+1) = nanmean(baselines);
        
        
        sf = sequence.schedule.spatial_frequencies;
        tf = sequence.schedule.temporal_frequencies;
        dir = sequence.schedule.directions;
        T = horzcat(responses.', sf, tf, dir);
        T = array2table(T, "VariableNames", ["responses" "sf" "tf" "dir"]);
        
        summary = groupsummary(T, {'sf', 'tf', 'dir'}, 'median', 'responses');
        summary = sortrows(summary, 'median_responses', 'descend');
        
        % Take the top 4 stimulus combinations
        rows = (T.sf == summary{1, 1} & T.tf == summary{1, 2} & T.dir == summary{1, 3}) | ...
               (T.sf == summary{2, 1} & T.tf == summary{2, 2} & T.dir == summary{2, 3}) | ...
               (T.sf == summary{3, 1} & T.tf == summary{3, 2} & T.dir == summary{3, 3}) | ...
               (T.sf == summary{4, 1} & T.tf == summary{4, 2} & T.dir == summary{4, 3});
    
        responses_to_preferred_stim = T{rows, 'responses'};
        
        median_clusters_responses(end+1) = nanmedian(responses_to_preferred_stim);
        median_clusters_baselines(end+1) = nanmedian(baselines);
        [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(baselines, responses_to_preferred_stim);
        
        preferred_sf_tf_dir = [preferred_sf_tf_dir; summary(1, columnNames)];
   end
end
preferred_sf_tf_dir.sf = categorical(preferred_sf_tf_dir.sf);
preferred_sf_tf_dir.tf = categorical(preferred_sf_tf_dir.tf);
preferred_sf_tf_dir.dir = categorical(preferred_sf_tf_dir.dir);

% Find groups for each variable
[~, ~, group_sf] = unique(preferred_sf_tf_dir.sf);
[~, ~, group_tf] = unique(preferred_sf_tf_dir.tf);
[~, ~, group_dir] = unique(preferred_sf_tf_dir.dir);

% Apply splitapply to count occurrences in each group for each variable
count_sf = splitapply(@numel, group_sf, group_sf);
count_tf = splitapply(@numel, group_tf, group_tf);
count_dir = splitapply(@numel, group_dir, group_dir);

% Get unique values of each variable
unique_sf = unique(preferred_sf_tf_dir.sf);
unique_tf = unique(preferred_sf_tf_dir.tf);
unique_dir = unique(preferred_sf_tf_dir.dir);

% Convert count arrays to categorical
count_sf_cat = categorical(count_sf);
count_tf_cat = categorical(count_tf);
count_dir_cat = categorical(count_dir);

% Combine unique values and counts into arrays for each variable
sf_counts = [unique_sf, count_sf_cat]
tf_counts = [unique_tf, count_tf_cat]
dir_counts = [unique_dir, count_dir_cat]



mean(mean_clusters_responses)
std(mean_clusters_responses)

mean(mean_clusters_baselines)
std(mean_clusters_baselines)


% Lineplot
figure(1)
hold on
for clust_i = 1 : length(mean_clusters_responses)
    scatter(1, mean_clusters_baselines(clust_i), scatterball_size(1), 'blue', 'o');
    scatter(2, mean_clusters_responses(clust_i), scatterball_size(1), 'blue', 'o');
    plot([1 2], [mean_clusters_baselines(clust_i), mean_clusters_responses(clust_i)], 'blue', 'LineWidth', 0.1);
end
scatter(1, mean(mean_clusters_baselines), scatterball_size(3), 'black');
scatter(2, mean(mean_clusters_responses), scatterball_size(3), 'black');
plot([1 2], [mean(mean_clusters_baselines), mean(mean_clusters_responses)], 'black', 'LineWidth', 5);
xlim([0 3])
ylim([-2 100])

[p] = signrank(median_clusters_baselines, median_clusters_responses)


% Unity plot
figure(2);
h_ax = subplot(1, 1, 1);

hold on;

fmt.xy_limits       = [0, 40];
fmt.tick_space      = 20;
fmt.line_order      = 'top';
fmt.xlabel          = 'baseline';
fmt.ylabel          = 'response';
fmt.include_inset   = false;
fmt.colour_by       = 'significance';

unity_plot_plot(h_ax, median_clusters_baselines, median_clusters_responses, direction, fmt);

title('Original');





