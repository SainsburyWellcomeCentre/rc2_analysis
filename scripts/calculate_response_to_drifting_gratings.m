% Calculate responses to drifting gratings 
%
% On the stimuli and the way baslines are defined:
% Drifting gratings are characterized by a given spatial frequency (sf),
% temporal frequency (tf), and direction (dir). For every new drifting grating stimulus, starting times 
% are found via the photodiode signal.
% SF_TF stimuli are delivered in batches of 5. Every stimulus lasts 2.5s. 
% Stimuli starting times are derived by changes in photodiode signal (starts).
% If an interval is longer than 2.5s, then it means that there were no drifting gratings
% and the interval can be taken to measure baseline response. Each baseline period is long ehough to 
% be divided in four parts. This division helps us achieve enough statistical power. 
%
% How are responses to drifting gratings calculated?
% 1. Stimuli sequence is loaded
% 2. Starting times are found 
% 3. The mean firing rates are stored for baseline and response intervals
% 4. Response intervals are matched to stimulus type (which sf, tf, and dir) via the sequence information
% 5. For each cluster, the median of these intervals is calculated for baseline and stimulus_interval
% 6. The 4 stimuli which triggered the most response are used to measure responsivenes against the baseline
% 7. Information about the most preferred stimulus is printed
% 8. Plots of mean and median responses are made

stim_sampled_points = 25000; % duration of one stimulus is assumed to be 2.5s

% Load identity of sf_tf stimuli
sequence = load('C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sf_tf\sf_tf_rc2_20200706.mat');
sf = sequence.schedule.spatial_frequencies;
tf = sequence.schedule.temporal_frequencies;
dir = sequence.schedule.directions;

% Load data
ctl = RC2Analysis();
probe_ids  = ctl.get_probe_ids('visual_flow', ...
                               'mismatch_nov20', ...
                               'darkness');
                               %'mismatch_jul21', ...  % no drifting gratings sessions in these datasets
                               %'mismatch_darkness_oct21', ...
                               

% initialize arrays                              
mean_clusters_responses = [];
mean_clusters_baselines = [];
median_clusters_responses = [];
median_clusters_baselines = [];
direction = [];

% initialise table to store drifting gratings preferences
columnNames = {'sf', 'tf', 'dir'};
preferred_sf_tf_dir = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'double'}, 'VariableNames', columnNames);

for probe_i = 1 : length(probe_ids)
    
    this_probe = probe_ids(probe_i)

    if strcmp(this_probe, 'CAA-1112872_rec1_rec1b_rec2_rec3') || strcmp(this_probe, 'CAA-1112874_rec1_rec2_rec3') ||  ...
       strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3')          || strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3')    || ...
       strcmp(this_probe, 'CAA-1110264_rec1_rec2')
        % skip this probes as the photodiode signal is corrupted or no sf_tf
        % session present
        continue
    end

    % load drifting gratings session
    data = ctl.load_formatted_data(this_probe{1});
    drifting_gratings_session = data.sessions{1, 3};
    
    % Get animal ID - useful to obtain custom parameters to find starting times from the photodiode signal
    animal_id = get_animal_id(this_probe)

    % find starting times of the stimuli from photodiode
    starts = photodiode_times(animal_id, drifting_gratings_session, "sf_tf");
    
    % for each cluster, calculate response and baseline firing rate for each sf/tf/dir combination
    clusters = data.VISp_clusters();
    for clust_i = 1 : length(clusters)
        responses = [];
        baselines = [];

        % Get the firing rate in the first four baseline intervals, before the first batch of stimuli starts
        for bs_i = 4 : -1 : 1 % bo backwards from the first start to identify the baseline intervals
            start_           = starts(1) - stim_sampled_points * bs_i;
            end_             = starts(1) - stim_sampled_points * (bs_i - 1);
            baselines(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(start_ : end_)));
        end
        
        % store lengths of drifting grating stimuli
        deltas = [];
        for i = 1 : length(starts) - 1
            % calculate the length of the current stimulus
            delta = starts(i + 1) - starts(i);

            % if it's small enough, it's a drifting grating
            if delta < stim_sampled_points * 2
                % append currect delta 
                deltas(end+1) = delta;
                % store response
                responses(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i) : starts(i+1))));
            else
                % If the interval is long, we have the last drifting grating stimulus and the baseline period.
                % The threshold applied to the photodiode signal cannot pick up the end of the last drifting grating,
                % therefore the mean_duration of the stimulus is used to estimate its end.
                mean_duration = cast(mean(deltas), "uint32");
                responses(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i): (starts(i) + mean_duration))));
                
                % store the four baseline chunks.
                for bs_i = 1 : 4
                    start_ = starts(i) + stim_sampled_points * bs_i;
                    end_ = starts(i) + stim_sampled_points * (bs_i + 1);
                    baselines(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(start_ : end_)));
                end
            end
        end
        
        % get the response to the last drifting grating.
        responses(end+1) = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i + 1) : (starts(i + 1) + mean_duration))));

        % calculate mean firing rate in response and baseline intervals
        mean_clusters_responses(end+1) = nanmean(responses);
        mean_clusters_baselines(end+1) = nanmean(baselines);
        
        % add responses and stimulus information (which sf / tf / dir combination) to the table
        T = horzcat(responses.', sf, tf, dir);
        T = array2table(T, "VariableNames", ["responses" "sf" "tf" "dir"]);
        
        % calculate median responses and sort
        summary = groupsummary(T, {'sf', 'tf', 'dir'}, 'median', 'responses');
        summary = sortrows(summary, 'median_responses', 'descend');
        
        % Take the top 4 stimulus combinations
        rows = (T.sf == summary{1, 1} & T.tf == summary{1, 2} & T.dir == summary{1, 3}) | ...
               (T.sf == summary{2, 1} & T.tf == summary{2, 2} & T.dir == summary{2, 3}) | ...
               (T.sf == summary{3, 1} & T.tf == summary{3, 2} & T.dir == summary{3, 3}) | ...
               (T.sf == summary{4, 1} & T.tf == summary{4, 2} & T.dir == summary{4, 3});
        responses_to_preferred_stim = T{rows, 'responses'};
        
        % Measure responsiveness to drifting gratings. 
        % Compare the four preferred responses and compare to baseline using signrank.
        median_clusters_responses(end+1) = nanmedian(responses_to_preferred_stim);
        median_clusters_baselines(end+1) = nanmedian(baselines);
        [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(baselines, responses_to_preferred_stim);
        
        % store absolute preferred response
        preferred_sf_tf_dir = [preferred_sf_tf_dir; summary(1, columnNames)];
   end
end

% Convert preferred stim array to categorical
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


% Display statistics
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





