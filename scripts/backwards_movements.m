% Measure responses to backward movements
% Compare firing rate two seconds before the start and two second after
% REQUIRES adding mvelez_ms_figures to PATH

experiment_groups       = {'visual_flow'};
% backward movements criteria
backwards_treshold = - 10; % cm/s
min_dur = 0.2;

ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

median_per_cluster_before = [];
median_per_cluster_after = [];
direction = [];

for ii = 1 : length(probe_ids)
    data        = ctl.load_formatted_data(probe_ids{ii});
    
    % Search for backward movements across all session
    % cannot do this on a trial by trial basis as some backward movements
    % could be outside of the trial
    four_way_protocol_session = data.sessions{1};
    
    velocity = filter_trace(four_way_protocol_session.stage);
    backwards_mask = velocity < backwards_treshold;

    % find when the backward movements starts
    starts = find(diff(backwards_mask) == 1);

    % Add two second before and two after for comparison
    one_second = four_way_protocol_session.fs; % how many datapoints in a second
    
    % Filter by VISp clusters
    clusters    = data.VISp_clusters();

    mean_fr_before = zeros(length(starts), length(clusters));
    mean_fr_after = zeros(length(starts), length(clusters));
    for start_i = 1: length(starts)
        % calculate starting and eding times of the analyisis window
        % Before: from -3s to -1s
        % After: from 0s to 2s (during backward movement)
        before_interval = starts(start_i) - (one_second * 3) : starts(start_i) - one_second;
        after_interval = starts(start_i): starts(start_i) + (one_second * 2);
        probe_time_before = four_way_protocol_session.probe_t(before_interval); % transformation to probe time is required
        probe_time_after = four_way_protocol_session.probe_t(after_interval);
        
        % Get also the whole firing rate for comparison
        whole_interval = starts(start_i) - (one_second * 3) : starts(start_i) + (one_second * 2);
        probe_time_whole = four_way_protocol_session.probe_t(whole_interval);
        
        cluster_fr = cell(1, length(clusters));
        for clust_i = 1 : length(clusters) 
            % get firing rates in the intervals of interest
            fr_before = clusters(clust_i).fr.get_convolution(probe_time_before);
            fr_after = clusters(clust_i).fr.get_convolution(probe_time_after);
            
            % get also whole firing rate
            cluster_fr{clust_i} = clusters(clust_i).fr.get_convolution(probe_time_whole);
            
            % take the mean in one backwards event for the window before
            % and after
            mean_fr_before(start_i, clust_i) = mean(fr_before);
            mean_fr_after(start_i, clust_i) = mean(fr_after);
        end
        
        % plot a single backward event if you wish
%         population_fr = mean(cat(2, cluster_fr{:}), 2);
%         figure(start_i);
%         hold on;
%         plot(velocity(whole_interval));
%         plot(population_fr);
    end
    
    % for each cluster, get the mean firing rate before and after
    for clust_i = 1 : length(clusters) 
         median_per_cluster_before(end+1) = median(mean_fr_before(:, clust_i));
         median_per_cluster_after(end+1) = median(mean_fr_after(:, clust_i));

         % if modulated, calculate the direction (excited / suppressed)
         [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(mean_fr_before(:, clust_i), mean_fr_after(:, clust_i));
    end
end

% make unity plot
figure(1);
h_ax = subplot(1, 1, 1);
hold on;
fmt.xy_limits       = [0, 60];
fmt.tick_space      = 20;
fmt.line_order      = 'top';
fmt.xlabel          = 'before';
fmt.ylabel          = 'after';
fmt.include_inset   = false;
fmt.colour_by       = 'significance';

unity_plot_plot(h_ax, median_per_cluster_before, median_per_cluster_after, direction, fmt);








