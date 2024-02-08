experiment_groups       = {'visual_flow'};
backwards_treshold = - 10;
min_dur = 0.2;

%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

median_per_cluster_before = [];
median_per_cluster_after = [];
direction = [];

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    
    % Create the windows
    four_way_protocol_session = data.sessions{1};
    
    velocity = filter_trace(four_way_protocol_session.stage);
    backwards_mask = velocity < backwards_treshold;

    starts = find(diff(backwards_mask) == 1);
    ends = find(diff(backwards_mask) == -1);

    % Add two second before and two after for comparison
    one_second = four_way_protocol_session.fs;
    
    clusters    = data.VISp_clusters();

    mean_fr_before = zeros(length(starts), length(clusters));
    mean_fr_after = zeros(length(starts), length(clusters));
    for start_i = 1: length(starts)
        before_interval = starts(start_i) - (one_second * 3) : starts(start_i) - one_second;
        after_interval = starts(start_i): starts(start_i) + (one_second * 2);
        probe_time_before = four_way_protocol_session.probe_t(before_interval);
        probe_time_after = four_way_protocol_session.probe_t(after_interval);
        
        whole_interval = starts(start_i) - (one_second * 3) : starts(start_i) + (one_second * 2);
        probe_time_whole = four_way_protocol_session.probe_t(whole_interval);
        
        cluster_fr = cell(1, length(clusters));
        for clust_i = 1 : length(clusters) 
            fr_before = clusters(clust_i).fr.get_convolution(probe_time_before);
            fr_after = clusters(clust_i).fr.get_convolution(probe_time_after);
            
            cluster_fr{clust_i} = clusters(clust_i).fr.get_convolution(probe_time_whole);
            
            % take the mean in one backwards event for the window before
            % and after
            mean_fr_before(start_i, clust_i) = mean(fr_before);
            mean_fr_after(start_i, clust_i) = mean(fr_after);
        end
        
        % plot
        population_fr = mean(cat(2, cluster_fr{:}), 2);
        figure(start_i);
        hold on;
        plot(velocity(whole_interval));
        plot(population_fr);
       
    end
    
    for clust_i = 1 : length(clusters) 
         median_per_cluster_before(end+1) = median(mean_fr_before(:, clust_i));
         median_per_cluster_after(end+1) = median(mean_fr_after(:, clust_i));
         [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(mean_fr_before(:, clust_i), mean_fr_after(:, clust_i));
    end
end


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








