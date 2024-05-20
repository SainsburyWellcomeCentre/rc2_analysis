%% same as bacward_movements.m, but with cluster loop outside starts loop

close all

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
  
    cluster_fr = cell(length(starts), length(clusters));
    cluster_spike_train = cell(length(starts), length(clusters));
    for clust_i = 1 : length(clusters) 
          
        for start_i = 1: length(starts)
            whole_interval = starts(start_i) - (one_second * 3) : starts(start_i) + (one_second * 2);
            probe_time_whole = four_way_protocol_session.probe_t(whole_interval);
            
            [train, timebase] = clusters(clust_i).fr.get_spiketrain(probe_time_whole);
            
            % plot
            figure(start_i);
            hold on;
            plot(velocity(whole_interval) / 10);
            plot(clusters(clust_i).fr.get_convolution(probe_time_whole) / 20 + 4);
            plot(train * 2 + 1);
        end

    end
end


% figure(1);
% h_ax = subplot(1, 1, 1);
% hold on;
% fmt.xy_limits       = [0, 60];
% fmt.tick_space      = 20;
% fmt.line_order      = 'top';
% fmt.xlabel          = 'before';
% fmt.ylabel          = 'after';
% fmt.include_inset   = false;
% fmt.colour_by       = 'significance';
% 
% unity_plot_plot(h_ax, median_per_cluster_before, median_per_cluster_after, direction, fmt);








