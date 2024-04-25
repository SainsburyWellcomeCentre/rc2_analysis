close all;

experiment_groups           = 'passive_same_luminance';
trial_types                 = {'VT', 'V', 'T_Vstatic'};

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);

figure(1);

for trial_type_i = 1 : length(trial_types)
    motion_fr_median_rec1 = [];
    stationary_fr_median_rec1 = [];
    direction_rec1 = [];
    
    for probe_i = 1 : length(probe_ids)
        data   = ctl.load_formatted_data(probe_ids{probe_i});
        clusters  = data.VISp_clusters();
        trials = data.get_trials_with_trial_group_label(trial_types{trial_type_i});
        
        motion_fr_rec1 = nan(length(trials), length(clusters));
        stationary_fr_rec1 = nan(length(trials), length(clusters));
        
        for trial_i = 1 : length(trials)
            
            trial  = trials{trial_i}.to_aligned;
                
            if contains(trial.session_id, 'rec1')
                
                stationary_mask = trial.stationary_mask;
                motion_mask = trial.motion_mask;

                for clust_i = 1 : length(clusters) 

                    fr = clusters(clust_i).fr.get_convolution(trial.probe_t);

                    motion_fr_rec1(trial_i, clust_i) = nanmean(fr(motion_mask));
                    stationary_fr_rec1(trial_i, clust_i) = nanmean(fr(stationary_mask));

                end 
            end
        end
        
        for clust_i = 1 : length(clusters) 
            motion_fr_median_rec1(end+1) = nanmedian(motion_fr_rec1(:, clust_i));
            stationary_fr_median_rec1(end+1) = nanmedian(stationary_fr_rec1(:, clust_i));
            [~, ~, ~, direction_rec1(end+1)] = compare_groups_with_signrank(stationary_fr_rec1(:, clust_i), motion_fr_rec1(:, clust_i));
        end
    end
   
    h_ax = subplot(1, 3, trial_type_i);
    hold on;

    fmt.xy_limits       = [0, 60];
    fmt.tick_space      = 20;
    fmt.line_order      = 'top';
    fmt.xlabel          = 'FR baseline';
    fmt.ylabel          = trial_types{trial_type_i};
    fmt.include_inset   = false;
    fmt.colour_by       = 'significance';

    unity_plot_plot(h_ax, stationary_fr_median_rec1, motion_fr_median_rec1, direction_rec1, fmt);
end

