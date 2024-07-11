close all;

experiment_groups           = 'passive_same_luminance_muscimol'; 
trial_types                 = {'VT', 'V', 'T_Vstatic'};

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);


figure(1);

for trial_type_i = 1 : length(trial_types)
    motion_fr_median_rec1 = [];
    stationary_fr_median_rec1 = [];
    direction_rec1 = [];

    motion_fr_median_rec2 = [];
    stationary_fr_median_rec2 = [];
    direction_rec2 = [];
    
    for probe_i = 1 : length(probe_ids)
        data   = ctl.load_formatted_data(probe_ids{probe_i});
        clusters  = data.VISp_clusters();
        trials = data.get_trials_with_trial_group_label(trial_types{trial_type_i});
        
        motion_fr_rec1 = nan(length(trials), length(clusters));
        stationary_fr_rec1 = nan(length(trials), length(clusters));
        motion_fr_rec2 = nan(length(trials), length(clusters));
        stationary_fr_rec2 = nan(length(trials), length(clusters));
        
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
            
            if contains(trial.session_id, 'rec2') & trial.trial_id < 31

                stationary_mask = trial.stationary_mask;
                motion_mask = trial.motion_mask;

                for clust_i = 1 : length(clusters) 

                    fr = clusters(clust_i).fr.get_convolution(trial.probe_t);

                    motion_fr_rec2(trial_i - 10, clust_i) = nanmean(fr(motion_mask));
                    stationary_fr_rec2(trial_i - 10, clust_i) = nanmean(fr(stationary_mask));

                end 
            end
        end
        
        for clust_i = 1 : length(clusters) 
            motion_fr_median_rec1(end+1) = nanmedian(motion_fr_rec1(:, clust_i));
            stationary_fr_median_rec1(end+1) = nanmedian(stationary_fr_rec1(:, clust_i));
            [~, ~, ~, direction_rec1(end+1)] = compare_groups_with_signrank(stationary_fr_rec1(:, clust_i), motion_fr_rec1(:, clust_i));
            
            motion_fr_median_rec2(end+1) = nanmedian(motion_fr_rec2(:, clust_i));
            stationary_fr_median_rec2(end+1) = nanmedian(stationary_fr_rec2(:, clust_i));
            [~, ~, ~, direction_rec2(end+1)] = compare_groups_with_signrank(stationary_fr_rec2(:, clust_i), motion_fr_rec2(:, clust_i));
            
        end
    end
   
    h_ax = subplot(4, 3, trial_type_i);
    hold on;

    fmt.xy_limits       = [0, 60];
    fmt.tick_space      = 20;
    fmt.line_order      = 'top';
    fmt.xlabel          = 'FR baseline';
    fmt.ylabel          = trial_types{trial_type_i};
    fmt.include_inset   = false;
    fmt.colour_by       = 'significance';

    unity_plot_plot(h_ax, stationary_fr_median_rec1, motion_fr_median_rec1, direction_rec1, fmt);

    h_ax = subplot(4, 3, trial_type_i + 3);
    hold on;

    fmt.xy_limits       = [0, 60];
    fmt.tick_space      = 20;
    fmt.line_order      = 'top';
    fmt.xlabel          = 'FR baseline';
    fmt.ylabel          = trial_types{trial_type_i};
    fmt.include_inset   = false;
    fmt.colour_by       = 'significance';

    unity_plot_plot(h_ax, stationary_fr_median_rec2, motion_fr_median_rec2, direction_rec2, fmt);


    % MI from significantly modulated cells
    h_ax = subplot(4, 3, trial_type_i + 6);
    hold on;
    modulation_index_rec1 = [];
    modulation_index_rec2 = [];

    for clust_i = 1 : 82
        modulation_index_rec1(end+1) = (motion_fr_median_rec1(clust_i) - stationary_fr_median_rec1(clust_i))...
            / (motion_fr_median_rec1(clust_i) + stationary_fr_median_rec1(clust_i));

        modulation_index_rec2(end+1) = (motion_fr_median_rec2(clust_i) - stationary_fr_median_rec2(clust_i))...
            / (motion_fr_median_rec2(clust_i) + stationary_fr_median_rec2(clust_i));

        if direction_rec1(clust_i) ~= 0
            if direction_rec2(clust_i) == 1
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'red', 'o');
            elseif direction_rec2(clust_i) == -1
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'blue', 'o');
            else 
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'black', 'o');
            end

            if direction_rec1(clust_i) == 1
                scatter(1, modulation_index_rec1(clust_i), scatterball_size(3), 'red', 'o');  
            elseif direction_rec1(clust_i) == -1
                scatter(1, modulation_index_rec1(clust_i), scatterball_size(3), 'blue', 'o');
            end
            plot([1 2], [modulation_index_rec1(clust_i), modulation_index_rec2(clust_i)], 'black');
        end
    end
    xlim([0 3]);
    ylim([-1.2 1.2]);
        

    only_responsive_rec1 = direction_rec1 ~= 0;
    avg_mi_rec1 = nanmean(modulation_index_rec1(only_responsive_rec1));
    std_mi_rec1 = nanstd(modulation_index_rec1(only_responsive_rec1));
    avg_mi_rec2  = nanmean(modulation_index_rec2(only_responsive_rec1));
    std_mi_rec2  = nanstd(modulation_index_rec2(only_responsive_rec1));
    [p] = signrank(modulation_index_rec1(only_responsive_rec1), modulation_index_rec2(only_responsive_rec1));
    
    sprintf('Trial type: %s, before: %.2f + %.2f, after: %.2f + %.2f, p: %.5f', trial_types{trial_type_i}, avg_mi_rec1, std_mi_rec1, avg_mi_rec2, std_mi_rec2, p)
    
    % MI from significantly and positively modulated cells
    h_ax = subplot(4, 3, trial_type_i + 9);
    hold on;
    modulation_index_rec1 = [];
    modulation_index_rec2 = [];

    for clust_i = 1 : 82
        modulation_index_rec1(end+1) = (motion_fr_median_rec1(clust_i) - stationary_fr_median_rec1(clust_i))...
            / (motion_fr_median_rec1(clust_i) + stationary_fr_median_rec1(clust_i));

        modulation_index_rec2(end+1) = (motion_fr_median_rec2(clust_i) - stationary_fr_median_rec2(clust_i))...
            / (motion_fr_median_rec2(clust_i) + stationary_fr_median_rec2(clust_i));

        if direction_rec1(clust_i) > 0
            if direction_rec2(clust_i) == 1
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'red', 'o');
            elseif direction_rec2(clust_i) == -1
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'blue', 'o');
            else 
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'black', 'o');
            end

            if direction_rec1(clust_i) == 1
                scatter(1, modulation_index_rec1(clust_i), scatterball_size(3), 'red', 'o');  
            end
            plot([1 2], [modulation_index_rec1(clust_i), modulation_index_rec2(clust_i)], 'black');
        end
    end
    xlim([0 3]);
    ylim([-1.2 1.2]);
        

    only_responsive_and_positive_rec1 = direction_rec1 > 0;
    avg_mi_rec1 = nanmean(modulation_index_rec1(only_responsive_and_positive_rec1));
    std_mi_rec1 = nanstd(modulation_index_rec1(only_responsive_and_positive_rec1));
    avg_mi_rec2  = nanmean(modulation_index_rec2(only_responsive_and_positive_rec1));
    std_mi_rec2  = nanstd(modulation_index_rec2(only_responsive_and_positive_rec1));
    [p] = signrank(modulation_index_rec1(only_responsive_and_positive_rec1), modulation_index_rec2(only_responsive_and_positive_rec1));
    
    sprintf('(Positive modulation) Trial type: %s, before: %.2f + %.2f, after: %.2f + %.2f, p: %.5f', trial_types{trial_type_i}, avg_mi_rec1, std_mi_rec1, avg_mi_rec2, std_mi_rec2, p)
end

