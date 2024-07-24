% Generate a figure with unity plots for the dataset passive_same_luminance_muscimol
% Compare stationary VS motion before and after muscimol
% Calculate modulation index before and after muscimol 
% and plot the paired data of the clusters that were positively modulated
% before muscimol injection

% Modulation index:
% MI = (motion - stationary) / (motion + stationary)

% Passive_protocol always viz is repeated twice (in rec1 and rec2)
% Rec1 is the session before muscimol application, rec2 after muscimol

% Requires unity_plot_plot from mvelez_ms_figures repository


% Initialization
experiment_groups           = 'passive_same_luminance_muscimol'; 
trial_types                 = {'VT', 'V', 'T_Vstatic'};

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);

% Make a figure with unity plots and modulation index
figure(1);

% Count the total number of clusters across probes
n_clusters = 0;

% Loop through conditions
for trial_type_i = 1 : length(trial_types)
    % Initialize empty arrays
    motion_fr_median_rec1 = [];
    stationary_fr_median_rec1 = [];
    direction_rec1 = [];

    motion_fr_median_rec2 = [];
    stationary_fr_median_rec2 = [];
    direction_rec2 = [];
    
    % Loop across probes
    for probe_i = 1 : length(probe_ids)
        % Get clusters and trials data
        data   = ctl.load_formatted_data(probe_ids{probe_i});
        clusters  = data.VISp_clusters();
        if trial_type_i == 1
            % On the first time we observe a condition,
            % count the cluster number
            n_clusters = n_clusters + length(clusters);
        end
        trials = data.get_trials_with_trial_group_label(trial_types{trial_type_i});
        
        % Initialzie matrices
        motion_fr_rec1 = nan(length(trials), length(clusters));
        stationary_fr_rec1 = nan(length(trials), length(clusters));
        motion_fr_rec2 = nan(length(trials), length(clusters));
        stationary_fr_rec2 = nan(length(trials), length(clusters));
        
        % Loop across trials
        for trial_i = 1 : length(trials)
            % Get the aligned trial for given index
            trial  = trials{trial_i}.to_aligned;
            
            % Focus on pre-muscimol dataset
            if contains(trial.session_id, 'rec1')
                % Get the stationary and motion masks to filter the clusters' firing rate
                stationary_mask = trial.stationary_mask;
                motion_mask = trial.motion_mask;
                
                % Loop across clusters 
                for clust_i = 1 : length(clusters) 
                    % Get the convolved firing rate
                    fr = clusters(clust_i).fr.get_convolution(trial.probe_t);

                    % Filter the firing rate
                    motion_fr_rec1(trial_i, clust_i) = nanmean(fr(motion_mask));
                    stationary_fr_rec1(trial_i, clust_i) = nanmean(fr(stationary_mask));

                end 
            end
            
            % Focus on post-muscimol dataset
            if contains(trial.session_id, 'rec2') & trial.trial_id < 31
                % Get the stationary and motion masks to filter the clusters' firing rate
                stationary_mask = trial.stationary_mask;
                motion_mask = trial.motion_mask;

                % Loop across clusters
                for clust_i = 1 : length(clusters) 
                    % Get the convolved firing rate
                    fr = clusters(clust_i).fr.get_convolution(trial.probe_t);

                    % Filter the firing rate
                    motion_fr_rec2(trial_i - 10, clust_i) = nanmean(fr(motion_mask));
                    stationary_fr_rec2(trial_i - 10, clust_i) = nanmean(fr(stationary_mask));

                end 
            end
        end
        
        % Loop across clusters to calculate the mean of the filtered firing rate
        % and perform signrank test
        for clust_i = 1 : length(clusters) 
            motion_fr_median_rec1(end+1) = nanmedian(motion_fr_rec1(:, clust_i));
            stationary_fr_median_rec1(end+1) = nanmedian(stationary_fr_rec1(:, clust_i));
            [~, ~, ~, direction_rec1(end+1)] = compare_groups_with_signrank(stationary_fr_rec1(:, clust_i), motion_fr_rec1(:, clust_i));
            
            motion_fr_median_rec2(end+1) = nanmedian(motion_fr_rec2(:, clust_i));
            stationary_fr_median_rec2(end+1) = nanmedian(stationary_fr_rec2(:, clust_i));
            [~, ~, ~, direction_rec2(end+1)] = compare_groups_with_signrank(stationary_fr_rec2(:, clust_i), motion_fr_rec2(:, clust_i));            
        end
    end
   
    % Make unity plot of stationary VS motion before muscimol
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

    % Make unity plot of stationary VS motion after muscimol
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


    % Calculate modulation index (MI) from significantly modulated cells before and after muscimol
    % And plot it with a paired plot
    % On the left pre-muscimol, on the right post-muscimol
    % Plot only clusters that were modulated before muscimol

    h_ax = subplot(4, 3, trial_type_i + 6);

    hold on;
    modulation_index_rec1 = [];
    modulation_index_rec2 = [];

    % For each cluster
    for clust_i = 1 : n_clusters
        % Calculate MI
        modulation_index_rec1(end+1) = (motion_fr_median_rec1(clust_i) - stationary_fr_median_rec1(clust_i))...
            / (motion_fr_median_rec1(clust_i) + stationary_fr_median_rec1(clust_i));

        modulation_index_rec2(end+1) = (motion_fr_median_rec2(clust_i) - stationary_fr_median_rec2(clust_i))...
            / (motion_fr_median_rec2(clust_i) + stationary_fr_median_rec2(clust_i));

        % Start plotting
        
        % Filter for clusters that were modulated only before muscimol
        if direction_rec1(clust_i) ~= 0
            % Left side: pre-muscimol
            % Make a circle for each cluster. 
            % Red: excited
            % Blue: inhibited
            if direction_rec1(clust_i) == 1
                scatter(1, modulation_index_rec1(clust_i), scatterball_size(3), 'red', 'o');  
            elseif direction_rec1(clust_i) == -1
                scatter(1, modulation_index_rec1(clust_i), scatterball_size(3), 'blue', 'o');
            end

            % Right side: post-muscimol
            % Make a circle for each cluster. 
            % Red: excited
            % Blue: inhibited
            % Black: no significant modulation 
            if direction_rec2(clust_i) == 1
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'red', 'o');
            elseif direction_rec2(clust_i) == -1
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'blue', 'o');
            else 
                scatter(2, modulation_index_rec2(clust_i), scatterball_size(3), 'black', 'o');
            end

            % Plot lines connecting the circles
            plot([1 2], [modulation_index_rec1(clust_i), modulation_index_rec2(clust_i)], 'black');
        end
    end
    xlim([0 3]);
    ylim([-1.2 1.2]);
        

    % Calculate some statics on the MI of modulated clusters and print them
    only_responsive_rec1 =  direction_rec1 ~= 0;
    n = sum(only_responsive_rec1(~isnan(only_responsive_rec1)));
    
    avg_mi_rec1 = nanmean(modulation_index_rec1(only_responsive_rec1));
    std_mi_rec1 = nanstd(modulation_index_rec1(only_responsive_rec1));
    avg_mi_rec2  = nanmean(modulation_index_rec2(only_responsive_rec1));
    std_mi_rec2  = nanstd(modulation_index_rec2(only_responsive_rec1));
    [p] = signrank(modulation_index_rec1(only_responsive_rec1), modulation_index_rec2(only_responsive_rec1));
    
    sprintf('(Responsive) Trial type: %s, before: %.2f + %.2f sem, after: %.2f + %.2f sem, p: %.5f, n: %.0f', trial_types{trial_type_i}, avg_mi_rec1, std_mi_rec1 / sqrt(n), avg_mi_rec2, std_mi_rec2 / sqrt(n), p, n)
    
    % MI from significantly and positively modulated cells
    % Same MI calculation and same plotting strategy as above
    h_ax = subplot(4, 3, trial_type_i + 9);
    hold on;
    modulation_index_rec1 = [];
    modulation_index_rec2 = [];

    for clust_i = 1 : n_clusters
        modulation_index_rec1(end+1) = (motion_fr_median_rec1(clust_i) - stationary_fr_median_rec1(clust_i))...
            / (motion_fr_median_rec1(clust_i) + stationary_fr_median_rec1(clust_i));

        modulation_index_rec2(end+1) = (motion_fr_median_rec2(clust_i) - stationary_fr_median_rec2(clust_i))...
            / (motion_fr_median_rec2(clust_i) + stationary_fr_median_rec2(clust_i));

        % Filter for only positive modulated clusters
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
        
    % Filter for only positive modulated clusters
    only_responsive_rec1 =  direction_rec1 > 0;
    n = sum(only_responsive_rec1(~isnan(only_responsive_rec1)));
    
    avg_mi_rec1 = nanmean(modulation_index_rec1(only_responsive_rec1));
    std_mi_rec1 = nanstd(modulation_index_rec1(only_responsive_rec1));
    avg_mi_rec2  = nanmean(modulation_index_rec2(only_responsive_rec1));
    std_mi_rec2  = nanstd(modulation_index_rec2(only_responsive_rec1));
    [p] = signrank(modulation_index_rec1(only_responsive_rec1), modulation_index_rec2(only_responsive_rec1));
    
    sprintf('(Responsive and positive) Trial type: %s, before: %.2f + %.2f sem, after: %.2f + %.2f sem, p: %.5f, n: %.0f', trial_types{trial_type_i}, avg_mi_rec1, std_mi_rec1 / sqrt(n), avg_mi_rec2, std_mi_rec2 / sqrt(n), p, n)
end

