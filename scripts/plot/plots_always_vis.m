% Generate a figure with unity plots for the dataset passive_same_luminance
% Compare stationary VS motion

% Requires unity_plot_plot from mvelez_ms_figures repository

% Initialization
experiment_groups           = 'passive_same_luminance';
trial_types                 = {'VT', 'V', 'T_Vstatic'};

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);

% Make a figure with unity plots
figure(1);

% Loop across conditions
for trial_type_i = 1 : length(trial_types)
    motion_fr_median_rec1 = [];
    stationary_fr_median_rec1 = [];
    direction_rec1 = [];
    
    % Loop across probes
    for probe_i = 1 : length(probe_ids)
        % Get clusters and trials data
        data   = ctl.load_formatted_data(probe_ids{probe_i});
        clusters  = data.VISp_clusters();
        trials = data.get_trials_with_trial_group_label(trial_types{trial_type_i});
        
        % Initialize empty matrices
        motion_fr_rec1 = nan(length(trials), length(clusters));
        stationary_fr_rec1 = nan(length(trials), length(clusters));
        
        % Loop across trials
        for trial_i = 1 : length(trials)
            % Get the aligned trial for given index
            trial  = trials{trial_i}.to_aligned;

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
        
        % Loop across clusters to calculate the mean of the filtered firing rate
        % and perform signrank test
        for clust_i = 1 : length(clusters) 
            motion_fr_median_rec1(end+1) = nanmedian(motion_fr_rec1(:, clust_i));
            stationary_fr_median_rec1(end+1) = nanmedian(stationary_fr_rec1(:, clust_i));
            [~, ~, ~, direction_rec1(end+1)] = compare_groups_with_signrank(stationary_fr_rec1(:, clust_i), motion_fr_rec1(:, clust_i));
        end
    end
   
    % Make a unity plot
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

