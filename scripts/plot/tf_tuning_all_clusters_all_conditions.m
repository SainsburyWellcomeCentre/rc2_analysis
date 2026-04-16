% Plot temporal frequency tuning curves for all cortical (VISp) clusters in all
% experiments. This script analyzes TF tuning by pooling data across all three
% TF batches (different VR gains) after transforming velocity to temporal frequency.
%
% TF transformation: TF = velocity * gain, where gain depends on batch:
%   Batch 1: 30 cm/s = 1 Hz  (gain = 1/30)
%   Batch 2: 30 cm/s = 2 Hz  (gain = 2/30)
%   Batch 3: 30 cm/s = 4 Hz  (gain = 4/30)
%
% For each cluster, creates tuning curves showing how firing rate varies with
% temporal frequency of the visual stimulus. Fits an asymmetric Gaussian function
% to the data: R = R_max * exp(-(x - x_max)^2 / sigma_x), where sigma_x differs
% for x < x_max (sigma_minus) and x >= x_max (sigma_plus). Includes shuffling 
% controls for significance testing.
%
% If save_figs is true, one PDF is created for each probe recording containing
% an A4 page for each cluster with TF tuning curves for multiple conditions.
%
% See also: speed_tuning_all_clusters_all_conditions, create_tf_tuning_tables, AsymmetricGaussianFit

%% Configuration
experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};
restricted              = false;  % If true, only include clusters tuned in both VT and V

save_figs               = true;
overwrite               = true;
figure_dir              = {'tf_tuning_curves'};

%% Collect data

% Main controller object
ctl                     = RC2Analysis();

% Get all probe_ids occurring in the specified experiment groups
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

% Setup directory to save figures
ctl.setup_figures(figure_dir, save_figs);

% Arrays and cell arrays to store information for each experiment,
% VISp cluster, and trial group type
tuning          = {};
p_svm           = [];
direction       = [];
probe_id        = {};
cluster_id      = [];
cluster_region  = {};

% ANOVA results across TF bins
anova_p         = [];  % p-values from one-way ANOVA
anova_F         = [];  % F-statistics
anova_df        = [];  % degrees of freedom (stored as [df_between, df_within])

% V+T data storage (T_Vstatic + V deviations from baseline)
tuning_VplusT   = {};

% Two-way ANOVA results between VT and V+T
anova2_VT_VplusT_p = [];  % p-value for condition effect
anova2_VT_VplusT_F = [];  % F-statistic for condition effect

n_clusters = 0;

% Loop over experiments
for ii = 1 : length(probe_ids)

    % Load data for the experiment
    data        = ctl.load_formatted_data(probe_ids{ii});

    % Return array of Cluster objects containing info about the VISp clusters
    clusters    = data.VISp_clusters();

    n_clusters = n_clusters + length(clusters);

    % Load TF tuning curves for all clusters from precalculated .mat files
    % Will be cell array of length equal to 'trial_group_labels', each
    % entry is itself a structure array with each entry the tuning curve
    % info for a cluster
    tbl = data.ctl.load_tf_tuning_curves(data.probe_id);
    
    % Extract tuning curves for each trial group
    tuning_curves = cell(1, length(trial_group_labels));
    for kk = 1:length(trial_group_labels)
        % Find matching trial group in loaded data
        group_idx = cellfun(@(x)(isequal(x, trial_group_labels{kk})), tbl.trial_groups);
        if any(group_idx)
            tuning_curves{kk} = tbl.tuning_curves{group_idx};
        else
            tuning_curves{kk} = [];
        end
    end
    
    % Loop over all VISp clusters
    for jj = 1 : length(clusters)
        
        % Loop over trial group types
        for kk = 1 : length(trial_group_labels)
            
            fprintf('Collecting for %i/%i, %i/%i, %i/%i\n', ...
                    kk, length(trial_group_labels), ...
                    jj, length(clusters), ...
                    ii, length(probe_ids));
            
            % If tuning_curves{kk} is empty there were no trials of the
            % specified type in the experiment
            if isempty(tuning_curves{kk})
                
                % Store nan's and empty array 
                p_svm(ii, jj, kk) = nan;
                direction(ii, jj, kk) = nan;
                tuning{ii}{jj}{kk} = [];
                anova_p(ii, jj, kk) = nan;
                anova_F(ii, jj, kk) = nan;
                anova_df(ii, jj, kk, 1:2) = [nan, nan];
            else
                
                % For this cluster, and this collection of trial types
                % return statistics of whether activity during stationary
                % and motion periods was significantly different
                [~, p_svm(ii, jj, kk), direction(ii, jj, kk)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{kk});
                
                % Find the tuning curve information for the current cluster
                % in 'tuning_curves' and re-store it in 'tuning'
                cluster_idx = [tuning_curves{kk}(:).cluster_id] == clusters(jj).id;
                tuning{ii}{jj}{kk} = tuning_curves{kk}(cluster_idx);
                
                % Perform one-way ANOVA across TF bins
                if ~isempty(tuning{ii}{jj}{kk}) && isfield(tuning{ii}{jj}{kk}, 'tuning')
                    % Get firing rates for each TF bin (average across repetitions)
                    tf_data = tuning{ii}{jj}{kk}.tuning;  % rows = TF bins, cols = repetitions
                    
                    % Prepare data for ANOVA: create groups and values
                    groups = [];
                    values = [];
                    for bin_idx = 1:size(tf_data, 1)
                        bin_values = tf_data(bin_idx, :);
                        bin_values = bin_values(~isnan(bin_values));  % Remove NaNs
                        groups = [groups; repmat(bin_idx, length(bin_values), 1)];
                        values = [values; bin_values(:)];
                    end
                    
                    % Run one-way ANOVA if we have data
                    if ~isempty(values) && length(unique(groups)) > 1
                        [p, tbl, ~] = anova1(values, groups, 'off');
                        anova_p(ii, jj, kk) = p;
                        anova_F(ii, jj, kk) = tbl{2, 5};  % F-statistic
                        anova_df(ii, jj, kk, 1) = tbl{2, 3};  % df between groups
                        anova_df(ii, jj, kk, 2) = tbl{3, 3};  % df within groups
                    else
                        anova_p(ii, jj, kk) = nan;
                        anova_F(ii, jj, kk) = nan;
                        anova_df(ii, jj, kk, 1:2) = [nan, nan];
                    end
                else
                    anova_p(ii, jj, kk) = nan;
                    anova_F(ii, jj, kk) = nan;
                    anova_df(ii, jj, kk, 1:2) = [nan, nan];
                end
            end
        end
        
        % Calculate V+T data: T_Vstatic + (V - V_baseline)
        % V+T combines the T_Vstatic tuning with deviations from V baseline
        if ~isempty(tuning{ii}{jj}{2}) && ~isempty(tuning{ii}{jj}{3})
            % Get V data
            V_tuning = tuning{ii}{jj}{2}.tuning;  % [n_bins x n_trials]
            V_baseline = nanmean(tuning{ii}{jj}{2}.stationary_fr);  % scalar
            
            % Get T_Vstatic data
            T_tuning = tuning{ii}{jj}{3}.tuning;  % [n_bins x n_trials]
            
            % Calculate V deviations from baseline for each bin
            V_deviations = nanmean(V_tuning, 2) - V_baseline;  % [n_bins x 1]
            
            % Add deviations to each trial of T_Vstatic
            VplusT_tuning = T_tuning + repmat(V_deviations, 1, size(T_tuning, 2));
            
            % Store V+T data in same format as other tuning data
            tuning_VplusT{ii}{jj}.tuning = VplusT_tuning;
            tuning_VplusT{ii}{jj}.cluster_id = clusters(jj).id;
            tuning_VplusT{ii}{jj}.stationary_fr = tuning{ii}{jj}{2}.stationary_fr;  % Use V baseline
            tuning_VplusT{ii}{jj}.bin_centers = tuning{ii}{jj}{2}.bin_centers;  % Copy bin_centers from V
            % Note: No shuffled data for V+T as it's a derived computation
            tuning_VplusT{ii}{jj}.shuffled = [];  % Empty - no shuffling for derived data
            
            % Run two-way ANOVA between VT and V+T
            % Rows = conditions (VT, V+T), Columns = TF bins
            VT_tuning = tuning{ii}{jj}{1}.tuning;  % [n_bins x n_trials]
            
            % Find minimum number of trials
            n_trials = min(size(VT_tuning, 2), size(VplusT_tuning, 2));
            
            % Create matrix for anova2: [condition1_trials; condition2_trials]
            % Each row is a trial, grouped by condition
            anova2_data = [VT_tuning(:, 1:n_trials)'; VplusT_tuning(:, 1:n_trials)'];
            
            % drop rows with NaNs
            anova2_data = anova2_data(~any(isnan(anova2_data), 2), :);
            
            % Update n_trials after removing NaN rows
            % After removing rows, we have equal trials from each condition
            n_trials_clean = size(anova2_data, 1) / 2;
            
            try
                [p, tbl, ~] = anova2(anova2_data, n_trials_clean, 'off');
                % p(1) = condition effect (VT vs V+T)
                % p(2) = bin effect
                % p(3) = interaction
                anova2_VT_VplusT_p(ii, jj) = p(1);
                anova2_VT_VplusT_F(ii, jj) = tbl{2, 5};  % F-statistic for condition
            catch
                anova2_VT_VplusT_p(ii, jj) = nan;
                anova2_VT_VplusT_F(ii, jj) = nan;
            end
        else
            tuning_VplusT{ii}{jj} = [];
            anova2_VT_VplusT_p(ii, jj) = nan;
            anova2_VT_VplusT_F(ii, jj) = nan;
        end
        
        % Extra info to plot on the figure
        probe_id{ii, jj} = probe_ids{ii};
        cluster_id(ii, jj) = clusters(jj).id;
        cluster_region{ii, jj} = clusters(jj).region_str;
    end
end

fprintf('\nTotal clusters processed: %d\n\n', n_clusters);

%% Plotting

% Loop over experiments
for ii = 1 : length(probe_ids)  
    
    % Loop over VISp clusters in that experiment
    for jj = 1 : length(tuning{ii})
        
        % Create A4 figure and object controlling the array of axes on the figure
        h_fig                   = ctl.figs.a4figure('landscape');
        plot_array              = PlotArray(2, 4);  % 2 rows, 4 columns for 8 subplots
        plot_array.ax_size_cm   = 5.0;  % Larger axes for better visibility
        plot_array.nx_total     = 4;  % 4 columns in underlying grid
        plot_array.ny_total     = 2;  % 2 rows in underlying grid
        
        % Cell arrays which will contain objects controlling the main
        % tuning curve plot and the shuffled tuning curve histogram info
        tuning_curve_plot       = {};
        shuff_histogram         = {};
        
        % Loop over trial group types
        for kk = 1 : length(trial_group_labels)
            
            % If there was no information in an experiment for a particular
            % trial type then skip this trial type
            if isempty(tuning{ii}{jj}{kk})
                continue
            end
            
            % Find where the axis should go and create it
            pos         = plot_array.get_position(2*kk-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            % Object controlling the tuning curve plot (uses asymmetric Gaussian fit)
            tuning_curve_plot{kk} = TFTuningCurvePlot(h_ax);
            
            % Determine the colour of the plot from the R² shuffle test results
            if tuning{ii}{jj}{kk}.shuffled.p < 0.05
                % Significant tuning from shuffle test
                main_col = [0.85, 0.32, 0.1];
            else
                % No significant tuning
                main_col = 'k';
            end
            
            % Plot the tuning curve and give axis a title with model, p-values, and ANOVA
            tuning_curve_plot{kk}.main_col = main_col;
            tuning_curve_plot{kk}.plot(tuning{ii}{jj}{kk});
            
            % Create title with condition name, model name, R-squared on mean, and ANOVA p-value
            if ~isnan(anova_p(ii, jj, kk))
                % Get model info from shuffled structure
                model_name = tuning{ii}{jj}{kk}.shuffled.best_model;
                
                % Get R-squared on mean (pre-calculated during model selection)
                if isfield(tuning{ii}{jj}{kk}.shuffled.best_model_info, 'rsq_on_mean')
                    rsq_mean = tuning{ii}{jj}{kk}.shuffled.best_model_info.rsq_on_mean;
                else
                    rsq_mean = nan;
                end
                
                % Format R-squared
                if ~isnan(rsq_mean)
                    rsq_str = sprintf('R^2=%.3f', rsq_mean);
                else
                    rsq_str = 'R^2=N/A';
                end
                
                % Format ANOVA p-value
                if anova_p(ii, jj, kk) < 0.001
                    anova_p_str = sprintf('p_{ANOVA}<0.001');
                else
                    anova_p_str = sprintf('p_{ANOVA}=%.3f', anova_p(ii, jj, kk));
                end
                
                title_str = sprintf('%s [%s]\n%s, %s', ...
                    trial_group_labels{kk}, model_name, rsq_str, anova_p_str);
            else
                title_str = trial_group_labels{kk};
            end
            title(gca, title_str, 'interpreter', 'none', 'fontsize', 8);
            
            % Axis labels
            if kk == 1
                tuning_curve_plot{kk}.xlabel('Temporal Frequency (Hz)');
                tuning_curve_plot{kk}.ylabel('Firing rate (Hz)');
            elseif kk == 2
                tuning_curve_plot{kk}.xlabel('Temporal Frequency (Hz)');
            elseif kk == 3
                tuning_curve_plot{kk}.xlabel('Equivalent bins');
            end
             
            % Find where the shuffled histogram axis should go and create it 
            pos         = plot_array.get_position(2*kk);
            h_ax2        = axes('units', 'centimeters', 'position', pos);
            
            % Object controlling the shuffled histogram plot
            shuff_histogram{kk} = TuningCurveHistogram(h_ax2);
            shuff_histogram{kk}.use_log_scale = true;  % Use log scale for model selection
            shuff_histogram{kk}.plot(tuning{ii}{jj}{kk});
        end
        
        % Plot V+T data (4th plot in position 8: row 2, col 2)
        if ~isempty(tuning_VplusT{ii}{jj})
            % Find where the axis should go and create it
            pos = plot_array.get_position(2*4-1);  % Position 7 (but we want 8 for col 2)
            pos = plot_array.get_position(2*4);  % Position 8
            % Adjust to make it same size as tuning plots (odd positions)
            pos_ref = plot_array.get_position(2*1-1);
            pos(3:4) = pos_ref(3:4);  % Copy width and height
            h_ax_VplusT = axes('units', 'centimeters', 'position', pos);
            hold(h_ax_VplusT, 'on');
            
            % Calculate mean and SEM for V+T tuning
            VplusT_mean = nanmean(tuning_VplusT{ii}{jj}.tuning, 2);
            VplusT_sem = nanstd(tuning_VplusT{ii}{jj}.tuning, [], 2) ./ sqrt(sum(~isnan(tuning_VplusT{ii}{jj}.tuning), 2));
            
            % Get bin centers (x-axis values)
            if isfield(tuning_VplusT{ii}{jj}, 'bin_centers')
                x_vals = tuning_VplusT{ii}{jj}.bin_centers;
            else
                x_vals = 1:length(VplusT_mean);
            end
            
            % Determine color from two-way ANOVA result
            if anova2_VT_VplusT_p(ii, jj) < 0.05
                main_col = [0.85, 0.32, 0.1];
            else
                main_col = 'k';
            end
            
            % Plot V+T tuning curve with SEM
            h_ebar = errorbar(h_ax_VplusT, x_vals, VplusT_mean, VplusT_sem);
            h_ebar.Marker = 'o';
            h_ebar.MarkerSize = 4;
            h_ebar.Color = main_col;
            h_ebar.LineStyle = 'none';
            h_ebar.LineWidth = 1.5;
            h_ebar.CapSize = 0;
            
            % Plot stationary baseline
            stat_mean = nanmean(tuning_VplusT{ii}{jj}.stationary_fr);
            plot(h_ax_VplusT, xlim(h_ax_VplusT), [stat_mean, stat_mean], '--', 'Color', [0.5, 0.5, 0.5]);
            
            % Title with two-way ANOVA p-value
            if ~isnan(anova2_VT_VplusT_p(ii, jj))
                title_str = sprintf('V+T (p=%.4f)', anova2_VT_VplusT_p(ii, jj));
            else
                title_str = 'V+T';
            end
            title(h_ax_VplusT, title_str, 'interpreter', 'none');
            
            % Axis labels
            xlabel(h_ax_VplusT, 'Temporal Frequency (Hz)');
            
            % Store axis limits for later alignment
            tuning_curve_plot{4}.xmin = min(x_vals);
            tuning_curve_plot{4}.xmax = max(x_vals);
            tuning_curve_plot{4}.ymax = max(VplusT_mean + VplusT_sem);
            
            % Store axis handle for later xlim/ylim adjustments
            tuning_curve_plot{4}.ax = h_ax_VplusT;
            
            % Note: No shuffled histogram for V+T since it's a derived computation
            % (not original data with shuffling statistics)
        end
        
        % Loop over all axes and find the min X, max X and max Y values
        mx              = inf;
        Mx              = -inf;
        My              = -inf;
        
        for kk = 1 : length(tuning_curve_plot)
            % Check if this is V+T (kk=4) or regular trial group (kk<=3)
            if kk == 4
                if isempty(tuning_VplusT{ii}{jj})
                    continue
                end
            else
                if isempty(tuning{ii}{jj}{kk})
                    continue
                end
            end
            mx = min(tuning_curve_plot{kk}.xmin, mx);
            Mx = max(tuning_curve_plot{kk}.xmax, Mx);
            My = max(tuning_curve_plot{kk}.ymax, My);
        end
        
        % Again loop over axes and set the X and Y limits to the same value
        for kk = 1 : length(tuning_curve_plot)
            % Check if this is V+T (kk=4) or regular trial group (kk<=3)
            if kk == 4
                if isempty(tuning_VplusT{ii}{jj})
                    continue
                end
                % For V+T, set limits directly on the axis
                xlim(tuning_curve_plot{kk}.ax, [mx, Mx]);
                ylim(tuning_curve_plot{kk}.ax, [0, My]);
            else
                if isempty(tuning{ii}{jj}{kk})
                    continue
                end
                % For regular plots, use the TuningCurvePlot methods
                tuning_curve_plot{kk}.xlim([mx, Mx]);
                tuning_curve_plot{kk}.ylim([0, My]);
            end
        end

        % Give whole figure a title (probe_id, cluster # and cluster region)
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                probe_id{ii, jj}, ...
                cluster_id(ii, jj), ...
                cluster_region{ii, jj}));
        
        % Temporarily save the figure to a .pdf to join later
        ctl.figs.save_fig_to_join();
    end
    
    % Join temporarily saved figures and save to a new .pdf
    fname = sprintf('%s_TF_tuning.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    
    % Delete temporarily saved .pdfs and clear MATLAB figure windows
    ctl.figs.clear_figs();
end

%% Plot average across clusters

% Get an average tuning curve for each trial type (including V+T)
trial_type_tuning_store = cell(1, length(trial_group_labels) + 1);  % +1 for V+T
stationary_store = cell(1, length(trial_group_labels) + 1);

% Loop over experiments
for ii = 1 : length(probe_ids)  
    
    % Loop over VISp clusters in that experiment
    for jj = 1 : length(tuning{ii})
        
        % Loop over trial group types
        for kk = 1 : length(trial_group_labels)
            
            % If there was no information in an experiment for a particular
            % trial type then skip this trial type
            if isempty(tuning{ii}{jj}{kk})
                continue
            end
            
            % If restricted mode, only include clusters tuned in both VT and V
            if restricted
                if (direction(ii, jj, 1) == 0) || (direction(ii, jj, 2) == 0)
                    continue
                end
            end
            
            % Average tuning curve for this cluster
            avg_tuning = nanmean(tuning{ii}{jj}{kk}.tuning, 2);
            avg_stationary = mean(tuning{ii}{jj}{kk}.stationary_fr);
            
            trial_type_tuning_store{kk} = [trial_type_tuning_store{kk}, avg_tuning];
            stationary_store{kk} = [stationary_store{kk}, avg_stationary];
        end
        
        % Add V+T data to averages
        if ~isempty(tuning_VplusT{ii}{jj})
            % Apply same restriction if enabled
            if restricted
                if (direction(ii, jj, 1) == 0) || (direction(ii, jj, 2) == 0)
                    continue
                end
            end
            
            avg_VplusT_tuning = nanmean(tuning_VplusT{ii}{jj}.tuning, 2);
            avg_VplusT_stationary = mean(tuning_VplusT{ii}{jj}.stationary_fr);
            
            trial_type_tuning_store{4} = [trial_type_tuning_store{4}, avg_VplusT_tuning];
            stationary_store{4} = [stationary_store{4}, avg_VplusT_stationary];
        end
    end
end

% Create A4 figure and object controlling the array of axes on the figure
h_fig                   = ctl.figs.a4figure('landscape');
plot_array              = PlotArray(2, 4);  % 2 rows, 4 columns for 8 subplots
plot_array.ax_size_cm   = 5.0;  % Larger axes for better visibility
plot_array.nx_total     = 4;  % 4 columns in underlying grid
plot_array.ny_total     = 2;  % 2 rows in underlying grid

% Cell arrays which will contain objects controlling the main
% tuning curve plot and the shuffled tuning curve histogram info
h_ax                    = {};

% Loop over trial group types (including V+T)
for kk = 1 : length(trial_group_labels) + 1

    % Find where the axis should go and create it
    pos         = plot_array.get_position(2*kk-1);
    h_ax{kk}    = axes('units', 'centimeters', 'position', pos);
    hold on;
    
    % Get mean and std across clusters for each TF bin
    avg_fr = mean(trial_type_tuning_store{kk}, 2);
    sd_fr = std(trial_type_tuning_store{kk}, [], 2);
    n = sum(~isnan(trial_type_tuning_store{kk}), 2);
    
    % Plot the errorbar
    h_ebar = errorbar(h_ax{kk}, avg_fr, sd_fr./sqrt(n));
    
    h_ebar.Marker = 'o';
    h_ebar.MarkerSize = 3;
    h_ebar.Color = 'k';
    h_ebar.CapSize = 0;
    
    % Mean and std of firing rate across clusters for stationary data
    avg_stat_fr = mean(stationary_store{kk}, 2);
    sd_stat_fr = std(stationary_store{kk}, [], 2);
    n_stat = sum(~isnan(stationary_store{kk}), 2);
    
    % Plot the stationary data
    stat_col = [0.5, 0.5, 0.5];
    scatter(h_ax{kk}, 0, avg_stat_fr, 5, stat_col);
    line(h_ax{kk}, [0, 0], [-1, 1]*sd_stat_fr/n_stat, 'color', stat_col)
    
    % Display the trial type as title
    if kk <= length(trial_group_labels)
        title(gca, trial_group_labels{kk}, 'interpreter', 'none', 'fontsize', 8);
    else
        title(gca, 'V+T', 'interpreter', 'none', 'fontsize', 8);
    end
    
    % Axis labels on last axis
    if kk == length(trial_group_labels) + 1
        xlabel('TF Bin #');
        ylabel('Firing rate (Hz)');
    end
    
    % Print the number of clusters
    n_bins = size(trial_type_tuning_store{kk}, 1);
    n_clusters_plotted = size(trial_type_tuning_store{kk}, 2);
    text(h_ax{kk}, n_bins, 0, sprintf('n = %i', n_clusters_plotted), ...
        'horizontalalignment', 'right', 'verticalalignment', 'bottom');
end

% Loop over all axes and find the min X, max X and max Y values
mx              = inf;
Mx              = -inf;
My              = -inf;

for kk = 1 : length(trial_group_labels) + 1
    mx = min(h_ax{kk}.XLim(1), mx);
    Mx = max(h_ax{kk}.XLim(2), Mx);
    My = max(h_ax{kk}.YLim(2), My);
end

% Again loop over axes and set the X and Y limits to the same value
for kk = 1 : length(trial_group_labels) + 1
    h_ax{kk}.XLim = [-1, Mx];
    h_ax{kk}.YLim = [0, My];
    h_ax{kk}.Box = 'off';
end

% Figure title
if restricted
    FigureTitle(h_fig, 'Averaged TF tuning curves (restricted to tuned clusters)');
else
    FigureTitle(h_fig, 'Averaged TF tuning curves');
end

% Save to a new .pdf
if restricted
    ctl.figs.save_fig('average_tf_tuning_each_condition_restricted.pdf');
else
    ctl.figs.save_fig('average_tf_tuning_each_condition.pdf');
end

%% Save model fit information to CSV

fprintf('\nSaving model fit information to CSV...\n');

% Prepare data for CSV export
model_csv_data = {};
row_idx = 1;

for ii = 1:length(probe_ids)
    for jj = 1:size(cluster_id, 2)
        if jj > length(tuning{ii})
            continue;
        end
        
        for kk = 1:length(trial_group_labels)
            if ~isempty(tuning{ii}{jj}{kk})
                tc = tuning{ii}{jj}{kk};
                
                % Basic identifiers
                model_csv_data{row_idx, 1} = probe_ids{ii};
                model_csv_data{row_idx, 2} = cluster_id(ii, jj);
                model_csv_data{row_idx, 3} = cluster_region{ii, jj};
                model_csv_data{row_idx, 4} = trial_group_labels{kk};
                
                % Model selection info
                model_csv_data{row_idx, 5} = tc.shuffled.best_model;
                
                % Get R² value
                if isfield(tc.shuffled, 'rsq')
                    model_csv_data{row_idx, 6} = tc.shuffled.rsq;
                elseif isfield(tc.shuffled, 'best_model_info')
                    model_csv_data{row_idx, 6} = tc.shuffled.best_model_info.rsq;
                else
                    model_csv_data{row_idx, 6} = nan;
                end
                
                % Calculate R² on mean (not stored, must be calculated)
                % R² on mean = how well the fitted model explains the mean tuning curve
                if isfield(tc, 'tuning') && ~isempty(tc.tuning) && isfield(tc, 'bin_centers')
                    % Get mean firing rate for each bin
                    mean_fr = nanmean(tc.tuning, 2);
                    
                    % Evaluate fitted model at bin centers
                    x_data = tc.bin_centers(:);
                    beta = tc.shuffled.best_model_info.beta;
                    
                    % Calculate predicted values based on model type
                    switch tc.shuffled.best_model
                        case 'linear'
                            y_pred = polyval(beta, x_data);
                        case 'quadratic'
                            y_pred = polyval(beta, x_data);
                        case 'cubic'
                            y_pred = polyval(beta, x_data);
                        case 'gaussian'
                            y_pred = beta(1) * exp(-(x_data - beta(2)).^2 / (2*beta(3)^2)) + beta(4);
                        case 'asymmetric_gaussian'
                            sigma = (x_data < beta(2)) * beta(3) + (x_data >= beta(2)) * beta(4);
                            y_pred = beta(1) * exp(-((x_data - beta(2)).^2) ./ (2 * sigma.^2));
                        case 'sigmoid'
                            y_pred = beta(1) ./ (1 + exp(-beta(2) * (x_data - beta(3)))) + beta(4);
                        otherwise
                            y_pred = nan(size(x_data));
                    end
                    
                    % Calculate R² on mean
                    ss_res = sum((mean_fr - y_pred).^2);
                    ss_tot = sum((mean_fr - nanmean(mean_fr)).^2);
                    rsq_on_mean = 1 - (ss_res / ss_tot);
                    model_csv_data{row_idx, 7} = rsq_on_mean;
                else
                    model_csv_data{row_idx, 7} = nan;
                end
                
                % Get BIC
                if isfield(tc.shuffled.best_model_info, 'bic')
                    model_csv_data{row_idx, 8} = tc.shuffled.best_model_info.bic;
                else
                    model_csv_data{row_idx, 8} = nan;
                end
                
                % Get shuffled statistics
                if isfield(tc.shuffled, 'rsq_shuff')
                    shuff_rsq = tc.shuffled.rsq_shuff;
                elseif isfield(tc.shuffled, 'fit_metric_shuff')
                    shuff_rsq = tc.shuffled.fit_metric_shuff;
                else
                    shuff_rsq = [];
                end
                
                if ~isempty(shuff_rsq)
                    model_csv_data{row_idx, 9} = prctile(shuff_rsq, 95);
                    model_csv_data{row_idx, 10} = tc.shuffled.p;
                    model_csv_data{row_idx, 11} = model_csv_data{row_idx, 6} > model_csv_data{row_idx, 9};
                else
                    model_csv_data{row_idx, 9} = nan;
                    model_csv_data{row_idx, 10} = nan;
                    model_csv_data{row_idx, 11} = false;
                end
                
                % Fit parameters (model-specific)
                beta = tc.shuffled.best_model_info.beta;
                switch tc.shuffled.best_model
                    case 'linear'
                        model_csv_data{row_idx, 12} = beta(1);  % slope
                        model_csv_data{row_idx, 13} = beta(2);  % intercept
                        model_csv_data{row_idx, 14} = nan;
                        model_csv_data{row_idx, 15} = nan;
                        model_csv_data{row_idx, 16} = nan;
                        model_csv_data{row_idx, 17} = nan;
                        if isfield(tc.shuffled, 'r')
                            model_csv_data{row_idx, 18} = tc.shuffled.r;
                        else
                            model_csv_data{row_idx, 18} = nan;
                        end
                    case 'gaussian'
                        model_csv_data{row_idx, 12} = nan;
                        model_csv_data{row_idx, 13} = nan;
                        model_csv_data{row_idx, 14} = beta(1);  % amplitude
                        model_csv_data{row_idx, 15} = beta(2);  % center (mu)
                        model_csv_data{row_idx, 16} = beta(3);  % sigma
                        model_csv_data{row_idx, 17} = nan;
                        model_csv_data{row_idx, 18} = nan;
                    case 'asymmetric_gaussian'
                        model_csv_data{row_idx, 12} = nan;
                        model_csv_data{row_idx, 13} = nan;
                        model_csv_data{row_idx, 14} = beta(1);  % amplitude (R_max)
                        model_csv_data{row_idx, 15} = beta(2);  % center (x_max)
                        model_csv_data{row_idx, 16} = beta(3);  % sigma_minus
                        model_csv_data{row_idx, 17} = beta(4);  % sigma_plus
                        model_csv_data{row_idx, 18} = nan;
                    case 'sigmoid'
                        model_csv_data{row_idx, 12} = beta(2);  % slope (k)
                        model_csv_data{row_idx, 13} = beta(4);  % baseline
                        model_csv_data{row_idx, 14} = beta(1);  % amplitude
                        model_csv_data{row_idx, 15} = beta(3);  % x50 (midpoint)
                        model_csv_data{row_idx, 16} = nan;
                        model_csv_data{row_idx, 17} = nan;
                        model_csv_data{row_idx, 18} = nan;
                    otherwise
                        model_csv_data{row_idx, 12} = nan;
                        model_csv_data{row_idx, 13} = nan;
                        model_csv_data{row_idx, 14} = nan;
                        model_csv_data{row_idx, 15} = nan;
                        model_csv_data{row_idx, 16} = nan;
                        model_csv_data{row_idx, 17} = nan;
                        model_csv_data{row_idx, 18} = nan;
                end
                
                % Additional info
                model_csv_data{row_idx, 19} = size(tc.tuning, 2);  % n_trials
                model_csv_data{row_idx, 20} = nanmean(tc.stationary_fr);  % mean_stationary_fr
                
                % ANOVA results
                model_csv_data{row_idx, 21} = anova_F(ii, jj, kk);
                model_csv_data{row_idx, 22} = anova_p(ii, jj, kk);
                model_csv_data{row_idx, 23} = anova_df(ii, jj, kk, 1);
                model_csv_data{row_idx, 24} = anova_df(ii, jj, kk, 2);
                
                row_idx = row_idx + 1;
            end
        end
    end
end

% Create table and save to CSV
if ~isempty(model_csv_data)
    model_csv_table = cell2table(model_csv_data, ...
        'VariableNames', {'probe_id', 'cluster_id', 'region', 'condition', ...
                         'best_model', 'rsq', 'rsq_on_mean', 'bic', ...
                         'shuffled_rsq_95th_pct', 'p_value_shuffle', 'is_significant', ...
                         'param_slope_or_k', 'param_intercept_or_baseline', ...
                         'param_amplitude', 'param_center', 'param_sigma_or_sigma_minus', ...
                         'param_sigma_plus', 'param_r_correlation', ...
                         'n_trials', 'mean_stationary_fr', ...
                         'anova_F', 'anova_p', 'anova_df_between', 'anova_df_within'});
    
    % Save in the figures directory
    model_csv_filename = fullfile(ctl.figs.curr_dir, 'tf_tuning_model_fits.csv');
    writetable(model_csv_table, model_csv_filename);
    fprintf('Model fit information saved to: %s\n', model_csv_filename);
    fprintf('Total rows: %d\n', size(model_csv_table, 1));
end

% Save ANOVA results to CSV file (keeping original for backward compatibility)
fprintf('\nSaving ANOVA-only results to CSV...\n');

% Prepare data for CSV export
csv_data = {};
row_idx = 1;

for ii = 1:length(probe_ids)
    for jj = 1:size(cluster_id, 2)
        if jj > length(tuning{ii})
            continue;
        end
        
        for kk = 1:length(trial_group_labels)
            if ~isempty(tuning{ii}{jj}{kk})
                csv_data{row_idx, 1} = probe_ids{ii};
                csv_data{row_idx, 2} = cluster_id(ii, jj);
                csv_data{row_idx, 3} = cluster_region{ii, jj};
                csv_data{row_idx, 4} = trial_group_labels{kk};
                csv_data{row_idx, 5} = anova_F(ii, jj, kk);
                csv_data{row_idx, 6} = anova_p(ii, jj, kk);
                csv_data{row_idx, 7} = anova_df(ii, jj, kk, 1);
                csv_data{row_idx, 8} = anova_df(ii, jj, kk, 2);
                row_idx = row_idx + 1;
            end
        end
    end
end

% Create table and save to CSV
if ~isempty(csv_data)
    csv_table = cell2table(csv_data, ...
        'VariableNames', {'probe_id', 'cluster_id', 'region', 'condition', ...
                         'F_statistic', 'p_value', 'df_between', 'df_within'});
    
    % Save in the figures directory
    csv_filename = fullfile(ctl.figs.curr_dir, 'tf_tuning_anova_results.csv');
    writetable(csv_table, csv_filename);
    fprintf('ANOVA results saved to: %s\n', csv_filename);
end

fprintf('\n=== TF tuning analysis complete ===\n');
fprintf('Individual cluster figures: <probe_id>_TF_tuning.pdf\n');
fprintf('Average figure: average_tf_tuning_each_condition.pdf\n');
fprintf('Model fits CSV: tf_tuning_model_fits.csv\n');
fprintf('ANOVA results CSV: tf_tuning_anova_results.csv\n');
