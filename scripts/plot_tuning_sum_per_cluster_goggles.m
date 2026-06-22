%% plot_tuning_sum_per_cluster_goggles.m
% VF + T = VT additive-tuning test, per cluster — goggles rebuild.
%
% For each of the top-20 Speed+TF-selected clusters, test whether the VT 
% tuning curve is the additive sum of the two unimodal contributions:
%   * Speed axis: observed VT speed-tuning vs
%       Model A (gain+offset): VT ~ a·T + b         [black]
%       Model B (additive):     VT ~ a·T + b·V + c   [red]
%   * TF axis: same approach
%
% Usage:
%   plot_tuning_sum_per_cluster_goggles()

function plot_tuning_sum_per_cluster_goggles()
    
    % Paths
    CSV_PATH = 'D:\mvelez\formatted_data\csvs\glm_model_comparison.csv';
    SPEED_TUNING_DIR = 'D:\mvelez\formatted_data\csvs\tuning_curves';
    TF_TUNING_DIR = 'D:\mvelez\formatted_data\csvs\tf_tuning_curves';
    ACCEL_TUNING_DIR = 'D:\mvelez\formatted_data\csvs\acceleration_tuning_curves';
    SF_TUNING_DIR = 'D:\mvelez\formatted_data\csvs\sf_tuning_curves';
    
    % Initialize controller for figure management
    ctl = RC2Analysis();
    save_figs = true;
    overwrite = true;
    
    % Parameters
    PROBE_FILTER = {'CAA-1124370_rec1_rec2_rec3', 'CAA-1124371_rec1_rec2_rec3'};
    
    % Select clusters with Speed as main variable
    clusters_speed = select_clusters_by_variable(CSV_PATH, 'Speed', PROBE_FILTER);
    fprintf('Selected %d clusters with Speed as main variable\n', height(clusters_speed));
    
    % Select clusters with TF as main variable
    clusters_tf = select_clusters_by_variable(CSV_PATH, 'TF', PROBE_FILTER);
    fprintf('Selected %d clusters with TF as main variable\n', height(clusters_tf));
    
    % Select clusters with Acceleration as main variable
    clusters_accel = select_clusters_by_variable(CSV_PATH, 'Acceleration', PROBE_FILTER);
    fprintf('Selected %d clusters with Acceleration as main variable\n', height(clusters_accel));
    
    % Select clusters with SF as main variable
    clusters_sf = select_clusters_by_variable(CSV_PATH, 'SF', PROBE_FILTER);
    fprintf('Selected %d clusters with SF as main variable\n', height(clusters_sf));
    
    % Load tuning data for all unique probes
    all_probes = unique([clusters_speed.probe_id; clusters_tf.probe_id; clusters_accel.probe_id; clusters_sf.probe_id]);
    tuning_data = struct();
    for i = 1:length(all_probes)
        probe = all_probes{i};
        fprintf('Loading tuning data for %s...\n', probe);
        
        spd_file = fullfile(SPEED_TUNING_DIR, [probe '.mat']);
        tf_file = fullfile(TF_TUNING_DIR, [probe '.mat']);
        accel_file = fullfile(ACCEL_TUNING_DIR, [probe '.mat']);
        sf_file = fullfile(SF_TUNING_DIR, [probe '.mat']);
        
        if exist(spd_file, 'file')
            tuning_data.(sanitize_fieldname(probe)).speed = load(spd_file);
        end
        if exist(tf_file, 'file')
            tuning_data.(sanitize_fieldname(probe)).tf = load(tf_file);
        end
        if exist(accel_file, 'file')
            tuning_data.(sanitize_fieldname(probe)).acceleration = load(accel_file);
        end
        if exist(sf_file, 'file')
            tuning_data.(sanitize_fieldname(probe)).sf = load(sf_file);
        end
    end
    
    % Process Speed-selected clusters
    results_speed = table();
    if ~isempty(clusters_speed)
        fprintf('\n=== Processing Speed-selected clusters ===\n');
        figure_dir = {'glm', 'exploration', 'subspace_population_goggles', 'speed_selected'};
        ctl.setup_figures(figure_dir, save_figs);
        results_speed = render_individual_clusters('Speed', clusters_speed, tuning_data, ctl);
        ctl.figs.join_figs('speed_clusters_merged.pdf', overwrite);
        ctl.figs.clear_figs();
        
        % Save results CSV
        out_csv_speed = fullfile(ctl.path_config.figure_dir, figure_dir{:}, 'speed_clusters_two_models_r2.csv');
        writetable(results_speed, out_csv_speed);
        fprintf('Saved Speed results: %s\n', out_csv_speed);
        
        % Print statistics
        fprintf('\nSpeed-selected clusters statistics:\n');
        print_statistics(results_speed);
    end
    
    % Process TF-selected clusters
    results_tf = table();
    if ~isempty(clusters_tf)
        fprintf('\n=== Processing TF-selected clusters ===\n');
        figure_dir = {'glm', 'exploration', 'subspace_population_goggles', 'tf_selected'};
        ctl.setup_figures(figure_dir, save_figs);
        results_tf = render_individual_clusters('TF', clusters_tf, tuning_data, ctl);
        ctl.figs.join_figs('tf_clusters_merged.pdf', overwrite);
        ctl.figs.clear_figs();
        
        % Save results CSV
        out_csv_tf = fullfile(ctl.path_config.figure_dir, figure_dir{:}, 'tf_clusters_two_models_r2.csv');
        writetable(results_tf, out_csv_tf);
        fprintf('Saved TF results: %s\n', out_csv_tf);
        
        % Print statistics
        fprintf('\nTF-selected clusters statistics:\n');
        print_statistics(results_tf);
    end
    
    % Process Acceleration-selected clusters
    results_accel = table();
    if ~isempty(clusters_accel)
        fprintf('\n=== Processing Acceleration-selected clusters ===\n');
        figure_dir = {'glm', 'exploration', 'subspace_population_goggles', 'accel_selected'};
        ctl.setup_figures(figure_dir, save_figs);
        results_accel = render_individual_clusters('Acceleration', clusters_accel, tuning_data, ctl);
        ctl.figs.join_figs('accel_clusters_merged.pdf', overwrite);
        ctl.figs.clear_figs();
        
        % Save results CSV
        out_csv_accel = fullfile(ctl.path_config.figure_dir, figure_dir{:}, 'accel_clusters_two_models_r2.csv');
        writetable(results_accel, out_csv_accel);
        fprintf('Saved Acceleration results: %s\n', out_csv_accel);
        
        % Print statistics
        fprintf('\nAcceleration-selected clusters statistics:\n');
        print_statistics(results_accel);
    end
    
    % Process SF-selected clusters
    results_sf = table();
    if ~isempty(clusters_sf)
        fprintf('\n=== Processing SF-selected clusters ===\n');
        figure_dir = {'glm', 'exploration', 'subspace_population_goggles', 'sf_selected'};
        ctl.setup_figures(figure_dir, save_figs);
        results_sf = render_individual_clusters('SF', clusters_sf, tuning_data, ctl);
        ctl.figs.join_figs('sf_clusters_merged.pdf', overwrite);
        ctl.figs.clear_figs();
        
        % Save results CSV
        out_csv_sf = fullfile(ctl.path_config.figure_dir, figure_dir{:}, 'sf_clusters_two_models_r2.csv');
        writetable(results_sf, out_csv_sf);
        fprintf('Saved SF results: %s\n', out_csv_sf);
        
        % Print statistics
        fprintf('\nSF-selected clusters statistics:\n');
        print_statistics(results_sf);
    end
    
    % Create combined beta coefficient scatter plot
    if ~isempty(results_speed) || ~isempty(results_tf) || ~isempty(results_accel) || ~isempty(results_sf)
        fprintf('\n=== Creating combined beta scatter plot ===\n');
        figure_dir = {'glm', 'exploration', 'subspace_population_goggles'};
        ctl.setup_figures(figure_dir, save_figs);
        plot_beta_scatter_combined(results_speed, results_tf, results_accel, results_sf, ctl);
        ctl.figs.join_figs('all_clusters_beta_scatter.pdf', overwrite);
        ctl.figs.clear_figs();
    end
    
    fprintf('\n=== Analysis complete ===\n');
end

%% select_clusters_by_variable: Select clusters with specific variable as main
function clusters = select_clusters_by_variable(csv_path, variable, probe_filter)
    % variable: 'Speed' or 'TF'
    % Returns all clusters that have the specified variable selected
    
    % Read CSV
    opts = detectImportOptions(csv_path, 'TextType', 'string');
    df = readtable(csv_path, opts);
    
    % Filter for specified probes
    if ~isempty(probe_filter)
        probe_mask = false(height(df), 1);
        for i = 1:length(probe_filter)
            probe_mask = probe_mask | strcmp(df.probe_id, probe_filter{i});
        end
        df = df(probe_mask, :);
    end
    
    % Filter for clusters with time_selected_vars starting with the specified variable
    if ~ismember('time_selected_vars', df.Properties.VariableNames)
        error('Column time_selected_vars not found in CSV');
    end
    
    % Fill missing values with empty string
    sv = df.time_selected_vars;
    sv(ismissing(sv)) = "";
    
    % Select clusters where time_selected_vars starts with the variable (main effect)
    has_variable = startsWith(sv, variable);
    df_filtered = df(has_variable, :);
    
    if isempty(df_filtered)
        warning('No clusters found with %s as main variable', variable);
        clusters = table();
        return;
    end
    
    % Sort by delta (best performers first)
    if ismember('time_delta_selected_vs_null', df_filtered.Properties.VariableNames)
        df_filtered = sortrows(df_filtered, 'time_delta_selected_vs_null', 'descend');
    end
    
    clusters = df_filtered(:, {'probe_id', 'cluster_id'});
    
    fprintf('Selected %d clusters with %s:\n', height(clusters), variable);
    for i = 1:min(10, height(clusters))  % Show first 10
        fprintf('  %s / %d\n', clusters.probe_id{i}, clusters.cluster_id(i));
    end
    if height(clusters) > 10
        fprintf('  ... and %d more\n', height(clusters) - 10);
    end
end

%% mean_tuning: Get mean tuning curve for one condition/cluster/variable
function [centres, mean_curve] = mean_tuning(tuning_data, probe_id, cluster_id, condition, var_type)
    % var_type: 'Speed', 'TF', 'Acceleration', or 'SF'
    % condition: 'VT', 'T_Vstatic', 'V'
    
    probe_field = sanitize_fieldname(probe_id);
    if ~isfield(tuning_data, probe_field)
        centres = [];
        mean_curve = [];
        return;
    end
    
    % Select appropriate data structure
    if strcmp(var_type, 'Speed')
        D = tuning_data.(probe_field).speed;
    elseif strcmp(var_type, 'TF')
        D = tuning_data.(probe_field).tf;
    elseif strcmp(var_type, 'Acceleration')
        if ~isfield(tuning_data.(probe_field), 'acceleration')
            centres = [];
            mean_curve = [];
            return;
        end
        D = tuning_data.(probe_field).acceleration;
    elseif strcmp(var_type, 'SF')
        if ~isfield(tuning_data.(probe_field), 'sf')
            centres = [];
            mean_curve = [];
            return;
        end
        D = tuning_data.(probe_field).sf;
    else
        centres = [];
        mean_curve = [];
        return;
    end
    
    % Find condition
    cond_idx = find(strcmp(D.trial_groups, condition), 1);
    if isempty(cond_idx)
        centres = [];
        mean_curve = [];
        return;
    end
    
    tc_array = D.tuning_curves{cond_idx};
    if isempty(tc_array)
        centres = [];
        mean_curve = [];
        return;
    end
    
    % Find cluster
    tc_cluster_ids = arrayfun(@(x) double(x.cluster_id), tc_array);
    tc_idx = find(tc_cluster_ids == cluster_id, 1);
    if isempty(tc_idx)
        centres = [];
        mean_curve = [];
        return;
    end
    
    tc = tc_array(tc_idx);
    centres = double(tc.bin_centers(:))';
    
    % Mean over trials (dimension 2)
    mean_curve = nanmean(tc.tuning, 2)';  % tuning is (n_bins x n_trials)
end

%% calc_r2: Calculate R² between observed and predicted
function r2_val = calc_r2(y, pred)
    ok = isfinite(y) & isfinite(pred);
    if sum(ok) < 2
        r2_val = NaN;
        return;
    end
    
    y_valid = y(ok);
    pred_valid = pred(ok);
    
    ss_res = sum((y_valid - pred_valid).^2);
    ss_tot = sum((y_valid - mean(y_valid)).^2);
    
    if ss_tot > 0
        r2_val = 1.0 - ss_res / ss_tot;
    else
        r2_val = NaN;
    end
end

%% fit_model: Linear regression with multiple columns
function [coef, pred, r2_val] = fit_model(cols, y)
    % cols: cell array of column vectors
    % Add intercept
    design = [cell2mat(cols), ones(length(y), 1)];
    
    % Filter out rows with NaN
    ok = all(isfinite(design), 2) & isfinite(y);
    
    if sum(ok) < size(design, 2) + 1
        coef = [];
        pred = [];
        r2_val = NaN;
        return;
    end
    
    % Least squares fit
    coef = design(ok, :) \ y(ok);
    
    % Predict on full grid
    pred = design * coef;
    
    % Calculate R²
    r2_val = calc_r2(y, pred);
end

%% models_for: Build condition curves and model predictions
function res = models_for(tuning_data, probe_id, cluster_id, axis)
    % axis: 'Speed', 'TF', 'Acceleration', or 'SF'
    
    % Get VT tuning (observed) on the axis grid
    [grid, vt] = mean_tuning(tuning_data, probe_id, cluster_id, 'VT', axis);
    
    % Get T and V curves (on their respective grids, but quantile-matched to VT)
    [~, t_curve] = mean_tuning(tuning_data, probe_id, cluster_id, 'T_Vstatic', 'Speed');
    [~, v_curve] = mean_tuning(tuning_data, probe_id, cluster_id, 'V', 'TF');
    
    if isempty(grid) || isempty(t_curve) || isempty(v_curve)
        res = [];
        return;
    end
    
    % Check alignment (all should have same number of bins due to quantile binning)
    n = length(grid);
    if length(t_curve) ~= n || length(v_curve) ~= n
        res = [];
        return;
    end
    
    % Select primary predictor based on axis
    if strcmp(axis, 'Speed') || strcmp(axis, 'Acceleration')
        primary = t_curve(:);   % T_Vstatic for Speed and Acceleration
    else % TF or SF
        primary = v_curve(:);   % V for TF and SF
    end
    
    vt = vt(:);
    
    % Fit Model A: gain+offset (one regressor + intercept)
    [coef_a, pred_a, r2_a] = fit_model({primary}, vt);
    
    % Fit Model B: additive (two regressors + intercept)
    % Always fit as [V, T] so beta1=V, beta2=T
    [coef_b, pred_b, r2_b] = fit_model({v_curve(:), t_curve(:)}, vt);
    
    if isempty(coef_a) || isempty(coef_b)
        res = [];
        return;
    end
    
    % Package results
    res = struct();
    res.grid = grid(:);
    res.vt = vt;
    res.t_curve = t_curve(:);
    res.v_curve = v_curve(:);
    res.predA = pred_a;
    res.r2A = r2_a;
    res.coefA = coef_a;  % [beta1, beta0]
    res.predB = pred_b;
    res.r2B = r2_b;
    res.coefB = coef_b;  % [beta1, beta2, beta0]
end

%% render_individual_clusters: Create individual PDFs for each cluster and merge
function results_table = render_individual_clusters(selection_type, clusters, tuning_data, ctl)
    % selection_type: 'Speed', 'TF', 'Acceleration', or 'SF' (which variable was used for selection)
    % Uses RC2Analysis controller for figure management and merging
    
    n_clusters = height(clusters);
    rows = {};
    
    fprintf('Rendering %d clusters...\n', n_clusters);
    
    for idx = 1:n_clusters
        probe_id = clusters.probe_id{idx};
        cluster_id = clusters.cluster_id(idx);
        
        fprintf('  [%d/%d] %s cl%d\n', idx, n_clusters, probe_id, cluster_id);
        
        % Create figure for this cluster using controller
        h_fig = ctl.figs.a4figure();
        set(h_fig, 'Position', [100, 100, 400, 400]);
        
        % Plot only the axis corresponding to selection_type
        res = models_for(tuning_data, probe_id, cluster_id, selection_type);
        if ~isempty(res)
            plot_single_cluster(res, probe_id, cluster_id, selection_type);
            axis_name = sprintf('%s_to_VT', selection_type);
            rows{end+1} = struct('probe_id', string(probe_id), 'cluster_id', cluster_id, ...
                'axis', string(axis_name), ...
                'r2_gain_offset', res.r2A, 'r2_additive', res.r2B, ...
                'gain_offset_beta1', res.coefA(1), 'gain_offset_beta0', res.coefA(2), ...
                'additive_betaV', res.coefB(1), 'additive_betaT', res.coefB(2), 'additive_beta0', res.coefB(3));
        else
            title(sprintf('%s - No data', selection_type));
        end
        
        % Overall title with beta coefficients
        beta_str = '';
        if ~isempty(res)
            beta_str = sprintf('β_V=%.2f β_T=%.2f', res.coefB(1), res.coefB(2));
        end
        sgtitle(sprintf('%s selected: %s cluster %d | %s\nVF+T = VT additive test | green T · gold V · blue VT | black -- Model A · red — Model B', ...
            selection_type, probe_id, cluster_id, beta_str), 'FontSize', 10, 'FontWeight', 'bold');
        
        % Save figure to join queue
        ctl.figs.save_fig_to_join();
    end
    
    % Convert results to table
    if ~isempty(rows)
        results_table = struct2table(vertcat(rows{:}));
    else
        results_table = table();
    end
end

%% plot_single_cluster: Plot tuning curves and models for one cluster/axis
function plot_single_cluster(res, probe_id, cluster_id, axis)
    hold on;
    
    % Plot tuning curves (data)
    plot(res.grid, res.t_curve, '-', 'Color', [0 0.5 0], 'LineWidth', 1.5, ...
        'DisplayName', 'T');
    plot(res.grid, res.v_curve, '-', 'Color', [0.85 0.65 0.13], 'LineWidth', 1.5, ...
        'DisplayName', 'V');
    plot(res.grid, res.vt, 'o-', 'Color', [0 0.45 0.74], 'MarkerSize', 4, ...
        'LineWidth', 1.5, 'DisplayName', 'VT');
    
    % Plot model predictions
    plot(res.grid, res.predA, '--', 'Color', 'k', 'LineWidth', 2, ...
        'DisplayName', sprintf('A gain+offset r²=%.2f', res.r2A));
    plot(res.grid, res.predB, '-', 'Color', [0.64 0.08 0.18], 'LineWidth', 2.2, ...
        'DisplayName', sprintf('B additive r²=%.2f', res.r2B));
    
    % Formatting
    if strcmp(axis, 'Speed')
        xlabel_txt = 'speed (cm/s)';
    elseif strcmp(axis, 'TF')
        xlabel_txt = 'TF (Hz)';
    elseif strcmp(axis, 'Acceleration')
        xlabel_txt = 'acceleration (cm/s²)';
        % Make x-axis symmetric around 0
        current_xlim = xlim;
        max_abs = max(abs(current_xlim));
        xlim([-max_abs, max_abs]);
    elseif strcmp(axis, 'SF')
        xlabel_txt = 'SF (cyc/deg)';
        % Start from 0
        current_xlim = xlim;
        xlim([0, current_xlim(2)]);
    else
        xlabel_txt = axis;
    end
    
    title(sprintf('%s tuning', axis), 'FontSize', 10);
    xlabel(xlabel_txt, 'FontSize', 9);
    ylabel('FR (Hz)', 'FontSize', 9);
    set(gca, 'FontSize', 8);
    lgd = legend('FontSize', 7, 'Location', 'best');
    lgd.Box = 'off';
    hold off;
end

%% plot_beta_scatter_combined: Create scatter plot of beta coefficients with marginal histograms for all clusters
function plot_beta_scatter_combined(results_speed, results_tf, results_accel, results_sf, ctl)
    % Combine results from all 4 selection types (Speed, TF, Acceleration, SF)
    betaV_all = [];
    betaT_all = [];
    colors_all = [];
    
    % Color scheme: Speed=blue, TF=orange, Acceleration=green, SF=purple
    color_speed = [0 0.45 0.74];
    color_tf = [0.85 0.33 0.1];
    color_accel = [0.47 0.67 0.19];
    color_sf = [0.49 0.18 0.56];
    
    % Extract Speed data
    n_speed = 0;
    if ~isempty(results_speed)
        betaV_speed = results_speed.additive_betaV;
        betaT_speed = results_speed.additive_betaT;
        valid_speed = isfinite(betaV_speed) & isfinite(betaT_speed);
        betaV_all = [betaV_all; betaV_speed(valid_speed)];
        betaT_all = [betaT_all; betaT_speed(valid_speed)];
        n_speed = sum(valid_speed);
        colors_all = [colors_all; repmat(color_speed, n_speed, 1)];
    end
    
    % Extract TF data
    n_tf = 0;
    if ~isempty(results_tf)
        betaV_tf = results_tf.additive_betaV;
        betaT_tf = results_tf.additive_betaT;
        valid_tf = isfinite(betaV_tf) & isfinite(betaT_tf);
        betaV_all = [betaV_all; betaV_tf(valid_tf)];
        betaT_all = [betaT_all; betaT_tf(valid_tf)];
        n_tf = sum(valid_tf);
        colors_all = [colors_all; repmat(color_tf, n_tf, 1)];
    end
    
    % Extract Acceleration data
    n_accel = 0;
    if ~isempty(results_accel)
        betaV_accel = results_accel.additive_betaV;
        betaT_accel = results_accel.additive_betaT;
        valid_accel = isfinite(betaV_accel) & isfinite(betaT_accel);
        betaV_all = [betaV_all; betaV_accel(valid_accel)];
        betaT_all = [betaT_all; betaT_accel(valid_accel)];
        n_accel = sum(valid_accel);
        colors_all = [colors_all; repmat(color_accel, n_accel, 1)];
    end
    
    % Extract SF data
    n_sf = 0;
    if ~isempty(results_sf)
        betaV_sf = results_sf.additive_betaV;
        betaT_sf = results_sf.additive_betaT;
        valid_sf = isfinite(betaV_sf) & isfinite(betaT_sf);
        betaV_all = [betaV_all; betaV_sf(valid_sf)];
        betaT_all = [betaT_all; betaT_sf(valid_sf)];
        n_sf = sum(valid_sf);
        colors_all = [colors_all; repmat(color_sf, n_sf, 1)];
    end
    
    if length(betaV_all) < 2
        fprintf('Insufficient data for beta scatter plot\n');
        return;
    end
    
    % Create figure
    h_fig = ctl.figs.a4figure();
    set(h_fig, 'Position', [100, 100, 800, 800]);
    
    % Calculate quartiles
    quartiles_V = quantile(betaV_all, [0.25, 0.5, 0.75]);
    quartiles_T = quantile(betaT_all, [0.25, 0.5, 0.75]);
    
    % Define subplot positions [left bottom width height]
    pos_scatter = [0.15, 0.15, 0.65, 0.65];  % Main scatter plot
    pos_hist_top = [0.15, 0.82, 0.65, 0.15];  % Top histogram
    pos_hist_right = [0.82, 0.15, 0.15, 0.65];  % Right histogram
    
    % Main scatter plot
    subplot('Position', pos_scatter);
    hold on;
    
    % Set axis limits: same range for both axes, including all data
    min_val = min([min(betaV_all), min(betaT_all)]);
    max_val = max([max(betaV_all), max(betaT_all)]);
    xlim([min_val, max_val]);
    ylim([min_val, max_val]);
    
    % Unity line
    plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Plot scatter points with separate handles for legend
    handles = [];
    labels = {};
    
    idx_start = 1;
    if n_speed > 0
        idx_end = idx_start + n_speed - 1;
        h = scatter(betaV_all(idx_start:idx_end), betaT_all(idx_start:idx_end), 50, ...
            color_speed, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
        handles(end+1) = h;
        labels{end+1} = sprintf('Speed (n=%d)', n_speed);
        idx_start = idx_end + 1;
    end
    
    if n_tf > 0
        idx_end = idx_start + n_tf - 1;
        h = scatter(betaV_all(idx_start:idx_end), betaT_all(idx_start:idx_end), 50, ...
            color_tf, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
        handles(end+1) = h;
        labels{end+1} = sprintf('TF (n=%d)', n_tf);
        idx_start = idx_end + 1;
    end
    
    if n_accel > 0
        idx_end = idx_start + n_accel - 1;
        h = scatter(betaV_all(idx_start:idx_end), betaT_all(idx_start:idx_end), 50, ...
            color_accel, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
        handles(end+1) = h;
        labels{end+1} = sprintf('Acceleration (n=%d)', n_accel);
        idx_start = idx_end + 1;
    end
    
    if n_sf > 0
        idx_end = idx_start + n_sf - 1;
        h = scatter(betaV_all(idx_start:idx_end), betaT_all(idx_start:idx_end), 50, ...
            color_sf, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.6);
        handles(end+1) = h;
        labels{end+1} = sprintf('SF (n=%d)', n_sf);
    end
    
    xlabel('\beta_V (Visual TF coefficient)', 'FontSize', 10);
    ylabel('\beta_T (Speed coefficient)', 'FontSize', 10);
    title(sprintf('All Clusters: Additive Model Coefficients (n=%d)', length(betaV_all)), 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
    
    % Add legend
    if ~isempty(handles)
        lgd = legend(handles, labels, 'Location', 'northwest', 'FontSize', 8);
        lgd.ItemTokenSize = [8, 8];
        lgd.Box = 'off';
    end
    
    axis square;
    hold off;
    
    % Top histogram (betaV)
    subplot('Position', pos_hist_top);
    hold on;
    
    edges = linspace(min_val, max_val, 30);
    
    % Plot histograms for each type
    idx_start = 1;
    if n_speed > 0
        idx_end = idx_start + n_speed - 1;
        histogram(betaV_all(idx_start:idx_end), edges, 'FaceColor', color_speed, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        idx_start = idx_end + 1;
    end
    if n_tf > 0
        idx_end = idx_start + n_tf - 1;
        histogram(betaV_all(idx_start:idx_end), edges, 'FaceColor', color_tf, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        idx_start = idx_end + 1;
    end
    if n_accel > 0
        idx_end = idx_start + n_accel - 1;
        histogram(betaV_all(idx_start:idx_end), edges, 'FaceColor', color_accel, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        idx_start = idx_end + 1;
    end
    if n_sf > 0
        idx_end = idx_start + n_sf - 1;
        histogram(betaV_all(idx_start:idx_end), edges, 'FaceColor', color_sf, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
    
    % Quartile lines
    ylims_hist = ylim;
    for q = quartiles_V
        plot([q, q], ylims_hist, 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
    xlim([min_val, max_val]);
    set(gca, 'XTick', [], 'FontSize', 8);
    ylabel('Count', 'FontSize', 8);
    box on;
    hold off;
    
    % Right histogram (betaT)
    subplot('Position', pos_hist_right);
    hold on;
    
    % Plot histograms for each type
    idx_start = 1;
    if n_speed > 0
        idx_end = idx_start + n_speed - 1;
        histogram(betaT_all(idx_start:idx_end), edges, 'Orientation', 'horizontal', ...
            'FaceColor', color_speed, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        idx_start = idx_end + 1;
    end
    if n_tf > 0
        idx_end = idx_start + n_tf - 1;
        histogram(betaT_all(idx_start:idx_end), edges, 'Orientation', 'horizontal', ...
            'FaceColor', color_tf, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        idx_start = idx_end + 1;
    end
    if n_accel > 0
        idx_end = idx_start + n_accel - 1;
        histogram(betaT_all(idx_start:idx_end), edges, 'Orientation', 'horizontal', ...
            'FaceColor', color_accel, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        idx_start = idx_end + 1;
    end
    if n_sf > 0
        idx_end = idx_start + n_sf - 1;
        histogram(betaT_all(idx_start:idx_end), edges, 'Orientation', 'horizontal', ...
            'FaceColor', color_sf, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
    
    % Quartile lines
    xlims_hist = xlim;
    for q = quartiles_T
        plot(xlims_hist, [q, q], 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
    ylim([min_val, max_val]);
    set(gca, 'YTick', [], 'FontSize', 8);
    xlabel('Count', 'FontSize', 8);
    box on;
    hold off;
    
    % Save figure to join queue
    ctl.figs.save_fig_to_join();
    
    fprintf('Created combined beta scatter plot with %d total clusters\n', length(betaV_all));
end

%% print_statistics: Print summary statistics for results
function print_statistics(results)
    if isempty(results)
        fprintf('  No results to summarize\n');
        return;
    end
    
    % All results now have single axis per selection
    beats = sum(results.r2_additive > results.r2_gain_offset);
    fprintf('  r² gain+offset %.3f → additive %.3f  (B>A %d/%d)\n', ...
        median(results.r2_gain_offset, 'omitnan'), ...
        median(results.r2_additive, 'omitnan'), beats, height(results));
end

%% Utility functions
function name = sanitize_fieldname(str)
    % Convert string to valid MATLAB field name
    name = matlab.lang.makeValidName(char(str));
end
