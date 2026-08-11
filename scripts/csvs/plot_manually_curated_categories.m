% Plot manually curated categories
%
% This script reads the manually_curated.csv file and plots the fitted
% tuning curves for clusters in each of the 3 categories.
%
% The CSV has alternating columns: Category 1, cluster_id, Category 2, cluster_id, Category 3, cluster_id

%% Configuration

% Path to manually curated CSV
figure_dir = {'tf_tuning_curves'};
curated_csv_filename = 'manually_curated.csv';

%% Load data

ctl = RC2Analysis();

% Setup directory
ctl.setup_figures(figure_dir, true);  % Enable figure saving

% Load manually curated CSV
curated_csv_path = fullfile(ctl.figs.curr_dir, curated_csv_filename);
if ~exist(curated_csv_path, 'file')
    error('Manually curated CSV file not found: %s', curated_csv_path);
end

fprintf('Loading manually curated data from:\n  %s\n\n', curated_csv_path);

% Read with options to preserve column names exactly as they are
opts = detectImportOptions(curated_csv_path);
curated_data = readtable(curated_csv_path, opts);

% Display column names to verify structure
fprintf('CSV columns: %s\n\n', strjoin(curated_data.Properties.VariableNames, ', '));

%% Parse the CSV structure

categories = cell(1, 3);
category_names = {'Category 1', 'Category 2', 'Category 3'};

for cat_idx = 1:3
    cat_name = sprintf('Category%d', cat_idx);
    
    % Find the probe_id column (Category 1, Category 2, or Category 3)
    probe_col_idx = find(strcmp(curated_data.Properties.VariableNames, cat_name));
    
    if isempty(probe_col_idx)
        warning('Could not find column for %s', cat_name);
        categories{cat_idx} = table();
        continue;
    end
    
    % The cluster_id column should be right after the probe_id column
    cluster_col_idx = probe_col_idx + 1;
    
    if cluster_col_idx > width(curated_data)
        warning('Could not find cluster_id column for %s', cat_name);
        categories{cat_idx} = table();
        continue;
    end
    
    % Extract probe_id and cluster_id
    probe_ids = curated_data{:, probe_col_idx};
    cluster_ids = curated_data{:, cluster_col_idx};
    
    % Handle different data types
    if iscell(probe_ids)
        % Remove empty cells
        valid_mask = ~cellfun(@isempty, probe_ids);
    else
        % Numeric or other - check for NaN
        valid_mask = ~isnan(probe_ids);
        probe_ids = cellstr(num2str(probe_ids));
    end
    
    % Also check cluster_ids are valid
    if ~isnumeric(cluster_ids)
        cluster_ids = str2double(cluster_ids);
    end
    valid_mask = valid_mask & ~isnan(cluster_ids);
    
    categories{cat_idx} = table(probe_ids(valid_mask), cluster_ids(valid_mask), ...
        'VariableNames', {'probe_id', 'cluster_id'});
    
    fprintf('Category %d: %d clusters\n', cat_idx, height(categories{cat_idx}));
end

fprintf('\n');

%% Load tuning data for all unique probes

% Get all unique probe_ids across all categories
all_probe_ids = {};
for cat_idx = 1:3
    if ~isempty(categories{cat_idx})
        all_probe_ids = [all_probe_ids; categories{cat_idx}.probe_id];
    end
end
unique_probes = unique(all_probe_ids);

fprintf('Loading tuning data for %d unique probes...\n', length(unique_probes));

% Build lookup structure: probe_id -> cluster_id -> condition -> tuning_data
tuning_lookup = struct();

for ii = 1:length(unique_probes)
    probe_id = unique_probes{ii};
    
    % Load tuning curves for this probe
    tbl = ctl.load_tf_tuning_curves(probe_id);
    
    % Store in lookup structure
    tuning_lookup.(matlab.lang.makeValidName(probe_id)) = struct();
    
    for kk = 1:length(tbl.trial_groups)
        trial_group = tbl.trial_groups{kk};
        tuning_curves = tbl.tuning_curves{kk};
        
        for jj = 1:length(tuning_curves)
            tc = tuning_curves(jj);
            cluster_key = sprintf('cluster_%d', tc.cluster_id);
            
            if ~isfield(tuning_lookup.(matlab.lang.makeValidName(probe_id)), cluster_key)
                tuning_lookup.(matlab.lang.makeValidName(probe_id)).(cluster_key) = struct();
            end
            
            tuning_lookup.(matlab.lang.makeValidName(probe_id)).(cluster_key).(trial_group) = tc;
        end
    end
end

fprintf('Tuning data loaded.\n\n');

%% Plot each category

fprintf('Creating category visualization figures...\n');

for cat_idx = 1:3
    cat_name = category_names{cat_idx};
    cat_clusters = categories{cat_idx};
    n_cat_clusters = height(cat_clusters);
    
    fprintf('  Plotting %s (%d clusters)\n', cat_name, n_cat_clusters);
    
    if n_cat_clusters == 0
        fprintf('    Skipping - no clusters in this category.\n');
        continue;
    end
    
    % Generate distinct colors for each cluster
    colors = lines(n_cat_clusters);
    
    % Create figure
    h_fig = ctl.figs.a4figure('landscape');
    
    % Create subplots: V (left), VT (middle), and T_Vstatic (right)
    subplot(1, 3, 1);
    h_ax_V = gca;
    hold on;
    title('Condition V');
    xlabel('Temporal Frequency (Hz)');
    ylabel('Firing Rate (Hz)');
    
    subplot(1, 3, 2);
    h_ax_VT = gca;
    hold on;
    title('Condition VT');
    xlabel('Temporal Frequency (Hz)');
    ylabel('Firing Rate (Hz)');
    
    subplot(1, 3, 3);
    h_ax_Tvstatic = gca;
    hold on;
    title('Condition T\_Vstatic');
    xlabel('Temporal Frequency (Hz)');
    ylabel('Firing Rate (Hz)');
    
    % Track x-axis range for axis alignment
    all_x_min = inf;
    all_x_max = -inf;
    
    % Plot individual cluster curves with unique colors
    for ii = 1:n_cat_clusters
        probe_id = cat_clusters.probe_id{ii};
        cluster_id = cat_clusters.cluster_id(ii);
        
        probe_key = matlab.lang.makeValidName(probe_id);
        cluster_key = sprintf('cluster_%d', cluster_id);
        
        % Check if data exists
        if ~isfield(tuning_lookup, probe_key) || ...
           ~isfield(tuning_lookup.(probe_key), cluster_key)
            fprintf('    Warning: No tuning data for %s cluster %d\n', probe_id, cluster_id);
            continue;
        end
        
        cluster_data = tuning_lookup.(probe_key).(cluster_key);
        cluster_color = colors(ii, :);
        
        % Plot V condition
        if isfield(cluster_data, 'V')
            tc_V = cluster_data.V;
            [x_fit, y_fit] = evaluate_fitted_curve(tc_V);
            
            if ~isempty(y_fit)
                plot(h_ax_V, x_fit, y_fit, 'Color', cluster_color, 'LineWidth', 1.5);
                all_x_min = min([all_x_min, min(x_fit)]);
                all_x_max = max([all_x_max, max(x_fit)]);
            end
        end
        
        % Plot VT condition
        if isfield(cluster_data, 'VT')
            tc_VT = cluster_data.VT;
            [x_fit, y_fit] = evaluate_fitted_curve(tc_VT);
            
            if ~isempty(y_fit)
                plot(h_ax_VT, x_fit, y_fit, 'Color', cluster_color, 'LineWidth', 1.5);
                all_x_min = min([all_x_min, min(x_fit)]);
                all_x_max = max([all_x_max, max(x_fit)]);
            end
        end
        
        % Plot T_Vstatic condition
        if isfield(cluster_data, 'T_Vstatic')
            tc_Tvstatic = cluster_data.T_Vstatic;
            [x_fit, y_fit] = evaluate_fitted_curve(tc_Tvstatic);
            
            if ~isempty(y_fit)
                plot(h_ax_Tvstatic, x_fit, y_fit, 'Color', cluster_color, 'LineWidth', 1.5);
                all_x_min = min([all_x_min, min(x_fit)]);
                all_x_max = max([all_x_max, max(x_fit)]);
            end
        end
    end
    
    % Make axes equal
    if isfinite(all_x_min) && isfinite(all_x_max)
        xlim(h_ax_V, [all_x_min, all_x_max]);
        xlim(h_ax_VT, [all_x_min, all_x_max]);
        xlim(h_ax_Tvstatic, [all_x_min, all_x_max]);
        
        y_max = max([ylim(h_ax_V), ylim(h_ax_VT), ylim(h_ax_Tvstatic)]);
        ylim(h_ax_V, [0, y_max]);
        ylim(h_ax_VT, [0, y_max]);
        ylim(h_ax_Tvstatic, [0, y_max]);
    end
    
    % Add title with category and cluster count
    sgtitle(sprintf('%s | n = %d clusters', cat_name, n_cat_clusters), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % Save figure
    fig_filename = sprintf('tf_tuning_manually_curated_category%d.pdf', cat_idx);
    full_path = fullfile(ctl.figs.curr_dir, fig_filename);
    if exist(full_path, 'file')
        delete(full_path);
    end
    ctl.figs.save_fig(fig_filename);
    fprintf('    Saved: %s\n', full_path);
end

fprintf('\nAll category figures saved.\n');
fprintf('=== Analysis complete ===\n');


%% Helper function to evaluate fitted curve (copied from compare_vt_vs_tvstatic_filtered.m)

function [x_fit, y_fit] = evaluate_fitted_curve(tuning)
    % Evaluate the fitted model for a tuning curve
    
    % Get x range
    if isfield(tuning, 'bin_centers') && ~isempty(tuning.bin_centers)
        x_min = min(tuning.bin_centers);
        x_max = max(tuning.bin_centers);
        x_fit = linspace(x_min, x_max, 100);
    else
        x_fit = [];
        y_fit = [];
        return;
    end
    
    % Get model info
    if ~isfield(tuning, 'shuffled') || ~isfield(tuning.shuffled, 'best_model')
        y_fit = [];
        return;
    end
    
    model_type = tuning.shuffled.best_model;
    
    if isfield(tuning.shuffled, 'best_model_info')
        params = tuning.shuffled.best_model_info.beta;
    elseif isfield(tuning.shuffled, 'beta')
        params = tuning.shuffled.beta;
    else
        y_fit = [];
        return;
    end
    
    % Evaluate model
    switch model_type
        case 'linear'
            y_fit = polyval(params, x_fit);
            
        case 'quadratic'
            y_fit = polyval(params, x_fit);
            
        case 'cubic'
            y_fit = polyval(params, x_fit);
            
        case 'gaussian'
            % params = [amplitude, mu, sigma, baseline]
            if length(params) >= 4
                y_fit = params(1) * exp(-(x_fit - params(2)).^2 / (2*params(3)^2)) + params(4);
            else
                y_fit = [];
            end
            
        case 'asymmetric_gaussian'
            % params = [R_max, x_max, sigma_minus, sigma_plus]
            if length(params) >= 4
                % Use ModelSelectionTuning class if available
                try
                    y_fit = ModelSelectionTuning.asymmetric_gaussian_static(params, x_fit);
                catch
                    % Fallback implementation
                    sigma = (x_fit < params(2)) * params(3) + (x_fit >= params(2)) * params(4);
                    y_fit = params(1) * exp(-((x_fit - params(2)).^2) ./ (2 * sigma.^2));
                end
            else
                y_fit = [];
            end
            
        case 'sigmoid'
            % params = [amplitude, k, x0, baseline]
            if length(params) >= 4
                y_fit = params(1) ./ (1 + exp(-params(2) * (x_fit - params(3)))) + params(4);
            else
                y_fit = [];
            end
            
        otherwise
            y_fit = [];
    end
end
