% analyze_tf_model_distribution.m
% Script to analyze the distribution of selected models in TF tuning curves
% Counts how many instances of each model are chosen as best across all
% CSV files (one model per condition per cluster)

clear; clc;

ctl = RC2Analysis();

% Directory containing the CSV files
csv_dir = fullfile(ctl.path_config.formatted_data_dir, 'csvs', 'tf_tuning_curves');

% Get all model_comparison CSV files
csv_files = dir(fullfile(csv_dir, '*_model_comparison.csv'));

fprintf('Found %d CSV files to analyze\n\n', length(csv_files));

% Initialize storage for all best models
all_best_models = {};
all_conditions = {};
all_clusters = {};
all_probe_ids = {};
all_r_squared = [];
all_p_values = [];

% Loop through each CSV file
for i = 1:length(csv_files)
    csv_path = fullfile(csv_dir, csv_files(i).name);
    fprintf('Processing: %s\n', csv_files(i).name);
    
    % Read the CSV file
    data = readtable(csv_path);
    
    % Filter for best models (is_best_model == 1)
    best_models = data(data.is_best_model == 1, :);
    
    % Get corresponding MAT file for p-values
    mat_filename = strrep(csv_files(i).name, '_model_comparison.csv', '.mat');
    mat_path = fullfile(csv_dir, mat_filename);
    
    if exist(mat_path, 'file')
        % Load the MAT file
        mat_data = load(mat_path);
        
        % Extract p-values for each best model
        for j = 1:height(best_models)
            cluster_id = best_models.cluster_id(j);
            condition = best_models.trial_group{j};
            
            % Find matching entry in MAT file
            % Assuming the structure has fields like results.cluster_id, results.trial_group, results.p_value
            if isfield(mat_data, 'results')
                results = mat_data.results;
                % Find matching cluster and condition
                match_idx = find([results.cluster_id] == cluster_id & strcmp({results.trial_group}, condition), 1);
                if ~isempty(match_idx)
                    p_val = results(match_idx).p_value;
                else
                    p_val = NaN;
                end
            else
                p_val = NaN;
            end
            
            % Store the results
            all_best_models = [all_best_models; best_models.model_name(j)];
            all_conditions = [all_conditions; best_models.trial_group(j)];
            all_clusters = [all_clusters; num2cell(best_models.cluster_id(j))];
            all_probe_ids = [all_probe_ids; best_models.probe_id(j)];
            all_r_squared = [all_r_squared; best_models.rsq_on_mean(j)];
            all_p_values = [all_p_values; p_val];
        end
    else
        fprintf('  Warning: MAT file not found: %s\n', mat_filename);
        % Store with NaN p-values
        all_best_models = [all_best_models; best_models.model_name];
        all_conditions = [all_conditions; best_models.trial_group];
        all_clusters = [all_clusters; num2cell(best_models.cluster_id)];
        all_probe_ids = [all_probe_ids; best_models.probe_id];
        all_r_squared = [all_r_squared; best_models.rsq_on_mean];
        all_p_values = [all_p_values; nan(height(best_models), 1)];
    end
    
    fprintf('  Found %d best models\n', height(best_models));
end

fprintf('\n=== SUMMARY ===\n');
fprintf('Total number of selected models: %d\n\n', length(all_best_models));

% Count occurrences of each model
unique_models = unique(all_best_models);
model_counts = zeros(length(unique_models), 1);

for i = 1:length(unique_models)
    model_counts(i) = sum(strcmp(all_best_models, unique_models{i}));
end

% Sort by count (descending)
[model_counts_sorted, sort_idx] = sort(model_counts, 'descend');
unique_models_sorted = unique_models(sort_idx);

% Calculate percentages
model_percentages = (model_counts_sorted / length(all_best_models)) * 100;

% Display results
fprintf('Model Distribution:\n');
fprintf('%-25s %10s %10s\n', 'Model', 'Count', 'Percentage');
fprintf('%s\n', repmat('-', 1, 50));
for i = 1:length(unique_models_sorted)
    fprintf('%-25s %10d %9.2f%%\n', unique_models_sorted{i}, ...
            model_counts_sorted(i), model_percentages(i));
end

% Create a bar plot
figure('Position', [100, 100, 800, 600]);
bar(model_counts_sorted);
set(gca, 'XTickLabel', unique_models_sorted, 'XTick', 1:length(unique_models_sorted));
xtickangle(45);
ylabel('Count');
title(sprintf('Distribution of Best-Fit Models (N=%d)', length(all_best_models)));
grid on;

% Add percentage labels on bars
for i = 1:length(model_counts_sorted)
    text(i, model_counts_sorted(i), sprintf('%.1f%%', model_percentages(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 9);
end

% Also show distribution by condition
fprintf('\n=== Distribution by Condition ===\n');
unique_conditions = unique(all_conditions);
for c = 1:length(unique_conditions)
    condition_mask = strcmp(all_conditions, unique_conditions{c});
    condition_models = all_best_models(condition_mask);
    
    fprintf('\nCondition: %s (N=%d)\n', unique_conditions{c}, sum(condition_mask));
    
    condition_unique_models = unique(condition_models);
    for m = 1:length(condition_unique_models)
        count = sum(strcmp(condition_models, condition_unique_models{m}));
        pct = (count / length(condition_models)) * 100;
        fprintf('  %-20s: %4d (%5.1f%%)\n', condition_unique_models{m}, count, pct);
    end
end

% Statistical test for over-representation
% Chi-square goodness of fit test (null hypothesis: equal distribution)
fprintf('\n=== Statistical Test for Over-Representation ===\n');
expected_count = length(all_best_models) / length(unique_models);
expected_counts = repmat(expected_count, size(model_counts_sorted));
chi2stat = sum((model_counts_sorted - expected_counts).^2 ./ expected_counts);
df = length(unique_models) - 1;
p = 1 - chi2cdf(chi2stat, df);
fprintf('Chi-square test for equal distribution:\n');
fprintf('  Chi-square statistic: %.4f\n', chi2stat);
fprintf('  Degrees of freedom: %d\n', df);
fprintf('  p-value: %.4e\n', p);
if p < 0.05
    fprintf('  Result: Significant over-representation detected (reject null hypothesis)\n');
else
    fprintf('  Result: No significant over-representation (fail to reject null hypothesis)\n');
end

% Create scatterplot of R-squared vs p-value
figure('Position', [100, 100, 900, 700]);
hold on;

% Define colors for each condition
condition_colors = struct();
condition_colors.('VT') = [0.2, 0.6, 0.8];  % Blue
condition_colors.('V') = [0.9, 0.4, 0.2];   % Orange/Red
condition_colors.('T_Vstatic') = [0.3, 0.7, 0.3];  % Green

% Get unique conditions and plot each separately for legend
unique_conditions_scatter = unique(all_conditions);
legend_handles = [];

for c = 1:length(unique_conditions_scatter)
    condition = unique_conditions_scatter{c};
    condition_mask = strcmp(all_conditions, condition);
    
    % Get color for this condition
    if isfield(condition_colors, condition)
        color = condition_colors.(condition);
    else
        % Default color if condition not found
        color = [0.5, 0.5, 0.5];
    end
    
    % Plot scatter
    h = scatter(all_p_values(condition_mask), all_r_squared(condition_mask), ...
                60, color, 'filled', 'MarkerFaceAlpha', 0.6, ...
                'MarkerEdgeColor', color*0.7, 'MarkerEdgeAlpha', 0.8);
    legend_handles(c) = h;
end

% Add significance threshold line
sig_threshold = 0.05;
yline(0, 'k--', 'LineWidth', 1.5);  % R-squared = 0 line
xline(sig_threshold, 'r--', 'LineWidth', 2, 'Label', sprintf('p = %.3f', sig_threshold), ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 11);

% Labels and formatting
xlabel('p-value', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('R^2', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Model Fit Quality Across All Probes (N=%d)', length(all_best_models)), ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', 11);

% Add legend
legend(legend_handles, unique_conditions_scatter, 'Location', 'best', ...
       'FontSize', 11, 'Box', 'on');

% Set axis limits for better visualization
if ~isempty(all_p_values) && ~all(isnan(all_p_values))
    max_p = max(all_p_values(~isnan(all_p_values)));
    if max_p > 0
        xlim([0, max_p*1.05]);
    else
        xlim([0, 1]);
    end
else
    xlim([0, 1]);
end

if ~isempty(all_r_squared) && ~all(isnan(all_r_squared))
    min_r = min(all_r_squared(~isnan(all_r_squared)));
    max_r = max(all_r_squared(~isnan(all_r_squared)));
    if max_r > min_r
        ylim([min_r*1.1, max_r*1.05]);
    else
        ylim([min_r - 0.1, max_r + 0.1]);
    end
end

hold off;

% Print summary statistics for the scatter plot
fprintf('\n=== R-squared and p-value Statistics ===\n');
fprintf('R-squared: Mean = %.4f, Median = %.4f, Range = [%.4f, %.4f]\n', ...
        mean(all_r_squared), median(all_r_squared), min(all_r_squared), max(all_r_squared));
fprintf('p-value: Mean = %.4f, Median = %.4f, Range = [%.4f, %.4f]\n', ...
        mean(all_p_values), median(all_p_values), min(all_p_values), max(all_p_values));
fprintf('Significant models (p < %.3f): %d (%.1f%%)\n', ...
        sig_threshold, sum(all_p_values < sig_threshold), ...
        100*sum(all_p_values < sig_threshold)/length(all_p_values));

fprintf('\nAnalysis complete!\n');
