% Compare VT vs T_Vstatic for filtered clusters
%
% This script analyzes temporal frequency tuning differences between VT and 
% T_Vstatic conditions, but ONLY for clusters that showed specific tuning
% patterns in the V vs VT comparison:
%   - VT_only: Tuned to VT but not V
%   - both_different_model: Tuned to both V and VT with different best models
%   - both_different_params: Tuned to both V and VT with same model but different parameters
%
% The script:
% 1. Loads tf_tuning_comparison_across_conditions.csv to filter clusters
% 2. Loads tf_tuning_model_fits.csv for the raw model fit data
% 3. Compares VT vs T_Vstatic for the filtered clusters
% 4. Saves classification results to CSV
% 5. Generates visualization plots (VT left, T_Vstatic right)
%
% See also: compare_tf_tuning_across_conditions

%% Configuration

% Path to CSV files (relative to figures directory)
figure_dir = {'tf_tuning_curves'};
filter_csv_filename = 'tf_tuning_comparison_across_conditions.csv';
model_fits_csv_filename = 'tf_tuning_model_fits.csv';

% Categories to include from V vs VT comparison
filter_categories = {'VT_only', 'both_different_model', 'both_different_params', "both_similar"};

% Significance threshold
p_threshold = 0.05;

% Threshold for "significantly different" parameters (as fraction of parameter value)
param_diff_threshold = 0.1;  % 10% difference considered significant

%% Load data

ctl = RC2Analysis();

% Setup directory
ctl.setup_figures(figure_dir, true);  % Enable figure saving

% Load filter CSV (V vs VT comparison results)
filter_csv_path = fullfile(ctl.figs.curr_dir, filter_csv_filename);
if ~exist(filter_csv_path, 'file')
    error('Filter CSV file not found: %s\nPlease run compare_tf_tuning_across_conditions.m first.', filter_csv_path);
end

fprintf('Loading filter data from:\n  %s\n\n', filter_csv_path);
filter_data = readtable(filter_csv_path);

% Load model fits CSV
model_fits_csv_path = fullfile(ctl.figs.curr_dir, model_fits_csv_filename);
if ~exist(model_fits_csv_path, 'file')
    error('Model fits CSV file not found: %s\nPlease run tf_tuning_all_clusters_all_conditions.m first.', model_fits_csv_path);
end

fprintf('Loading model fit data from:\n  %s\n\n', model_fits_csv_path);
model_fits = readtable(model_fits_csv_path);

%% Filter clusters based on V vs VT comparison

fprintf('Filtering clusters...\n');

% Filter for focus clusters: V significant + V R² > 0.65
if ~ismember('is_focus_cluster', filter_data.Properties.VariableNames)
    error('is_focus_cluster column not found. Please re-run compare_tf_tuning_across_conditions.m');
end

focus_mask = logical(filter_data.is_focus_cluster);

fprintf('  Total clusters in original data: %d\n', height(filter_data));
fprintf('  Focus clusters (V significant + V R² > 0.65): %d\n', sum(focus_mask));

% Apply focus filter
filtered_clusters = filter_data(focus_mask, :);
n_filtered_clusters = height(filtered_clusters);

fprintf('  Clusters for VT vs T_Vstatic analysis: %d\n\n', n_filtered_clusters);

if n_filtered_clusters == 0
    error('No clusters found in the specified categories.');
end

%% Filter model fits to keep only VT and T_Vstatic for filtered clusters

fprintf('Filtering model fits...\n');

% Create mask for filtered clusters and conditions
model_fits_mask = false(height(model_fits), 1);
for ii = 1:n_filtered_clusters
    probe_id = filtered_clusters.probe_id{ii};
    cluster_id = filtered_clusters.cluster_id(ii);
    
    cluster_match = strcmp(model_fits.probe_id, probe_id) & ...
                    model_fits.cluster_id == cluster_id;
    condition_match = strcmp(model_fits.condition, 'VT') | ...
                      strcmp(model_fits.condition, 'T_Vstatic');
    
    model_fits_mask = model_fits_mask | (cluster_match & condition_match);
end

model_fits_filtered = model_fits(model_fits_mask, :);

fprintf('  Model fit rows after filtering: %d\n\n', height(model_fits_filtered));

%% Analyze VT vs T_Vstatic for filtered clusters

fprintf('Processing filtered clusters...\n');

% Initialize results table
results = filtered_clusters(:, {'probe_id', 'cluster_id', 'region'});
results.V_VT_category = filtered_clusters.tuning_category;  % Keep original category
results.VT_is_tuned = false(n_filtered_clusters, 1);
results.T_Vstatic_is_tuned = false(n_filtered_clusters, 1);
results.VT_p_value = nan(n_filtered_clusters, 1);
results.T_Vstatic_p_value = nan(n_filtered_clusters, 1);
results.VT_best_model = cell(n_filtered_clusters, 1);
results.T_Vstatic_best_model = cell(n_filtered_clusters, 1);
results.VT_rsq = nan(n_filtered_clusters, 1);
results.T_Vstatic_rsq = nan(n_filtered_clusters, 1);
results.model_agreement_VT_Tvstatic = false(n_filtered_clusters, 1);
results.param_diff_metric = nan(n_filtered_clusters, 1);  % Model-specific parameter difference
results.params_are_similar = false(n_filtered_clusters, 1);
results.tuning_category = cell(n_filtered_clusters, 1);

for ii = 1:n_filtered_clusters
    probe_id = results.probe_id{ii};
    cluster_id = results.cluster_id(ii);
    
    % Get rows for this cluster
    cluster_rows = strcmp(model_fits_filtered.probe_id, probe_id) & ...
                   model_fits_filtered.cluster_id == cluster_id;
    cluster_data = model_fits_filtered(cluster_rows, :);
    
    % Extract data for VT and T_Vstatic
    conditions = {'VT', 'T_Vstatic'};
    for kk = 1:length(conditions)
        cond = conditions{kk};
        cond_row = strcmp(cluster_data.condition, cond);
        
        if any(cond_row)
            row_data = cluster_data(cond_row, :);
            
            % Tuning status
            results.([cond '_is_tuned'])(ii) = row_data.is_significant;
            results.([cond '_p_value'])(ii) = row_data.p_value_shuffle;
            results.([cond '_best_model']){ii} = row_data.best_model{1};
            results.([cond '_rsq'])(ii) = row_data.rsq;
        else
            results.([cond '_best_model']){ii} = 'N/A';
        end
    end
    
    % Check model agreement between VT and T_Vstatic
    if ~strcmp(results.VT_best_model{ii}, 'N/A') && ~strcmp(results.T_Vstatic_best_model{ii}, 'N/A')
        results.model_agreement_VT_Tvstatic(ii) = strcmp(results.VT_best_model{ii}, results.T_Vstatic_best_model{ii});
        
        % If same model and both tuned, compare parameters using model-specific logic
        if results.model_agreement_VT_Tvstatic(ii) && results.VT_is_tuned(ii) && results.T_Vstatic_is_tuned(ii)
            VT_data = cluster_data(strcmp(cluster_data.condition, 'VT'), :);
            Tvstatic_data = cluster_data(strcmp(cluster_data.condition, 'T_Vstatic'), :);
            
            % Get model-specific parameters and compare
            model_type = results.VT_best_model{ii};
            [param_diff, params_similar] = compare_model_params(model_type, VT_data, Tvstatic_data, param_diff_threshold);
            
            results.param_diff_metric(ii) = param_diff;
            results.params_are_similar(ii) = params_similar;
        end
    end
end

fprintf('Done processing filtered clusters.\n\n');

%% Categorize clusters based on VT vs T_Vstatic

fprintf('Categorizing clusters based on VT vs T_Vstatic...\n');

for ii = 1:n_filtered_clusters
    VT_tuned = results.VT_is_tuned(ii);
    Tvstatic_tuned = results.T_Vstatic_is_tuned(ii);
    
    if ~VT_tuned && ~Tvstatic_tuned
        results.tuning_category{ii} = 'neither';
    elseif VT_tuned && ~Tvstatic_tuned
        results.tuning_category{ii} = 'VT_only';
    elseif ~VT_tuned && Tvstatic_tuned
        results.tuning_category{ii} = 'T_Vstatic_only';
    elseif VT_tuned && Tvstatic_tuned
        % Both tuned - check if similar
        if ~results.model_agreement_VT_Tvstatic(ii)
            results.tuning_category{ii} = 'both_different_model';
        elseif results.params_are_similar(ii)
            results.tuning_category{ii} = 'both_similar';
        else
            results.tuning_category{ii} = 'both_different_params';
        end
    end
end

fprintf('Done categorizing.\n\n');

%% Print summary statistics

fprintf('=== SUMMARY STATISTICS (VT vs T_Vstatic for filtered clusters) ===\n\n');

fprintf('All focus clusters (V significant + V R² > 0.65): %d\n\n', n_filtered_clusters);

fprintf('Filtered cluster categories (from V vs VT comparison):\n');
v_vt_categories = unique(results.V_VT_category);
for ii = 1:length(v_vt_categories)
    cat = v_vt_categories{ii};
    n_cat = sum(strcmp(results.V_VT_category, cat));
    fprintf('  %s: %d (%.1f%%)\n', cat, n_cat, 100*n_cat/n_filtered_clusters);
end
fprintf('\n');

fprintf('Question 1: How many filtered clusters are tuned to VT?\n');
n_VT_tuned = sum(results.VT_is_tuned);
n_VT_not_tuned = sum(~results.VT_is_tuned);
fprintf('  Tuned to VT: %d (%.1f%%)\n', n_VT_tuned, 100*n_VT_tuned/n_filtered_clusters);
fprintf('  NOT tuned to VT: %d (%.1f%%)\n\n', n_VT_not_tuned, 100*n_VT_not_tuned/n_filtered_clusters);

fprintf('Question 2: Of those NOT tuned in VT, how many are tuned to T_Vstatic?\n');
not_VT_idx = ~results.VT_is_tuned;
n_not_VT_but_Tvstatic = sum(not_VT_idx & results.T_Vstatic_is_tuned);
fprintf('  NOT VT but tuned to T_Vstatic: %d (%.1f%% of non-VT-tuned)\n\n', ...
    n_not_VT_but_Tvstatic, 100*n_not_VT_but_Tvstatic/max(n_VT_not_tuned,1));

fprintf('Question 3: Of those tuned to VT, how many are also tuned to T_Vstatic?\n');
VT_idx = results.VT_is_tuned;
n_VT_and_Tvstatic = sum(VT_idx & results.T_Vstatic_is_tuned);
fprintf('  Tuned to both VT and T_Vstatic: %d (%.1f%% of VT-tuned)\n', ...
    n_VT_and_Tvstatic, 100*n_VT_and_Tvstatic/max(n_VT_tuned,1));

fprintf('  Of those tuned to both:\n');
both_idx = VT_idx & results.T_Vstatic_is_tuned;
n_same_model = sum(both_idx & results.model_agreement_VT_Tvstatic);
n_diff_model = sum(both_idx & ~results.model_agreement_VT_Tvstatic);
fprintf('    Same model: %d (%.1f%%)\n', n_same_model, 100*n_same_model/max(n_VT_and_Tvstatic,1));
fprintf('    Different model: %d (%.1f%%)\n', n_diff_model, 100*n_diff_model/max(n_VT_and_Tvstatic,1));

same_model_idx = both_idx & results.model_agreement_VT_Tvstatic;
n_similar_params = sum(same_model_idx & results.params_are_similar);
n_diff_params = sum(same_model_idx & ~results.params_are_similar);
fprintf('    Same model, similar parameters: %d (%.1f%% of same-model)\n', ...
    n_similar_params, 100*n_similar_params/max(n_same_model,1));
fprintf('    Same model, different parameters: %d (%.1f%% of same-model)\n\n', ...
    n_diff_params, 100*n_diff_params/max(n_same_model,1));

fprintf('Overall distribution by VT vs T_Vstatic category:\n');
categories = unique(results.tuning_category);
for ii = 1:length(categories)
    cat = categories{ii};
    n_cat = sum(strcmp(results.tuning_category, cat));
    fprintf('  %s: %d (%.1f%%)\n', cat, n_cat, 100*n_cat/n_filtered_clusters);
end
fprintf('\n');

fprintf('Cross-tabulation: V vs VT categories within each VT vs T_Vstatic category:\n');
for ii = 1:length(categories)
    cat = categories{ii};
    cat_mask = strcmp(results.tuning_category, cat);
    n_cat = sum(cat_mask);
    
    if n_cat == 0
        continue;
    end
    
    fprintf('  %s (n=%d):\n', cat, n_cat);
    
    % Break down by V vs VT category
    for jj = 1:length(v_vt_categories)
        v_vt_cat = v_vt_categories{jj};
        n_in_both = sum(cat_mask & strcmp(results.V_VT_category, v_vt_cat));
        if n_in_both > 0
            fprintf('    - from %s: %d (%.1f%%)\n', v_vt_cat, n_in_both, 100*n_in_both/n_cat);
        end
    end
    fprintf('\n');
end

%% Combined categorization across both comparisons

fprintf('=== COMBINED CATEGORIZATION (3 Main Categories) ===\n\n');

% Category 1: V = VT (same model)
cat1_mask = strcmp(results.V_VT_category, 'both_similar') | strcmp(results.V_VT_category, 'both_different_params');
n_cat1 = sum(cat1_mask);
n_cat1_similar = sum(strcmp(results.V_VT_category, 'both_similar'));
n_cat1_different = sum(strcmp(results.V_VT_category, 'both_different_params'));

fprintf('Category 1: V = VT (same model)\n');
fprintf('  Requirements: V sig + V R²>0.65, same model as VT\n');
fprintf('  Total: %d (%.1f%% of focus clusters)\n', n_cat1, 100*n_cat1/n_filtered_clusters);
fprintf('    Similar parameters: %d (%.1f%%)\n', n_cat1_similar, 100*n_cat1_similar/max(n_cat1,1));
fprintf('    Different parameters: %d (%.1f%%)\n\n', n_cat1_different, 100*n_cat1_different/max(n_cat1,1));

% Category 2: V ≠ VT, but VT = T_Vstatic (same model)
cat2_mask = ~cat1_mask & results.model_agreement_VT_Tvstatic;
n_cat2 = sum(cat2_mask);

% For those with same model, check if both are tuned to compare parameters
cat2_both_tuned = cat2_mask & results.VT_is_tuned & results.T_Vstatic_is_tuned;
n_cat2_similar = sum(cat2_both_tuned & results.params_are_similar);
n_cat2_different = sum(cat2_both_tuned & ~results.params_are_similar);

fprintf('Category 2: V ≠ VT, but VT = T_Vstatic (same model)\n');
fprintf('  Requirements: V sig + V R²>0.65, different model from VT, same model for VT and T_Vstatic\n');
fprintf('  Total: %d (%.1f%% of focus clusters)\n', n_cat2, 100*n_cat2/n_filtered_clusters);
fprintf('    Of these, both VT and T_Vstatic significant: %d (%.1f%%)\n', sum(cat2_both_tuned), 100*sum(cat2_both_tuned)/max(n_cat2,1));
fprintf('      Similar parameters: %d (%.1f%% of both-tuned)\n', n_cat2_similar, 100*n_cat2_similar/max(sum(cat2_both_tuned),1));
fprintf('      Different parameters: %d (%.1f%% of both-tuned)\n\n', n_cat2_different, 100*n_cat2_different/max(sum(cat2_both_tuned),1));

% Category 3: V ≠ VT, and VT ≠ T_Vstatic (different models)
cat3_mask = ~cat1_mask & ~results.model_agreement_VT_Tvstatic;
n_cat3 = sum(cat3_mask);

fprintf('Category 3: V ≠ VT, and VT ≠ T_Vstatic (different models)\n');
fprintf('  Requirements: V sig + V R²>0.65, different model from VT, different model for VT and T_Vstatic\n');
fprintf('  Total: %d (%.1f%% of focus clusters)\n\n', n_cat3, 100*n_cat3/n_filtered_clusters);

% Verification
fprintf('Verification: %d + %d + %d = %d (should equal %d)\n\n', ...
    n_cat1, n_cat2, n_cat3, n_cat1+n_cat2+n_cat3, n_filtered_clusters);

%% Save results to CSV

output_filename = fullfile(ctl.figs.curr_dir, 'tf_tuning_vt_vs_tvstatic_filtered.csv');
writetable(results, output_filename);

fprintf('Results saved to:\n  %s\n\n', output_filename);

%% Visualize fitted curves by category

fprintf('Creating category visualization figures...\n');
fprintf('Figures will be saved to: %s\n\n', ctl.figs.curr_dir);

% Get unique probe_ids to load tuning data
unique_probes = unique(results.probe_id);

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

fprintf('Loaded tuning data for %d probes.\n', length(unique_probes));

%% Plot detailed categories (V vs VT × VT vs T_Vstatic breakdown)

fprintf('\nCreating detailed category visualization figures...\n');

% Define VT vs T_Vstatic categories to plot
vt_tvstatic_categories = {'neither', 'VT_only', 'T_Vstatic_only', 'both_similar', 'both_different_params', 'both_different_model'};
vt_tvstatic_titles = {'Neither VT nor T_Vstatic', 'VT Only', 'T_Vstatic Only', ...
                      'Both Similar', 'Both (Different Parameters)', 'Both (Different Model)'};

% Get unique V vs VT categories present in the data
v_vt_categories = unique(results.V_VT_category);

% Plot each combination of V vs VT category and VT vs T_Vstatic category
for v_vt_cat_idx = 1:length(v_vt_categories)
    v_vt_cat = v_vt_categories{v_vt_cat_idx};
    
    fprintf('  Processing V vs VT category: %s\n', v_vt_cat);
    
    for cat_idx = 1:length(vt_tvstatic_categories)
        cat = vt_tvstatic_categories{cat_idx};
        cat_title = vt_tvstatic_titles{cat_idx};
        
        % Get clusters in both this VT vs T_Vstatic category AND the current V vs VT category
        cat_clusters = strcmp(results.tuning_category, cat) & strcmp(results.V_VT_category, v_vt_cat);
        n_cat_clusters = sum(cat_clusters);
        
        fprintf('    Plotting VT vs T_Vstatic category: %s (%d clusters)\n', cat, n_cat_clusters);
        
        if n_cat_clusters == 0
            fprintf('      Skipping - no clusters in this category.\n');
            continue;
        end
    
        % Generate distinct colors for each cluster in this category
        cluster_indices = find(cat_clusters);
        colors = lines(n_cat_clusters);
        color_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
        for idx = 1:n_cat_clusters
            color_map(cluster_indices(idx)) = colors(idx, :);
        end
        
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
        xlabel('Equivalent Bins');
        ylabel('Firing Rate (Hz)');
        
        % Track x-axis range for axis alignment
        all_x_min = inf;
        all_x_max = -inf;
        
        % Plot individual cluster curves with unique colors
        for ii = 1:n_filtered_clusters
            if ~cat_clusters(ii)
                continue;
            end
            
            probe_id = results.probe_id{ii};
            cluster_id = results.cluster_id(ii);
            
            probe_key = matlab.lang.makeValidName(probe_id);
            cluster_key = sprintf('cluster_%d', cluster_id);
            
            % Check if data exists
            if ~isfield(tuning_lookup, probe_key) || ...
               ~isfield(tuning_lookup.(probe_key), cluster_key)
                continue;
            end
            
            cluster_data = tuning_lookup.(probe_key).(cluster_key);
            cluster_color = color_map(ii);
            
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
        
        % Add title with both categories and cluster count
        sgtitle(sprintf('V vs VT: %s | VT vs T\\_Vstatic: %s | n = %d clusters', ...
            strrep(v_vt_cat, '_', '\\_'), strrep(cat_title, '_', '\\_'), n_cat_clusters), ...
            'FontSize', 12, 'FontWeight', 'bold');
        
        % Save figure (delete existing file first to allow overwriting)
        fig_filename = sprintf('tf_tuning_detailed_V_vt_%s_VT_tvstatic_%s.pdf', v_vt_cat, cat);
        full_path = fullfile(ctl.figs.curr_dir, fig_filename);
        if exist(full_path, 'file')
            delete(full_path);
        end
        ctl.figs.save_fig(fig_filename);
        fprintf('      Saved: %s\n', full_path);
    end
    
    fprintf('\n');
end

fprintf('Detailed category figures saved.\n\n');

%% Plot combined categories (3 main categories)

fprintf('Creating combined category visualization figures...\n');

% Define combined categories to plot (3 main categories)
combined_categories = {'cat1_V_eq_VT', 'cat2_V_neq_VT_VT_eq_Tvs', 'cat3_V_neq_VT_VT_neq_Tvs'};
combined_titles = {'Category 1: V = VT', 'Category 2: V ≠ VT, VT = T\_Vstatic', 'Category 3: V ≠ VT, VT ≠ T\_Vstatic'};

% Plot each combined category
for cat_idx = 1:length(combined_categories)
    cat = combined_categories{cat_idx};
    cat_title = combined_titles{cat_idx};
    
    % Create category mask
    if cat_idx == 1
        % Category 1: V = VT
        cat_mask = strcmp(results.V_VT_category, 'both_similar') | strcmp(results.V_VT_category, 'both_different_params');
    elseif cat_idx == 2
        % Category 2: V ≠ VT, but VT = T_Vstatic (same model)
        cat1_mask = strcmp(results.V_VT_category, 'both_similar') | strcmp(results.V_VT_category, 'both_different_params');
        cat_mask = ~cat1_mask & results.model_agreement_VT_Tvstatic;
    else
        % Category 3: V ≠ VT, and VT ≠ T_Vstatic (different models)
        cat1_mask = strcmp(results.V_VT_category, 'both_similar') | strcmp(results.V_VT_category, 'both_different_params');
        cat_mask = ~cat1_mask & ~results.model_agreement_VT_Tvstatic;
    end
    
    cat_clusters = cat_mask;
    n_cat_clusters = sum(cat_clusters);
        
    fprintf('  Plotting %s (%d clusters)\n', cat_title, n_cat_clusters);
        
        if n_cat_clusters == 0
            fprintf('    Skipping - no clusters in this category.\n');
            continue;
        end
    
        % Generate distinct colors for each cluster in this category
        cluster_indices = find(cat_clusters);
        colors = lines(n_cat_clusters);
        color_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
        for idx = 1:n_cat_clusters
            color_map(cluster_indices(idx)) = colors(idx, :);
        end
        
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
        xlabel('Equivalent Bins');
        ylabel('Firing Rate (Hz)');
        
        % Track x-axis range for axis alignment
        all_x_min = inf;
        all_x_max = -inf;
        
        % Plot individual cluster curves with unique colors
        for ii = 1:n_filtered_clusters
            if ~cat_clusters(ii)
                continue;
            end
            
            probe_id = results.probe_id{ii};
            cluster_id = results.cluster_id(ii);
            
            probe_key = matlab.lang.makeValidName(probe_id);
            cluster_key = sprintf('cluster_%d', cluster_id);
            
            % Check if data exists
            if ~isfield(tuning_lookup, probe_key) || ...
               ~isfield(tuning_lookup.(probe_key), cluster_key)
                continue;
            end
            
            cluster_data = tuning_lookup.(probe_key).(cluster_key);
            cluster_color = color_map(ii);
            
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
        sgtitle(sprintf('%s | n = %d clusters', cat_title, n_cat_clusters), ...
            'FontSize', 12, 'FontWeight', 'bold');
        
        % Save figure (delete existing file first to allow overwriting)
        fig_filename = sprintf('tf_tuning_focus_%s.pdf', cat);
        full_path = fullfile(ctl.figs.curr_dir, fig_filename);
        if exist(full_path, 'file')
            delete(full_path);
        end
        ctl.figs.save_fig(fig_filename);
        fprintf('    Saved: %s\n', full_path);
end

fprintf('Category figures saved.\n\n');

fprintf('=== Analysis complete ===\n');
fprintf('Output CSV: tf_tuning_vt_vs_tvstatic_filtered.csv\n');
fprintf('Detailed category figures: tf_tuning_detailed_V_vt_*_VT_tvstatic_*.pdf\n');
fprintf('Combined category figures: tf_tuning_focus_cat*.pdf\n');


%% Helper function to evaluate fitted curve

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


function [param_diff_metric, params_similar] = compare_model_params(model_type, data1, data2, threshold)
%%compare_model_params Compare model parameters based on model type
%
%   [DIFF, SIMILAR] = compare_model_params(MODEL_TYPE, DATA1, DATA2, THRESHOLD)
%   Returns a metric of parameter difference and whether parameters are similar.

    switch model_type
        case {'linear', 'quadratic', 'cubic'}
            % For polynomials, check if we have individual coefficient columns
            % If not available, we can't do a meaningful comparison
            % (the CSV may not store polynomial coefficients separately)
            
            % Try to compare R² values as a proxy
            if ~isnan(data1.rsq) && ~isnan(data2.rsq)
                param_diff_metric = abs(data1.rsq - data2.rsq);
                params_similar = param_diff_metric < threshold;
            else
                param_diff_metric = nan;
                params_similar = false;
            end
            
        case 'gaussian'
            % Gaussian: [amplitude, mu, sigma, baseline]
            if ~ismember('param_amplitude', data1.Properties.VariableNames) || ...
               ~ismember('param_center', data1.Properties.VariableNames) || ...
               ~ismember('param_sigma_or_sigma_minus', data1.Properties.VariableNames)
                param_diff_metric = nan;
                params_similar = false;
                return;
            end
            
            params1 = [data1.param_amplitude, data1.param_center, data1.param_sigma_or_sigma_minus];
            params2 = [data2.param_amplitude, data2.param_center, data2.param_sigma_or_sigma_minus];
            
            if any(isnan(params1)) || any(isnan(params2))
                param_diff_metric = nan;
                params_similar = false;
                return;
            end
            
            % Compare key parameters: center, amplitude, width
            center_diff = abs(params1(2) - params2(2)) / max(abs([params1(2), params2(2)]));
            amp_diff = abs(params1(1) - params2(1)) / max(abs([params1(1), params2(1)]));
            width_diff = abs(params1(3) - params2(3)) / max(abs([params1(3), params2(3)]));
            
            % Use max difference as metric
            param_diff_metric = max([center_diff, amp_diff, width_diff]);
            params_similar = all([center_diff, amp_diff, width_diff] < threshold);
            
        case 'asymmetric_gaussian'
            % Asymmetric Gaussian: [R_max, x_max, sigma_minus, sigma_plus]
            if ~ismember('param_amplitude', data1.Properties.VariableNames) || ...
               ~ismember('param_center', data1.Properties.VariableNames) || ...
               ~ismember('param_sigma_or_sigma_minus', data1.Properties.VariableNames) || ...
               ~ismember('param_sigma_plus', data1.Properties.VariableNames)
                param_diff_metric = nan;
                params_similar = false;
                return;
            end
            
            params1 = [data1.param_amplitude, data1.param_center, data1.param_sigma_or_sigma_minus, data1.param_sigma_plus];
            params2 = [data2.param_amplitude, data2.param_center, data2.param_sigma_or_sigma_minus, data2.param_sigma_plus];
            
            if any(isnan(params1)) || any(isnan(params2))
                param_diff_metric = nan;
                params_similar = false;
                return;
            end
            
            % Compare key parameters
            center_diff = abs(params1(2) - params2(2)) / max(abs([params1(2), params2(2)]));
            amp_diff = abs(params1(1) - params2(1)) / max(abs([params1(1), params2(1)]));
            sigma_minus_diff = abs(params1(3) - params2(3)) / max(abs([params1(3), params2(3)]));
            sigma_plus_diff = abs(params1(4) - params2(4)) / max(abs([params1(4), params2(4)]));
            
            % Use max difference as metric
            param_diff_metric = max([center_diff, amp_diff, sigma_minus_diff, sigma_plus_diff]);
            params_similar = all([center_diff, amp_diff, sigma_minus_diff, sigma_plus_diff] < threshold);
            
        case 'sigmoid'
            % Sigmoid: [amplitude, k, x0, baseline]
            if ~ismember('param_amplitude', data1.Properties.VariableNames) || ...
               ~ismember('param_steepness', data1.Properties.VariableNames) || ...
               ~ismember('param_center', data1.Properties.VariableNames)
                param_diff_metric = nan;
                params_similar = false;
                return;
            end
            
            params1 = [data1.param_amplitude, data1.param_steepness, data1.param_center];
            params2 = [data2.param_amplitude, data2.param_steepness, data2.param_center];
            
            if any(isnan(params1)) || any(isnan(params2))
                param_diff_metric = nan;
                params_similar = false;
                return;
            end
            
            % Compare key parameters
            center_diff = abs(params1(3) - params2(3)) / max(abs([params1(3), params2(3)]));
            amp_diff = abs(params1(1) - params2(1)) / max(abs([params1(1), params2(1)]));
            steep_diff = abs(params1(2) - params2(2)) / max(abs([params1(2), params2(2)]));
            
            % Use max difference as metric
            param_diff_metric = max([center_diff, amp_diff, steep_diff]);
            params_similar = all([center_diff, amp_diff, steep_diff] < threshold);
            
        otherwise
            param_diff_metric = nan;
            params_similar = false;
    end
end
