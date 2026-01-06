% Plot pooled spatial firing rate heatmaps across all probes
% Groups clusters by 'a' and 'b' conditions (extracted from probe name before _rec1)
%
% This script:
%   1. Loads cached spatial analysis data from all probes
%   2. Extracts the condition ('a' or 'b') from each probe ID
%   3. Pools clusters across all probes of the same condition
%   4. Creates heatmap plots similar to the single-probe version but with pooled data
%
% REQUIRES:
%   - Cached spatial analysis files must exist in the figure directory
%   - Cache files follow naming: <probe_id>_spatial_analysis_cache.mat
%
% OUTPUT:
%   - Pooled heatmaps showing spatial firing patterns for all clusters
%     grouped by condition 'a' and 'b'

%%
% Prevent figures from popping up on screen
set(groot, 'DefaultFigureVisible', 'off');

% Configuration
experiment_groups        = {'ambient_light'};
save_figs                = true;
overwrite                = true;
figure_dir               = {'spatial_firing_rate_pooled', 'ambient_light'};

% Analysis parameters (must match those used in spatial_firing_rate_profile.m)
bin_size_cm = 2;
gauss_sigma_cm = 8;
position_threshold_cm = 90;

% Statistical analysis parameters (for determining significance)
minSpikes = 25;
minPeakRate = 1.0;
fieldFrac = 0.7;
minFieldBins = 5;
maxNumFields = 1;
nShuf = 100;
pThresh = 0.05;

% Initialize controller and get probe IDs
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));

% Determine the base cache directory (where individual probe caches are saved)
% This should be the same directory where spatial_firing_rate_profile.m saves its caches
base_cache_dir = fullfile(ctl.path_config.figure_dir, 'spatial_firing_rate', experiment_groups{:});

if ~isfolder(base_cache_dir)
    error('Cache directory not found: %s\nPlease run spatial_firing_rate_profile.m first.', base_cache_dir);
end

fprintf('Loading cached data from: %s\n', base_cache_dir);

% Initialize structures to store pooled data
pooled_data = struct();
pooled_data.a.cluster_ids = [];
pooled_data.a.probe_ids = {};
pooled_data.a.rate_smooth_long = [];
pooled_data.a.rate_smooth_short = [];
pooled_data.a.Q2_smooth_long = [];
pooled_data.a.Q2_smooth_short = [];
pooled_data.a.is_significant_long = [];
pooled_data.a.is_significant_short = [];

pooled_data.b.cluster_ids = [];
pooled_data.b.probe_ids = {};
pooled_data.b.rate_smooth_long = [];
pooled_data.b.rate_smooth_short = [];
pooled_data.b.Q2_smooth_long = [];
pooled_data.b.Q2_smooth_short = [];
pooled_data.b.is_significant_long = [];
pooled_data.b.is_significant_short = [];

% Load data from each probe and pool by condition
for pid = 1:length(probe_ids)
    probe_id = probe_ids{pid};
    fprintf('  Processing probe %d/%d: %s\n', pid, length(probe_ids), probe_id);
    
    % Extract condition ('a' or 'b') from probe name
    % Expected format: <animal_id><condition>_rec1
    % e.g., CAA-1123432a_rec1 -> condition = 'a'
    match = regexp(probe_id, '([ab])_rec1', 'tokens');
    
    if isempty(match)
        fprintf('    Warning: Could not extract condition from probe ID %s, skipping\n', probe_id);
        continue;
    end
    
    condition = match{1}{1};  % Extract 'a' or 'b'
    fprintf('    Detected condition: %s\n', condition);
    
    % Load cache file
    cache_filename = sprintf('%s_spatial_analysis_cache.mat', probe_id);
    cache_filepath = fullfile(base_cache_dir, cache_filename);
    
    if ~exist(cache_filepath, 'file')
        fprintf('    Warning: Cache file not found: %s, skipping\n', cache_filename);
        continue;
    end
    
    fprintf('    Loading cache: %s\n', cache_filename);
    cached = load(cache_filepath);
    
    % Extract relevant data
    n_clusters = cached.n_clusters;
    cluster_ids = cached.cluster_ids;
    
    % Get rate data
    rate_smooth_long = cached.all_Q2_rate_smooth_by_group.long;  % Use median (Q2)
    rate_smooth_short = cached.all_Q2_rate_smooth_by_group.short;
    
    % Determine significance for each cluster in each condition
    is_significant_long = false(n_clusters, 1);
    is_significant_short = false(n_clusters, 1);
    
    if isfield(cached, 'spatial_tuning_stats') && ~isempty(cached.spatial_tuning_stats)
        stats = cached.spatial_tuning_stats;
        for c = 1:n_clusters
            cluster_field = sprintf('cluster_%d', cluster_ids(c));
            if isfield(stats, cluster_field)
                cluster_stats = stats.(cluster_field);
                if isfield(cluster_stats, 'long') && isfield(cluster_stats.long, 'isSpatiallyTuned')
                    is_significant_long(c) = cluster_stats.long.isSpatiallyTuned;
                end
                if isfield(cluster_stats, 'short') && isfield(cluster_stats.short, 'isSpatiallyTuned')
                    is_significant_short(c) = cluster_stats.short.isSpatiallyTuned;
                end
            end
        end
    end
    
    % Add to pooled data structure
    pooled_data.(condition).cluster_ids = [pooled_data.(condition).cluster_ids; cluster_ids'];
    pooled_data.(condition).probe_ids = [pooled_data.(condition).probe_ids; repmat({probe_id}, n_clusters, 1)];
    pooled_data.(condition).rate_smooth_long = [pooled_data.(condition).rate_smooth_long; rate_smooth_long];
    pooled_data.(condition).rate_smooth_short = [pooled_data.(condition).rate_smooth_short; rate_smooth_short];
    pooled_data.(condition).Q2_smooth_long = [pooled_data.(condition).Q2_smooth_long; rate_smooth_long];
    pooled_data.(condition).Q2_smooth_short = [pooled_data.(condition).Q2_smooth_short; rate_smooth_short];
    pooled_data.(condition).is_significant_long = [pooled_data.(condition).is_significant_long; is_significant_long];
    pooled_data.(condition).is_significant_short = [pooled_data.(condition).is_significant_short; is_significant_short];
    
    % Store rate_per_trial data for distribution comparison
    if ~isfield(pooled_data.(condition), 'rate_per_trial_long')
        pooled_data.(condition).rate_per_trial_long = {};
        pooled_data.(condition).rate_per_trial_short = {};
    end
    % Ensure we're appending column cell arrays
    rate_per_trial_long_to_add = cached.all_rate_per_trial_by_group.long;
    rate_per_trial_short_to_add = cached.all_rate_per_trial_by_group.short;
    if size(rate_per_trial_long_to_add, 2) > 1
        rate_per_trial_long_to_add = rate_per_trial_long_to_add';
        rate_per_trial_short_to_add = rate_per_trial_short_to_add';
    end
    pooled_data.(condition).rate_per_trial_long = [pooled_data.(condition).rate_per_trial_long; rate_per_trial_long_to_add];
    pooled_data.(condition).rate_per_trial_short = [pooled_data.(condition).rate_per_trial_short; rate_per_trial_short_to_add];
    
    % Load TTG data if available
    ttg_cache_filename = sprintf('%s_ttg_cache.mat', probe_id);
    ttg_cache_filepath = fullfile(base_cache_dir, ttg_cache_filename);
    if exist(ttg_cache_filepath, 'file')
        ttg_cached = load(ttg_cache_filepath);
        if ~isfield(pooled_data.(condition), 'ttg_rates_long')
            pooled_data.(condition).ttg_rates_long = [];
            pooled_data.(condition).ttg_rates_short = [];
            pooled_data.(condition).ttg_norm_bin_centers = ttg_cached.ttg_norm_bin_centers;
        end
        pooled_data.(condition).ttg_rates_long = [pooled_data.(condition).ttg_rates_long; ttg_cached.all_ttg_rates.long];
        pooled_data.(condition).ttg_rates_short = [pooled_data.(condition).ttg_rates_short; ttg_cached.all_ttg_rates.short];
    end
    
    fprintf('    Added %d clusters to condition %s pool\n', n_clusters, condition);
end

% Print summary
fprintf('\n=== Pooled Data Summary ===\n');
fprintf('Region A: %d clusters from %d probes\n', ...
    length(pooled_data.a.cluster_ids), length(unique(pooled_data.a.probe_ids)));
fprintf('  - %d significant in long trials\n', sum(pooled_data.a.is_significant_long));
fprintf('  - %d significant in short trials\n', sum(pooled_data.a.is_significant_short));
fprintf('Region B: %d clusters from %d probes\n', ...
    length(pooled_data.b.cluster_ids), length(unique(pooled_data.b.probe_ids)));
fprintf('  - %d significant in long trials\n', sum(pooled_data.b.is_significant_long));
fprintf('  - %d significant in short trials\n', sum(pooled_data.b.is_significant_short));

% Get bin centers from one of the cached files (they should all be the same)
% Reload first cache to get bin centers
first_cache = load(fullfile(base_cache_dir, sprintf('%s_spatial_analysis_cache.mat', probe_ids{1})));
bin_centers_long = first_cache.all_bin_centers_by_group.long;
bin_centers_short = first_cache.all_bin_centers_by_group.short;

% Load distribution comparisons (KW tests) from saved spatial tuning fits files
fprintf('\n=== Loading distribution comparisons (KW test) from saved files ===\n');
conditions_temp = {'a', 'b'};
for cond_idx = 1:2
    condition = conditions_temp{cond_idx};
    if isempty(pooled_data.(condition).cluster_ids)
        continue;
    end
    
    n_clusters_cond = length(pooled_data.(condition).cluster_ids);
    pooled_data.(condition).dist_comparison = cell(n_clusters_cond, 1);
    pooled_data.(condition).probe_ids_for_clusters = pooled_data.(condition).probe_ids;
    
    fprintf('  %s: Loading KW results for %d clusters...\n', upper(condition), n_clusters_cond);
    
    % For each cluster, load its KW result from the corresponding probe's fits file
    for c = 1:n_clusters_cond
        probe_id = pooled_data.(condition).probe_ids{c};
        cluster_id = pooled_data.(condition).cluster_ids(c);
        
        % Load the spatial tuning fits file for this probe
        fits_filename = sprintf('%s_spatial_tuning_fits.mat', probe_id);
        fits_filepath = fullfile(base_cache_dir, fits_filename);
        
        if ~exist(fits_filepath, 'file')
            fprintf('    Warning: Fits file not found for probe %s, skipping cluster %d\n', probe_id, cluster_id);
            pooled_data.(condition).dist_comparison{c} = [];
            continue;
        end
        
        % Load and extract the KW result for this specific cluster
        try
            fits_data = load(fits_filepath);
            if isfield(fits_data, 'gaussian_fits') && isfield(fits_data.gaussian_fits, 'dist_comparison_results')
                dist_results = fits_data.gaussian_fits.dist_comparison_results;
                
                % Find the index of this cluster in the original probe data
                cluster_idx = find(fits_data.gaussian_fits.cluster_ids == cluster_id, 1);
                
                if ~isempty(cluster_idx) && cluster_idx <= length(dist_results)
                    pooled_data.(condition).dist_comparison{c} = dist_results{cluster_idx};
                else
                    pooled_data.(condition).dist_comparison{c} = [];
                end
            else
                fprintf('    Warning: No KW results in fits file for probe %s\n', probe_id);
                pooled_data.(condition).dist_comparison{c} = [];
            end
        catch ME
            fprintf('    Warning: Could not load KW result for cluster %d from probe %s: %s\n', cluster_id, probe_id, ME.message);
            pooled_data.(condition).dist_comparison{c} = [];
        end
    end
    
    % Count categorizations
    n_same_abs = 0;
    n_diff_abs = 0;
    n_same_rel = 0;
    n_diff_rel = 0;
    for c = 1:n_clusters_cond
        if ~isempty(pooled_data.(condition).dist_comparison{c})
            res = pooled_data.(condition).dist_comparison{c};
            if isfield(res, 'absolute') && ~isnan(res.absolute.same_distribution)
                if res.absolute.same_distribution
                    n_same_abs = n_same_abs + 1;
                else
                    n_diff_abs = n_diff_abs + 1;
                end
            end
            if isfield(res, 'relative') && ~isnan(res.relative.same_distribution)
                if res.relative.same_distribution
                    n_same_rel = n_same_rel + 1;
                else
                    n_diff_rel = n_diff_rel + 1;
                end
            end
        end
    end
    fprintf('    Absolute position: %d same, %d different\n', n_same_abs, n_diff_abs);
    fprintf('    Relative position: %d same, %d different\n', n_same_rel, n_diff_rel);
end

% Create pooled heatmaps for each region
conditions = {'a', 'b'};
condition_labels = {'Region A', 'Region B'};

for cond_idx = 1:length(conditions)
    condition = conditions{cond_idx};
    condition_label = condition_labels{cond_idx};
    
    if isempty(pooled_data.(condition).cluster_ids)
        fprintf('  No data for condition %s, skipping\n', condition);
        continue;
    end
    
    fprintf('\nCreating pooled heatmap for %s...\n', condition_label);
    
    % Get data for this condition
    rate_mat_long = pooled_data.(condition).Q2_smooth_long;
    rate_mat_short = pooled_data.(condition).Q2_smooth_short;
    is_sig_long = pooled_data.(condition).is_significant_long;
    is_sig_short = pooled_data.(condition).is_significant_short;
    n_clusters_cond = size(rate_mat_long, 1);
    
    % Sort clusters by peak position in long trials
    peak_positions = zeros(n_clusters_cond, 1);
    for c = 1:n_clusters_cond
        [~, max_idx] = max(rate_mat_long(c, :));
        peak_positions(c) = bin_centers_long(max_idx);
    end
    [~, sort_idx] = sort(peak_positions);
    
    % Apply sorting
    rate_mat_long = rate_mat_long(sort_idx, :);
    rate_mat_short = rate_mat_short(sort_idx, :);
    is_sig_long = is_sig_long(sort_idx);
    is_sig_short = is_sig_short(sort_idx);
    
    % Create figure with 4 subplots (2 rows × 2 columns)
    fig_pooled = ctl.figs.a4figure('landscape');
    
    % Define group info
    group_names = {'long', 'short'};
    group_labels = {'Long (0-120 cm)', 'Short (60-120 cm)'};
    
    % IMPORTANT: All clusters are already sorted by peak position in long trials
    % (from sort_idx applied earlier). We need to maintain this global order
    % when splitting into significant/non-significant groups.
    
    % For LONG trials: get the indices of sig/nonsig clusters in the current sorted order
    is_sig_long_logical = logical(is_sig_long(:));
    sig_indices_long = find(is_sig_long_logical);
    nonsig_indices_long = find(~is_sig_long_logical);
    
    % For SHORT trials: get the indices of sig/nonsig clusters in the current sorted order
    is_sig_short_logical = logical(is_sig_short(:));
    sig_indices_short = find(is_sig_short_logical);
    nonsig_indices_short = find(~is_sig_short_logical);
    
    % Now plot both long and short with consistent ordering
    for g = 1:2
        if g == 1
            rate_mat = rate_mat_long;
            sig_indices = sig_indices_long;
            nonsig_indices = nonsig_indices_long;
            bin_centers = bin_centers_long;
        else
            rate_mat = rate_mat_short;
            sig_indices = sig_indices_short;
            nonsig_indices = nonsig_indices_short;
            bin_centers = bin_centers_short;
        end
        
        % Min-max normalization per cluster
        min_vals = min(rate_mat, [], 2);
        max_vals = max(rate_mat, [], 2);
        norm_rate_mat = (rate_mat - min_vals) ./ (max_vals - min_vals);
        norm_rate_mat(isnan(norm_rate_mat) | isinf(norm_rate_mat)) = 0;
        
        % Extract significant and non-significant clusters while maintaining global sort order
        % The indices already reflect the sorted order from long trials
        sig_rate_mat = norm_rate_mat(sig_indices, :);
        nonsig_rate_mat = norm_rate_mat(nonsig_indices, :);
        n_sig = length(sig_indices);
        n_nonsig = length(nonsig_indices);
        
        % Significant clusters (colored) - Row 1
        subplot(2, 2, g, 'Parent', fig_pooled);
        if n_sig > 0
            imagesc(bin_centers, 1:n_sig, sig_rate_mat);
            colormap(gca, 'jet');
        else
            axis off;
            text(0.5, 0.5, 'No significant clusters', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        xlabel('Position (cm)');
        if g == 1
            ylabel(sprintf('Significant Clusters (n=%d)', n_sig));
        else
            set(gca, 'YTickLabel', []);
        end
        set(gca, 'YTick', []);
        
        % Set x-axis limits
        if strcmp(group_names{g}, 'short')
            xlim([60, 120]);
        else
            xlim([0, 120]);
        end
        title(group_labels{g});
        
        % Non-significant clusters (grayscale) - Row 2
        subplot(2, 2, g+2, 'Parent', fig_pooled);
        if n_nonsig > 0
            imagesc(bin_centers, 1:n_nonsig, nonsig_rate_mat);
            colormap(gca, gray(256));
        else
            axis off;
            text(0.5, 0.5, 'No non-significant clusters', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        xlabel('Position (cm)');
        if g == 1
            ylabel(sprintf('Non-Significant Clusters (n=%d)', n_nonsig));
        else
            set(gca, 'YTickLabel', []);
        end
        set(gca, 'YTick', []);
        
        % Set x-axis limits
        if strcmp(group_names{g}, 'short')
            xlim([60, 120]);
        else
            xlim([0, 120]);
        end
    end
    
    % Adjust subplot positions
    all_axes = findall(fig_pooled, 'Type', 'axes');
    for ax_idx = 1:length(all_axes)
        ax = all_axes(ax_idx);
        pos = get(ax, 'Position');
        
        if pos(1) > 0.5  % Right column
            pos(3) = pos(3) * 0.5;
            pos(1) = 0.55;
        else  % Left column
            pos(3) = pos(3) * 1.05;
        end
        
        pos(4) = pos(4) * 0.82;
        pos(2) = max(0.08, pos(2) - 0.03);
        set(ax, 'Position', pos);
    end
    
    % Add figure title
    FigureTitle(fig_pooled, sprintf('Pooled Heatmaps - %s (%d clusters)', ...
        condition_label, n_clusters_cond));
    
    if save_figs
        ctl.figs.save_fig_to_join();
    else
        close(fig_pooled);
    end
end

% Create filtered heatmaps based on KW categorization
fprintf('\n=== Creating filtered heatmaps by KW categorization ===\n');

filter_combinations = {
    {'same_relative', 'diff_absolute', 'Same Relative / Different Absolute\n(relative distance)'};
    {'diff_relative', 'same_absolute', 'Different Relative / Same Absolute\n(absolute distance)'}
};

% Create one combined figure with all filtered heatmaps (4 rows × 2 columns)
% Rows: Region A filter 1, Region A filter 2, Region B filter 1, Region B filter 2
% Columns: Long trials, Short trials
fig_combined = ctl.figs.a4figure('landscape');

row_idx = 0;  % Track which row we're on

for cond_idx = 1:length(conditions)
    condition = conditions{cond_idx};
    condition_label = condition_labels{cond_idx};
    
    if isempty(pooled_data.(condition).cluster_ids)
        continue;
    end
    
    n_clusters_cond = length(pooled_data.(condition).cluster_ids);
    
    for filt_idx = 1:length(filter_combinations)
        rel_filter = filter_combinations{filt_idx}{1};  % 'same_relative' or 'diff_relative'
        abs_filter = filter_combinations{filt_idx}{2};  % 'same_absolute' or 'diff_absolute'
        filter_label = filter_combinations{filt_idx}{3};
        
        % Determine which clusters match this filter
        filter_mask = false(n_clusters_cond, 1);
        for c = 1:n_clusters_cond
            if ~isempty(pooled_data.(condition).dist_comparison{c})
                res = pooled_data.(condition).dist_comparison{c};
                
                % Check relative position criterion
                rel_match = false;
                if isfield(res, 'relative') && ~isnan(res.relative.same_distribution)
                    if strcmp(rel_filter, 'same_relative')
                        rel_match = res.relative.same_distribution;
                    else  % 'diff_relative'
                        rel_match = ~res.relative.same_distribution;
                    end
                end
                
                % Check absolute position criterion
                abs_match = false;
                if isfield(res, 'absolute') && ~isnan(res.absolute.same_distribution)
                    if strcmp(abs_filter, 'same_absolute')
                        abs_match = res.absolute.same_distribution;
                    else  % 'diff_absolute'
                        abs_match = ~res.absolute.same_distribution;
                    end
                end
                
                filter_mask(c) = rel_match && abs_match;
            end
        end
        
        row_idx = row_idx + 1;  % Increment row counter
        
        n_filtered = sum(filter_mask);
        fprintf('  %s - %s: %d clusters\n', condition_label, strrep(filter_label, sprintf('\n'), ' '), n_filtered);
        
        if n_filtered == 0
            fprintf('    No clusters match this filter, skipping\n');
            % Create empty subplots for this row
            for g = 1:2
                subplot(4, 2, (row_idx-1)*2 + g, 'Parent', fig_combined);
                axis off;
                text(0.5, 0.5, sprintf('No clusters\n%s - %s', condition_label, strrep(filter_label, sprintf('\n'), ' ')), ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', 9);
            end
            continue;
        end
        
        % Get filtered data (already sorted by peak position)
        rate_mat_long_filt = pooled_data.(condition).Q2_smooth_long(filter_mask, :);
        rate_mat_short_filt = pooled_data.(condition).Q2_smooth_short(filter_mask, :);
        
        % Sort by peak position
        peak_positions_filt = zeros(n_filtered, 1);
        for c = 1:n_filtered
            [~, max_idx] = max(rate_mat_long_filt(c, :));
            peak_positions_filt(c) = bin_centers_long(max_idx);
        end
        [~, sort_idx_filt] = sort(peak_positions_filt);
        rate_mat_long_filt = rate_mat_long_filt(sort_idx_filt, :);
        rate_mat_short_filt = rate_mat_short_filt(sort_idx_filt, :);
        
        % Plot long and short in columns, filters in rows (2×2 layout)
        for g = 1:2
            if g == 1
                rate_mat_filt = rate_mat_long_filt;
                bin_centers = bin_centers_long;
                group_label = 'Long (0-120 cm)';
            else
                rate_mat_filt = rate_mat_short_filt;
                bin_centers = bin_centers_short;
                group_label = 'Short (60-120 cm)';
            end
            
            % Min-max normalization
            min_vals = min(rate_mat_filt, [], 2);
            max_vals = max(rate_mat_filt, [], 2);
            norm_rate_mat_filt = (rate_mat_filt - min_vals) ./ (max_vals - min_vals);
            norm_rate_mat_filt(isnan(norm_rate_mat_filt) | isinf(norm_rate_mat_filt)) = 0;
            
            % Subplot index: 4 rows (2 regions × 2 filters) and 2 columns (long/short)
            subplot(4, 2, (row_idx-1)*2 + g, 'Parent', fig_combined);
            imagesc(bin_centers, 1:n_filtered, norm_rate_mat_filt);
            colormap(gca, 'jet');
            xlabel('Position (cm)', 'FontSize', 8);
            if g == 1
                % Split filter_label into separate lines for better readability
                if contains(filter_label, 'Same Relative')
                    ylabel_text = sprintf('%s\nSame Relative /\nDifferent Absolute\n(n=%d)', condition_label, n_filtered);
                else
                    ylabel_text = sprintf('%s\nDifferent Relative /\nSame Absolute\n(n=%d)', condition_label, n_filtered);
                end
                ylabel(ylabel_text, 'Interpreter', 'none', 'FontSize', 7);
            else
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'YTick', []);
            set(gca, 'FontSize', 8);
            
            % Set x-axis limits
            if g == 2
                xlim([60, 120]);
            else
                xlim([0, 120]);
            end
            
            % Add column title only for first row
            if row_idx == 1
                title(group_label, 'FontSize', 10);
            end
        end
    end
end

% Adjust subplot positions to make short trial columns narrower
all_axes = findall(fig_combined, 'Type', 'axes');
for ax_idx = 1:length(all_axes)
    ax = all_axes(ax_idx);
    pos = get(ax, 'Position');
    
    % Determine if this is a short (right column) or long (left column) subplot
    if pos(1) > 0.5  % Right column (short trials)
        % Make width 50% of original
        pos(3) = pos(3) * 0.5;
        pos(1) = 0.6;  % Move closer to left column
    else  % Left column (long trials)
        % Slightly increase width
        pos(3) = pos(3) * 1.1;
    end
    
    % Reduce height slightly and adjust vertical spacing
    pos(4) = pos(4) * 0.85;
    
    set(ax, 'Position', pos);
end

% Add overall figure title
FigureTitle(fig_combined, 'Pooled Heatmaps - KW Distribution Comparisons (All Regions)');

if save_figs
    ctl.figs.save_fig_to_join();
else
    close(fig_combined);
end

% Create pooled TTG heatmaps if data is available
fprintf('\n=== Creating pooled TTG heatmaps ===\n');
has_ttg_data = false;
for cond_idx = 1:length(conditions)
    condition = conditions{cond_idx};
    if isfield(pooled_data.(condition), 'ttg_rates_long') && ~isempty(pooled_data.(condition).ttg_rates_long)
        has_ttg_data = true;
        break;
    end
end

if has_ttg_data
    fig_ttg_pooled = ctl.figs.a4figure('landscape');
    
    for cond_idx = 1:length(conditions)
        condition = conditions{cond_idx};
        condition_label = condition_labels{cond_idx};
        
        if ~isfield(pooled_data.(condition), 'ttg_rates_long') || isempty(pooled_data.(condition).ttg_rates_long)
            fprintf('  No TTG data for condition %s, skipping\n', condition);
            continue;
        end
        
        fprintf('  Creating TTG heatmap for %s...\n', condition_label);
        
        % Get TTG data for this condition
        ttg_rate_mat_long = pooled_data.(condition).ttg_rates_long;
        ttg_rate_mat_short = pooled_data.(condition).ttg_rates_short;
        ttg_bin_centers = pooled_data.(condition).ttg_norm_bin_centers;
        n_clusters_ttg = size(ttg_rate_mat_long, 1);
        
        % Sort clusters by peak position in long trials (TTG)
        peak_positions_ttg = zeros(n_clusters_ttg, 1);
        for c = 1:n_clusters_ttg
            [~, max_idx] = max(ttg_rate_mat_long(c, :));
            peak_positions_ttg(c) = ttg_bin_centers(max_idx);
        end
        [~, sort_idx_ttg] = sort(peak_positions_ttg);
        
        % Apply sorting
        ttg_rate_mat_long = ttg_rate_mat_long(sort_idx_ttg, :);
        ttg_rate_mat_short = ttg_rate_mat_short(sort_idx_ttg, :);
        
        % Plot long and short TTG in columns
        for g = 1:2
            if g == 1
                ttg_rate_mat = ttg_rate_mat_long;
                group_label = 'Long Trials';
            else
                ttg_rate_mat = ttg_rate_mat_short;
                group_label = 'Short Trials';
            end
            
            % Min-max normalization per cluster
            min_vals = min(ttg_rate_mat, [], 2);
            max_vals = max(ttg_rate_mat, [], 2);
            norm_ttg_rate_mat = (ttg_rate_mat - min_vals) ./ (max_vals - min_vals);
            norm_ttg_rate_mat(isnan(norm_ttg_rate_mat) | isinf(norm_ttg_rate_mat)) = 0;
            
            % Create subplot (2 rows for 2 conditions, 2 columns for long/short)
            subplot(2, 2, (cond_idx-1)*2 + g, 'Parent', fig_ttg_pooled);
            imagesc(ttg_bin_centers, 1:n_clusters_ttg, norm_ttg_rate_mat);
            colormap(gca, 'jet');
            xlabel('Normalized Time-to-Goal (%)', 'FontSize', 9);
            if g == 1
                ylabel(sprintf('%s\n(n=%d)', condition_label, n_clusters_ttg), 'FontSize', 9);
            else
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'YTick', []);
            set(gca, 'XDir', 'reverse');  % 100% on left, 0% on right
            set(gca, 'FontSize', 8);
            xlim([min(ttg_bin_centers), max(ttg_bin_centers)]);
            
            % Add column title only for first row
            if cond_idx == 1
                title(group_label, 'FontSize', 10);
            end
        end
        
        fprintf('    Created TTG heatmap for %s with %d clusters\n', condition_label, n_clusters_ttg);
    end
    
    % Adjust subplot positions to make them larger
    all_axes = findall(fig_ttg_pooled, 'Type', 'axes');
    for ax_idx = 1:length(all_axes)
        ax = all_axes(ax_idx);
        pos = get(ax, 'Position');
        
        % Increase width and height
        pos(3) = pos(3) * 1.15;  % Width
        pos(4) = pos(4) * 1.1;   % Height
        
        % Adjust positioning
        if pos(1) > 0.5  % Right column
            pos(1) = 0.52;
        else  % Left column
            pos(1) = pos(1) - 0.02;
        end
        pos(2) = pos(2) - 0.03;
        
        set(ax, 'Position', pos);
    end
    
    % Add overall figure title
    FigureTitle(fig_ttg_pooled, 'Pooled TTG Heatmaps (All Regions)');
    
    if save_figs
        ctl.figs.save_fig_to_join();
    else
        close(fig_ttg_pooled);
    end
    
    % Create filtered TTG heatmaps (peaks between 10% and 90% in long condition)
    fprintf('\n=== Creating filtered TTG heatmaps (peaks 10-90%%) ===\n');
    fig_ttg_filtered = ctl.figs.a4figure('landscape');
    
    for cond_idx = 1:length(conditions)
        condition = conditions{cond_idx};
        condition_label = condition_labels{cond_idx};
        
        if ~isfield(pooled_data.(condition), 'ttg_rates_long') || isempty(pooled_data.(condition).ttg_rates_long)
            fprintf('  No TTG data for condition %s, skipping\n', condition);
            continue;
        end
        
        fprintf('  Creating filtered TTG heatmap for %s...\n', condition_label);
        
        % Get TTG data for this condition
        ttg_rate_mat_long = pooled_data.(condition).ttg_rates_long;
        ttg_rate_mat_short = pooled_data.(condition).ttg_rates_short;
        ttg_bin_centers = pooled_data.(condition).ttg_norm_bin_centers;
        n_clusters_ttg = size(ttg_rate_mat_long, 1);
        
        % Find peak positions and filter for 10% <= peak <= 90%
        peak_positions_ttg = zeros(n_clusters_ttg, 1);
        for c = 1:n_clusters_ttg
            [~, max_idx] = max(ttg_rate_mat_long(c, :));
            peak_positions_ttg(c) = ttg_bin_centers(max_idx);
        end
        
        % Filter for peaks between 10% and 90%
        filter_mask = (peak_positions_ttg >= 10) & (peak_positions_ttg <= 90);
        n_filtered = sum(filter_mask);
        
        fprintf('    %s: %d/%d clusters with peaks between 10-90%%\n', condition_label, n_filtered, n_clusters_ttg);
        
        if n_filtered == 0
            fprintf('    No clusters match filter, skipping\n');
            % Create empty subplots for this row
            for g = 1:2
                subplot(2, 2, (cond_idx-1)*2 + g, 'Parent', fig_ttg_filtered);
                axis off;
                text(0.5, 0.5, sprintf('No clusters\n%s\n(10-90%% peaks)', condition_label), ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', 9);
            end
            continue;
        end
        
        % Apply filter
        ttg_rate_mat_long_filt = ttg_rate_mat_long(filter_mask, :);
        ttg_rate_mat_short_filt = ttg_rate_mat_short(filter_mask, :);
        peak_positions_filt = peak_positions_ttg(filter_mask);
        
        % Sort by peak position
        [~, sort_idx_filt] = sort(peak_positions_filt);
        ttg_rate_mat_long_filt = ttg_rate_mat_long_filt(sort_idx_filt, :);
        ttg_rate_mat_short_filt = ttg_rate_mat_short_filt(sort_idx_filt, :);
        
        % Plot long and short TTG in columns
        for g = 1:2
            if g == 1
                ttg_rate_mat_filt = ttg_rate_mat_long_filt;
                group_label = 'Long Trials';
            else
                ttg_rate_mat_filt = ttg_rate_mat_short_filt;
                group_label = 'Short Trials';
            end
            
            % Min-max normalization per cluster
            min_vals = min(ttg_rate_mat_filt, [], 2);
            max_vals = max(ttg_rate_mat_filt, [], 2);
            norm_ttg_rate_mat_filt = (ttg_rate_mat_filt - min_vals) ./ (max_vals - min_vals);
            norm_ttg_rate_mat_filt(isnan(norm_ttg_rate_mat_filt) | isinf(norm_ttg_rate_mat_filt)) = 0;
            
            % Create subplot (2 rows for 2 conditions, 2 columns for long/short)
            subplot(2, 2, (cond_idx-1)*2 + g, 'Parent', fig_ttg_filtered);
            imagesc(ttg_bin_centers, 1:n_filtered, norm_ttg_rate_mat_filt);
            colormap(gca, 'jet');
            xlabel('Normalized Time-to-Goal (%)', 'FontSize', 9);
            if g == 1
                ylabel(sprintf('%s\nPeaks 10-90%%\n(n=%d)', condition_label, n_filtered), 'FontSize', 8);
            else
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'YTick', []);
            set(gca, 'XDir', 'reverse');  % 100% on left, 0% on right
            set(gca, 'FontSize', 8);
            xlim([min(ttg_bin_centers), max(ttg_bin_centers)]);
            
            % Add column title only for first row
            if cond_idx == 1
                title(group_label, 'FontSize', 10);
            end
        end
        
        fprintf('    Created filtered TTG heatmap for %s with %d clusters\n', condition_label, n_filtered);
    end
    
    % Adjust subplot positions to make them larger
    all_axes = findall(fig_ttg_filtered, 'Type', 'axes');
    for ax_idx = 1:length(all_axes)
        ax = all_axes(ax_idx);
        pos = get(ax, 'Position');
        
        % Increase width and height
        pos(3) = pos(3) * 1.15;  % Width
        pos(4) = pos(4) * 1.1;   % Height
        
        % Adjust positioning
        if pos(1) > 0.5  % Right column
            pos(1) = 0.52;
        else  % Left column
            pos(1) = pos(1) - 0.02;
        end
        pos(2) = pos(2) - 0.03;
        
        set(ax, 'Position', pos);
    end
    
    % Add overall figure title
    FigureTitle(fig_ttg_filtered, 'Pooled TTG Heatmaps - Filtered (Peaks 10-90%)');
    
    if save_figs
        ctl.figs.save_fig_to_join();
    else
        close(fig_ttg_filtered);
    end
else
    fprintf('  No TTG data found. TTG cache files may not exist.\n');
    fprintf('  Re-run spatial_firing_rate_profile.m to generate TTG caches.\n');
end

% Join figures if requested
if save_figs
    fprintf('\nJoining figures into single PDF...\n');
    fname = 'pooled_heatmaps.pdf';
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end

fprintf('\n=== Pooled heatmap analysis complete ===\n');

% Restore default figure visibility
set(groot, 'DefaultFigureVisible', 'on');
