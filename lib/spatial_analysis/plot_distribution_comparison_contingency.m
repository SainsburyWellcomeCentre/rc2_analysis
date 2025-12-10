function fig = plot_distribution_comparison_contingency(dist_comparison_results, cluster_ids, probe_id)
% PLOT_DISTRIBUTION_COMPARISON_CONTINGENCY Create contingency table heatmap
%
%   fig = plot_distribution_comparison_contingency(dist_comparison_results, cluster_ids, probe_id)
%
%   Creates a 2×2 contingency table showing the relationship between
%   absolute and relative position distribution comparisons:
%       - Same/Same: Same distribution in both absolute and relative position
%       - Same/Different: Same in absolute, different in relative
%       - Different/Same: Different in absolute, same in relative
%       - Different/Different: Different in both absolute and relative
%
%   Inputs:
%       dist_comparison_results - Cell array of comparison results for each cluster
%       cluster_ids             - Array of cluster IDs
%       probe_id                - String identifier for the probe (for title)
%
%   Outputs:
%       fig - Figure handle

    n_clusters_total = length(cluster_ids);
    
    % Initialize contingency table
    % Rows: Absolute (1=Same, 2=Different)
    % Cols: Relative (1=Same, 2=Different)
    contingency = zeros(2, 2);
    
    % Classify each cluster (only spatially tuned clusters have non-empty results)
    abs_labels = cell(n_clusters_total, 1);
    rel_labels = cell(n_clusters_total, 1);
    n_spatially_tuned = 0;
    
    for c = 1:n_clusters_total
        comp = dist_comparison_results{c};
        
        if isempty(comp)
            % Skip non-spatially tuned clusters
            abs_labels{c} = 'NA';
            rel_labels{c} = 'NA';
            continue;
        end
        
        % Count this as a spatially tuned cluster
        n_spatially_tuned = n_spatially_tuned + 1;
        
        % Absolute position comparison
        if isfield(comp, 'absolute') && isstruct(comp.absolute) && ...
           isfield(comp.absolute, 'same_distribution') && ~isnan(comp.absolute.same_distribution)
            if comp.absolute.same_distribution
                abs_labels{c} = 'Same';
                abs_idx = 1;
            else
                abs_labels{c} = 'Different';
                abs_idx = 2;
            end
        else
            abs_labels{c} = 'NA';
            abs_idx = 0;
        end
        
        % Relative position comparison
        if isfield(comp, 'relative') && isstruct(comp.relative) && ...
           isfield(comp.relative, 'same_distribution') && ~isnan(comp.relative.same_distribution)
            if comp.relative.same_distribution
                rel_labels{c} = 'Same';
                rel_idx = 1;
            else
                rel_labels{c} = 'Different';
                rel_idx = 2;
            end
        else
            rel_labels{c} = 'NA';
            rel_idx = 0;
        end
        
        % Update contingency table
        if abs_idx > 0 && rel_idx > 0
            contingency(abs_idx, rel_idx) = contingency(abs_idx, rel_idx) + 1;
        end
    end
    
    % Create figure
    fig = figure('Position', [100, 100, 800, 700], 'Visible', 'off');
    
    % Plot contingency table as heatmap
    imagesc(contingency);
    colormap('gray');
    colorbar;
    
    % Set axis labels (larger font sizes)
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Same', 'Different'});
    set(gca, 'YTick', 1:2, 'YTickLabel', {'Same', 'Different'});
    xlabel('Relative Position', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Absolute Position', 'FontSize', 20, 'FontWeight', 'bold');
    
    % Add count text to each cell (larger font size)
    for i = 1:2
        for j = 1:2
            text(j, i, sprintf('%d', contingency(i, j)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 48, 'FontWeight', 'bold', 'Color', [0.2 0.5 1]);
        end
    end
    
    % Add title with probe ID (showing only spatially tuned clusters)
    title(sprintf('Distribution Comparison: %s (n=%d spatially tuned clusters)', probe_id, n_spatially_tuned), ...
          'FontSize', 16, 'FontWeight', 'bold');
    
    % Adjust axes for better display (larger tick labels)
    axis square;
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');
end
