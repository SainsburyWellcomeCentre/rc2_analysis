function plot_beta_coefficients_subplot(beta_matrix, sig_vector, is_linear_vector, subplot_title)
% PLOT_BETA_COEFFICIENTS_SUBPLOT Plot beta coefficients for multiple clusters in a subplot
%
% Inputs:
%   beta_matrix       - n_clusters x 3 matrix of beta coefficients [x², x¹, x⁰]
%   sig_vector        - n_clusters logical vector indicating significance
%   is_linear_vector  - n_clusters logical vector indicating linear fit
%   subplot_title     - Title for the subplot
%
% Creates a plot with:
%   - X-axis: beta0, beta1, beta2 (reversed order from polyfit output)
%   - Y-axis: beta values (log scale)
%   - Lines connecting beta values for each cluster
%   - Color coding: green (linear sig), purple (quadratic sig), gray (not sig)

n_clusters = size(beta_matrix, 1);
x_positions = [1, 2, 3];  % beta0, beta1, beta2

hold on;
for c = 1:n_clusters
    % Reverse order: polyfit gives [x², x¹, x⁰], we want [x⁰, x¹, x²]
    beta_values = beta_matrix(c, [3, 2, 1]);
    
    % Skip if all values are NaN
    if all(isnan(beta_values))
        continue;
    end
    
    % Set color based on significance and linearity
    if sig_vector(c)
        if is_linear_vector(c)
            line_color = [0, 0.6, 0.3];  % Green for significant linear
            marker_face = [0, 0.6, 0.3];
        else
            line_color = [0.5, 0, 0.5];  % Purple for significant quadratic
            marker_face = [0.5, 0, 0.5];
        end
    else
        line_color = [0.7, 0.7, 0.7];  % Gray for not significant
        marker_face = [0.7, 0.7, 0.7];
    end
    
    % Plot line connecting the points
    plot(x_positions, beta_values, '-o', ...
        'Color', line_color, ...
        'MarkerFaceColor', marker_face, ...
        'MarkerEdgeColor', line_color, ...
        'MarkerSize', 4, ...
        'LineWidth', 1);
end
hold off;

xlim([0.5, 3.5]);
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'\beta_0', '\beta_1', '\beta_2'});
set(gca, 'YScale', 'log');
ylabel('Beta value (log scale)');
title(subplot_title, 'Interpreter', 'none');
grid on;

end
