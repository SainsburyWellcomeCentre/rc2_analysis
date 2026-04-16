% plot_model_selection_examples.m
% Generates example curves for each model used in ModelSelectionTuning.

x = linspace(-5, 5, 200)';

fig = figure('Position', [100, 100, 1400, 550]);

% -------------------------------------------------------------------------
% 1. Linear: b0 + b1*x
% -------------------------------------------------------------------------
subplot(2, 3, 1);
beta_linear = [0.6, 6];   % [b1, b0] (polyval order) — min at x=-5: y=3
y_linear = polyval(beta_linear, x);
plot(x, y_linear, 'LineWidth', 2, 'Color', 'k');
xlabel('Stimulus value'); ylabel('Firing rate');
title('Linear');
subtitle('y = b_0 + b_1 x');
box off;

% -------------------------------------------------------------------------
% 2. Quadratic: b0 + b1*x + b2*x^2
% -------------------------------------------------------------------------
subplot(2, 3, 2);
beta_quad = [-0.3, 0, 10];   % [b2, b1, b0] concave-down, peak at x=0 — min at edges: y=2.5
y_quad = polyval(beta_quad, x);
plot(x, y_quad, 'LineWidth', 2, 'Color', 'k');
xlabel('Stimulus value'); ylabel('Firing rate');
title('Quadratic');
subtitle('y = b_0 + b_1 x + b_2 x^2');
box off;

% -------------------------------------------------------------------------
% 3. Cubic: b0 + b1*x + b2*x^2 + b3*x^3
% -------------------------------------------------------------------------
subplot(2, 3, 3);
beta_cubic = [0.02, 0, 0.3, 6];    % [b3, b2, b1, b0] S-shape — min at x=-5: y=2
y_cubic = polyval(beta_cubic, x);
plot(x, y_cubic, 'LineWidth', 2, 'Color', 'k');
xlabel('Stimulus value'); ylabel('Firing rate');
title('Cubic');
subtitle('y = b_0 + b_1 x + b_2 x^2 + b_3 x^3');
box off;

% -------------------------------------------------------------------------
% 4. Gaussian: A*exp(-(x-mu)^2 / (2*sigma^2)) + baseline
% -------------------------------------------------------------------------
subplot(2, 3, 4);
A = 5; mu = 0.5; sigma = 1.2; baseline = 1;
y_gauss = A * exp(-(x - mu).^2 / (2 * sigma^2)) + baseline;
plot(x, y_gauss, 'LineWidth', 2, 'Color', 'k');
xlabel('Stimulus value'); ylabel('Firing rate');
title('Gaussian');
subtitle('y = A exp(-(x-\mu)^2 / 2\sigma^2) + b');
box off;

% -------------------------------------------------------------------------
% 5. Asymmetric Gaussian: R_max*exp(-(x-x_max)^2 / sigma_x)
%    sigma_x = sigma_minus when x < x_max, sigma_plus when x >= x_max
% -------------------------------------------------------------------------
subplot(2, 3, 5);
R_max = 6; x_max = 0.5; sigma_minus = 0.8; sigma_plus = 3.5;
y_asym = zeros(size(x));
left  = x < x_max;
right = x >= x_max;
y_asym(left)  = R_max * exp(-(x(left)  - x_max).^2 / sigma_minus);
y_asym(right) = R_max * exp(-(x(right) - x_max).^2 / sigma_plus);
plot(x, y_asym, 'LineWidth', 2, 'Color', 'k');
xlabel('Stimulus value'); ylabel('Firing rate');
title('Asymmetric Gaussian');
subtitle('y = R_{max} exp(-(x-x_{max})^2 / \sigma_x),  \sigma_- \neq \sigma_+');
box off;

% -------------------------------------------------------------------------
% 6. Sigmoid: A / (1 + exp(-k*(x - x0))) + baseline
% -------------------------------------------------------------------------
subplot(2, 3, 6);
A_sig = 5; k = 1.5; x0 = 0; baseline_sig = 0.5;
y_sig = A_sig ./ (1 + exp(-k * (x - x0))) + baseline_sig;
plot(x, y_sig, 'LineWidth', 2, 'Color', 'k');
xlabel('Stimulus value'); ylabel('Firing rate');
title('Sigmoid');
subtitle('y = A / (1 + exp(-k(x - x_0))) + b');
box off;

sgtitle('ModelSelectionTuning — candidate models', 'FontSize', 14, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% Save as PDF
% -------------------------------------------------------------------------
output_path = fullfile(fileparts(mfilename('fullpath')), 'figures', 'model_selection_examples.pdf');
if ~exist(fileparts(output_path), 'dir')
    mkdir(fileparts(output_path));
end
exportgraphics(fig, output_path, 'ContentType', 'vector');
fprintf('Figure saved to: %s\n', output_path);
