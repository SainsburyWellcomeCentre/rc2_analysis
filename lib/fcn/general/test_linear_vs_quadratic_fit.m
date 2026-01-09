function is_linear = test_linear_vs_quadratic_fit(x_data, y_data, beta_quad, alpha)
% TEST_LINEAR_VS_QUADRATIC_FIT Perform F-test to determine if quadratic fit is significantly better than linear
%
% Inputs:
%   x_data    - Independent variable data (column vector)
%   y_data    - Dependent variable data (column vector)
%   beta_quad - Coefficients from quadratic fit [x², x¹, x⁰] (from polyfit)
%   alpha     - Significance level (default: 0.05)
%
% Outputs:
%   is_linear - True if quadratic doesn't significantly improve fit (p > alpha)
%               False if quadratic is significantly better
%
% Uses F-test for nested models to compare linear vs quadratic polynomial fits

if nargin < 4
    alpha = 0.05;
end

% Default: assume not linear (i.e., quadratic is better)
is_linear = false;

% Remove NaN values
valid_idx = ~isnan(x_data) & ~isnan(y_data);
x_data = x_data(valid_idx);
y_data = y_data(valid_idx);
n = length(y_data);

% Need at least 4 points to fit quadratic and test
if n <= 3
    return;
end

% Linear fit
beta_lin = polyfit(x_data, y_data, 1);
yfit_lin = polyval(beta_lin, x_data);
ss_res_lin = sum((y_data - yfit_lin).^2);

% Quadratic fit (already computed)
yfit_quad = polyval(beta_quad, x_data);
ss_res_quad = sum((y_data - yfit_quad).^2);

% F-test for nested models
df1 = 1;  % difference in parameters (3 - 2)
df2 = n - 3;  % residual degrees of freedom for quadratic

if ss_res_quad > 0 && df2 > 0
    F = ((ss_res_lin - ss_res_quad) / df1) / (ss_res_quad / df2);
    p_quadratic = 1 - fcdf(F, df1, df2);
    
    % If p > alpha, quadratic doesn't significantly improve fit
    is_linear = (p_quadratic > alpha);
end

end
