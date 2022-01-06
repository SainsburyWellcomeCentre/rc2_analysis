function [x_med, y_med, p, direction] = compare_groups_with_signrank(x, y)
% COMPARE_GROUPS_WITH_SIGNRANK Compares two groups of values with a
% signrank (+ extra information if medians of data are equal)
%
%   [X_MEDIAN, Y_MEDIAN, P, DIRECTION] = compare_groups_with_signrank(x, y)
%   takes two sets of paired values X and Y (which must be equal length)
%   and performs a Wilcoxon sign rank test on them.
%
%   The median of the X and Y data are returned in X_MEDIAN and Y_MEDIAN
%   (which ignores NaN values), and the p-value of the test is returned in
%   P.
%
%   The reason this function exists is because occassionally we see two
%   groups of values with median = 0, but which differ significantly. 
%   In this case, we also return a DIRECTION which indicates
%   whether the non-zero values of the Y group are greater than those of the X group.
%
%   To get the DIRECTION we perform the sign rank test twice, signrank(X,
%   Y) and signrank(Y, X). If the test statistic, T, is greater for
%   signrank(y, x) then the data in Y are considered larger in magnitude
%   than the data in X, and if T is greater for signrank(x, y) then the
%   data in X are considered larger in magnitude than in Y.

x_med = nanmedian(x);
y_med = nanmedian(y);
[p, ~, stats] = signrank(x, y);
[~, ~, stats_opp] = signrank(y, x);

% if the medians are equal, but there is a signficant
% difference between the two conditions
if p < 0.05
    
    % perform the test again but swap the conditions
    if x_med == y_med
        if stats_opp.signedrank > stats.signedrank
            direction = 1;
        elseif stats_opp.signedrank < stats.signedrank
            direction = -1;
        else
            error('Signed ranks are equal?');
        end
    elseif y_med > x_med
        direction = 1;
    elseif y_med < x_med
        direction = -1;
    end
else
    direction = 0;
end