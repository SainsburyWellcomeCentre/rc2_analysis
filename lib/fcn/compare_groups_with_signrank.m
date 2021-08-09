function [x_med, y_med, p, is_increase] = compare_groups_with_signrank(x, y)

x_med = nanmedian(x);
y_med = nanmedian(y);
[p, ~, stats] = signrank(x, y);

% if the medians are equal, but there is a signficant
% difference between the two conditions
if x_med == y_med && p < 0.05
    
    % perform the test again but swap the conditions
    [~, ~, stats_opp] = signrank(y, x);
    
    if stats_opp.signedrank > stats.signedrank
        is_increase = true;
    elseif stats_opp.signedrank < stats.signedrank
        is_increase = false;
    else
        error('Signed ranks are equal?');
    end
else
    is_increase = y_med > x_med;
end