function [x_med, y_med, p, direction] = compare_groups_with_signrank(x, y)

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