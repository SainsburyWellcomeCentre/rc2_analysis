function l5_from_tip = search_for_L5(lfp_power, channel_from_tip, search_above, search_below)

new_x = channel_from_tip(1):channel_from_tip(end);
lfp_power_interp = interp1(channel_from_tip, lfp_power, new_x, 'spline');

search_idx = find(new_x > search_above & new_x < search_below);

[~, max_idx] = max(lfp_power_interp(search_idx));

l5_from_tip = new_x(search_idx(1) + max_idx - 1);
