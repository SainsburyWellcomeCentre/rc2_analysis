function mask = compute_camera_motion_mask(global_motion, plot_checks)

VariableDefault('plot_checks', false);

min_dur = 0.2;

[counts, x]     = histcounts(global_motion, 'binwidth', 1);
counts          = smooth(counts, 100);
[~, imax]       = max(counts);

% repeat values below x(imax)
rep_noise        = [global_motion(global_motion < x(imax)), ...
    x(imax) + (x(imax) - global_motion(global_motion < x(imax)))];

%% fit gaussian to "noise"
[counts2, x2]       = histcounts(rep_noise, 'binwidth', 1);
x2              = (x2(1:end-1) + x2(2:end))/2;
f               = fit(x2(:), counts2(:), 'gauss1');

% determine threshold
threshold       = f.b1 + norminv(0.995)*f.c1;

%% plot fit and threshold
if plot_checks
    m               = global_motion;
    m_restricted    = m(m > (-0.5*1e4 + x(imax)) & m < (0.5*1e4 + x(imax)));
    [n_, x_]        = histcounts(m_restricted, 'binwidth', 10);
    
    % for fit
    t               = (-0.5*1e4 + x(imax)):(0.5*1e4 + x(imax));
    y               = max(n_)*exp(-((t-f.b1)/f.c1).^2);
    
    figure;
    histogram(m_restricted, 'binwidth', 10);
    hold on;
    plot(t, y, 'r')
    line(threshold * [1, 1], get(gca, 'ylim'), 'color', 'r');
end

%% transform
m                   = global_motion;
mask                = m < threshold;

cc = bwconncomp(mask);
s = cellfun(@(x)(x(1)), cc.PixelIdxList);
e = cellfun(@(x)(x(end)), cc.PixelIdxList);

dur = (e - s)/10e3;

rm_idx = dur < min_dur;

mask(cat(1, cc.PixelIdxList{rm_idx})) = false;
mask = ~mask;