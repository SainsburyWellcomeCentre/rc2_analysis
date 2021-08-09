clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% location of LFP recording
animal_id = 'CAA-1113932';
probe_suffix = 'rec1_rec2';
restrict_batches = true;
batches_to_use = 1:10;
search_above = 1200;
search_below = inf;
compute_clusters = true;
basedir = 'E:\mateoData_probe\janelia_pipeline'; % 'F:\steinmetz_probe';  % 'E:\mateoData_probe\janelia_pipeline';

% Build filename for LFP data...
lf_fname = fullfile(basedir, animal_id, ...
    sprintf('%s_%s_g0', animal_id, probe_suffix), ...
    sprintf('%s_%s_g0_imec0', animal_id, probe_suffix), ...
    sprintf('%s_%s_g0_t0.imec0.lf.bin', animal_id, probe_suffix));

% KS-dir
ks_dir = fullfile(basedir, animal_id, 'output', ...
    sprintf('catgt_%s_%s_g0', animal_id, probe_suffix), ...
    sprintf('%s_%s_g0_imec0', animal_id, probe_suffix), ...
    'imec0_ks2');

track_fname = 'track_0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the LFP on all channels
[lfp_power, p, channel_list] = get_lfp_power(lf_fname);


%% for some mice it is beneficial to restrict the sampling of the LFP
if restrict_batches
    lfp_power_avg = mean(bsxfun(@rdivide, p(:, batches_to_use), max(p(:, batches_to_use), [], 1)), 2);%mean(p(:, batches_to_use), 2);
end

% distance from tip of each channel
channel_from_tip = get_channel_distance_from_tip(channel_list-1);

% take every other channel
lfp_power_plot = lfp_power_avg(1:2:end);
channel_from_tip_plot = channel_from_tip(1:2:end);

% search for the L5 peak
mid_l5_ephys = search_for_L5(lfp_power_plot, channel_from_tip_plot, search_above, search_below);


% multiunit activity
if compute_clusters
    clusters = format_clusters(ks_dir);
end

% non-noise clusters
idx = [clusters(:).class] == "good";

% distance of non-noise clusters from tip
cluster_dist_from_tip = [clusters(idx).distance_from_probe_tip];

% load the probe track csv file
probe_track = read_track_csv(ks_dir, track_fname);

if ~isempty(probe_track)
    % use to find positions of the boundaries from probe tip
    [boundaries, id, str] = find_layer_boundaries(probe_track);
    
    % find the middle of L5 using the anatomy boundaries
    mid_l5_anatomy = mid_L5_from_tip(boundaries, str);
    
    % difference between L5 from ephys and anatomy
    delta_l5 = mid_l5_ephys - mid_l5_anatomy;
    
    % use delta to correct the position of the boundaries
    boundaries_corrected = boundaries + delta_l5;
    
else
    
    boundaries = [];
    boundaries_corrected = [];
    mid_l5_anatomy = nan;
    delta_l5 = nan;
    str = {};
    id = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% one figure
figure('position', [680, 84, 1100, 894]);

% left most plot axis
subplot(1, 3, 1); hold on;

% normalize by dividing by max LFP power seen
norm_lfp = lfp_power_plot/max(lfp_power_plot);

% draw the LFP power line
plot(bsxfun(@rdivide, p(1:2:end, batches_to_use), max(p(1:2:end, batches_to_use), [], 1)), channel_from_tip_plot, 'color', 'k', 'linewidth', 2); %p(1:2:end, :), 

% somee formatting
set(gca, 'ylim', [0, 2500], 'xlim', [0, 1.2])
ylabel('Distance from probe tip (\mum)')
xlabel('Norm. power (500-1250Hz)');
set(gca, 'plotboxaspectratio', [1, 3, 1])
box off;
title('LFP power');

% plot line for middle of L5 from ephys
line(get(gca, 'xlim'), mid_l5_ephys*[1, 1], 'color', 'r', 'linewidth', 2);
% text for distance of line from tip
text(min(get(gca, 'xlim')), mid_l5_ephys, sprintf('peak LFP: %i', mid_l5_ephys), ...
    'verticalalignment', 'bottom', 'horizontalalignment', 'left', 'color', 'r');

% plot line for middle of L5 from anatomy
line(get(gca, 'xlim'), mid_l5_anatomy*[1, 1], 'color', 'b', 'linewidth', 2);
% text for distance of line from tip
text(min(get(gca, 'xlim')), mid_l5_anatomy, sprintf('mid VISp5: %i', round(mid_l5_anatomy)), ...
    'verticalalignment', 'top', 'horizontalalignment', 'left', 'color', 'b');

% plot the boundaries
for i = 1 : length(boundaries)
    % boundary
    line(get(gca, 'xlim'), boundaries(i)*[1, 1], 'color', 'k', 'linestyle', '--');
    % text in middle of boundary
    if i < length(boundaries)
        y = sum(boundaries([i, i+1]))/2;
        text(max(get(gca, 'xlim')), y, str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    end
end



% middle plot axis
subplot(1, 3, 2); hold on;

% again plot LFP power line
plot(norm_lfp, channel_from_tip_plot, 'color', 'k', 'linewidth', 2); %p(1:2:end, :), 

% formatting
set(gca, 'ylim', [0, 2500], 'xlim', [0, 1.2])
ylabel('Distance from probe tip (\mum)')
xlabel('Norm. power (500-1250Hz)');
set(gca, 'plotboxaspectratio', [1, 3, 1])
box off
title('LFP power (boundaries corrected)');

% plot line for middle of L5 from ephys (also will be equal to mid L5 from
% corrected boundaries)
line(get(gca, 'xlim'), mid_l5_ephys*[1, 1], 'color', 'r', 'linewidth', 2);

% plot the corrected boundaries
for i = 1 : length(boundaries_corrected)
    line(get(gca, 'xlim'), boundaries_corrected(i)*[1, 1], 'color', 'k', 'linestyle', '--');
    if i < length(boundaries_corrected)
        y = sum(boundaries_corrected([i, i+1]))/2;
        text(max(get(gca, 'xlim')), y, str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    end
end

% right most plot axis
subplot(1, 3, 3); hold on;

% histogram of number of clusters found at different depths
histogram(cluster_dist_from_tip, 40, 'orientation', 'horizontal')

% some formatting
set(gca, 'plotboxaspectratio', [1, 3, 1])
set(gca, 'ylim', [0, 2500])
xlabel('# units')
box off
title('MUA (boundaries corrected)')

% plot line for middle of L5 from ephys
line(get(gca, 'xlim'), mid_l5_ephys*[1, 1], 'color', 'r', 'linewidth', 2);

% plot the corrected boundari
for i = 1 : length(boundaries_corrected)
    line(get(gca, 'xlim'), boundaries_corrected(i)*[1, 1], 'color', 'k', 'linestyle', '--');
    if i < length(boundaries_corrected)
        y = sum(boundaries_corrected([i, i+1]))/2;
        text(max(get(gca, 'xlim')), y, str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    end
end

% give the whole figure a title
ft = FigureTitle(gcf, sprintf('%s_%s LFP power', animal_id, probe_suffix));
ft.y_position = 0.98;

% for saving set the paper orientation to landscape
set(gcf, 'paperorientation', 'landscape');

% save the figure to the current directory...
print(sprintf('%s_%s_lfp_power.pdf', animal_id, probe_suffix), '-bestfit', '-dpdf');
