function plot_lfp_power(power_depths, power_profiles, cluster_depths, boundaries, region_str, anatomy_l5, ephys_l5)

x_limits = [0, 1.2];

% one figure
figure('position', [680, 84, 1100, 894]);

% left most plot axis
subplot(1, 3, 1); hold on;

% normalize by dividing by max LFP power seen
norm_power_profiles = power_profiles/max(power_profiles(:));

% draw the LFP power line
plot(bsxfun(@rdivide, power_profiles, max(power_profiles, [], 1)), power_depths, 'color', 'k', 'linewidth', 2); 

% some formatting
set(gca, 'ylim', [0, 2500], 'xlim', [0, 1.2])   

xlabel('Norm. power (500-1250Hz)');  
ylabel('Distance from probe tip (\mum)')

set(gca, 'plotboxaspectratio', [1, 3, 1])
box off;
title('LFP power');

% plot line for middle of L5 from ephys
line(get(gca, 'xlim'), ephys_l5*[1, 1], 'color', 'r', 'linewidth', 2);
% text for distance of line from tip
text(min(get(gca, 'xlim')), ephys_l5, sprintf('peak LFP: %i', mid_l5_ephys), ...
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
print(fullfile(ks_dir, 'tracks', 'lfp_power.pdf'), '-bestfit', '-dpdf');

% write parameters
if ~isfolder(fullfile(ks_dir, 'tracks'))
    mkdir(fullfile(ks_dir, 'tracks'))
end

% write offset
fid = fopen(fullfile(ks_dir, 'tracks', 'offset.txt'), 'w');
fprintf(fid, '%.0f', delta_l5);
fclose(fid);

% write parameters
fid = fopen(fullfile(ks_dir, 'tracks', 'lfp_params.txt'), 'w');
fprintf(fid, 'offset_s = %.2f\r\n', params.offset_s);
fprintf(fid, 'offset_s = %.2f\r\n', params.duration_s);
fprintf(fid, 'bad_channels = ');
for i = 1 : length(params.bad_channels)
    fprintf(fid, '%i,', params.bad_channels(i));
end
fprintf('\b\r\n');
fprintf('entire_recording = %i\r\n', strcmp(batches_to_use, 'all'));
fprintf('bin_fname = %s\r\n', lf_fname);
fprintf('search_above = %s\r\n', search_above);
fprintf(fid, 'batches_to_use = ');
for i = 1 : length(batches_to_use)
    fprintf(fid, '%i,', batches_to_use(i));
end
fprintf('\b\r\n');




