probe_ids = {'CAA-1110262_rec1_rec2_rec3';
             'CAA-1110263_restricted_rec1_rec2_rec3';
             'CAA-1110264_rec1_rec2';
             'CAA-1110265_restricted_rec1_rec2_rec3';
             'CAA-1112222_rec1_rec2_rec3';
             'CAA-1112223_rec1_rec2_rec3';
             'CAA-1112224_rec1_rec2_rec3';
             'CAA-1112416_rec1_rec2_rec3';
             'CAA-1112417_rec1_rec2_rec3';
             'CA_176_1_rec1_rec2_rec3';
             'CA_176_3_rec1_rec2_rec3';
             'CAA-1112872_rec1_rec1b_rec2_rec3';
             'CAA-1112874_rec1_rec2_rec3';
             'CAA-1113219_rec1_rec2_rec3';
             'CAA-1113220_rec1_rec2_rec3';
             'CAA-1113222_rec1_rec2_rec3';
             'CAA-1114977_rec1_rec2_rec3';
             'CAA-1114978_rec1_rec2';
             'CAA-1114979_rec1_rec2';
             'CAA-1114980_rec1_rec2'};

ctl         = RC2Preprocess();
ctl.setup_figures({'high_frequency_power', 'batches'}, true);

% raw_powers  = cell(1, length(probe_ids));

for ii = 1 : length(probe_ids)
    
    rec = ctl.load.spikeglx_ap_recording(probe_ids{ii});
    
    hf_power = HighFrequencyPowerProfile(rec, 0);
%     hf_power.compute_power_on_all_channels();
    
%     raw_powers{ii}  = hf_power.power_raw;
    
    probe_track     = ctl.load_track(probe_ids{ii}, 0);
    [boundaries, ~, region_str] = probe_track.region_boundaries();
    
    h_ax = {};
    figure('position', [150, 140, 1400, 850]);
    for jj = 1 : min(length(hf_power.sample_start_t), 20)
        h_ax{jj} = subplot(4, 5, jj);
        from_tip = hf_power.rec.channel_from_tip_um(0:200);
        idx = ismember(0:200, hf_power.rec.reference_channel_ids);
        power = raw_powers{ii}(:, jj);
        power(idx) = nan;
        plot(from_tip, power, 'color', 'k', 'linewidth', 0.75);
        for i = 1 : length(boundaries)-1
            % boundary
            line(boundaries(i)*[1, 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--', 'linewidth', 0.25);
            % text in middle of boundary
            if strcmp(region_str{i}, 'VISp5')
                line([1, 1]*sum(boundaries([i, i+1]))/2, get(gca, 'ylim'), 'color', 'b');
                text(sum(boundaries([i, i+1]))/2, max(get(gca, 'ylim')), region_str{i}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
            end
        end
        box off
        set(gca, 'plotboxaspectratio', [2, 1, 1]);
        title(sprintf('%i-%is', hf_power.sample_start_t(jj), hf_power.sample_start_t(jj) + hf_power.batch_duration_s));
    end
    
    FigureTitle(gcf, sprintf('%s', probe_ids{ii}))
    set(gcf, 'position', [16, 80, 1440, 890], 'papersize', [50, 40], 'renderer', 'painters');
    
    ctl.figs.save_fig(sprintf('%s.pdf', probe_ids{ii}), true);
end


%% overlay batches

for ii = 1 : length(probe_ids)
    
    figure('position', [150, 140, 1400, 850]);
    hold on;
    
    probe_track     = ctl.load_track(probe_ids{ii}, 0);
    [boundaries, ~, region_str] = probe_track.region_boundaries();
    
    
    for jj = 1 : min(size(raw_powers{ii}, 2), 20)
        
        from_tip = hf_power.rec.channel_from_tip_um(0:200);
        idx = ismember(0:200, hf_power.rec.reference_channel_ids);
        power = raw_powers{ii}(:, jj);
        power(idx) = nan;
        plot(from_tip, power, 'linewidth', 0.75);
    end
    
    plot(
    
    for i = 1 : length(boundaries)-1
        % boundary
        line(boundaries(i)*[1, 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--', 'linewidth', 0.25);
        % text in middle of boundary
        if strcmp(region_str{i}, 'VISp5')
            line([1, 1]*sum(boundaries([i, i+1]))/2, get(gca, 'ylim'), 'color', 'b');
            text(sum(boundaries([i, i+1]))/2, max(get(gca, 'ylim')), region_str{i}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
        end
    end
    box off
    set(gca, 'plotboxaspectratio', [2, 1, 1]);
    FigureTitle(gcf, sprintf('%s', probe_ids{ii}))
end

%% replot with bad channels removed for
bad_channel_ids = [26, 71, 97, 104, 123, 151, 155, 156, 167];
old_probe_ids = {'CA_176_1_rec1_rec2_rec3';
                 'CA_176_3_rec1_rec2_rec3'};

for ii = 1 : 2
    
     rec = ctl.load.spikeglx_ap_recording(probe_ids{ii});
    
    hf_power = HighFrequencyPowerProfile(rec, 0);

    probe_idx       = find(strcmp(probe_ids, old_probe_ids{ii}));
    probe_track     = ctl.load_track(old_probe_ids{ii}, 0);
    [boundaries, ~, region_str] = probe_track.region_boundaries();
    
    figure('position', [150, 140, 1400, 850]);
    
    for jj = 1 : min(length(hf_power.sample_start_t), 20)
        
        subplot(4, 5, jj);
        
        from_tip = hf_power.rec.channel_from_tip_um(0:200);
        idx = ismember(0:200, [hf_power.rec.reference_channel_ids; bad_channel_ids(:)]);
        power = raw_powers{probe_idx}(:, jj);
        power(idx) = nan;
        plot(from_tip, power, 'color', 'k', 'linewidth', 0.75);
        for i = 1 : length(boundaries)-1
            % boundary
            line(boundaries(i)*[1, 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--', 'linewidth', 0.25);
            % text in middle of boundary
            if strcmp(region_str{i}, 'VISp5')
                line([1, 1]*sum(boundaries([i, i+1]))/2, get(gca, 'ylim'), 'color', 'b');
                text(sum(boundaries([i, i+1]))/2, max(get(gca, 'ylim')), region_str{i}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
            end
        end
        
        box off
        set(gca, 'plotboxaspectratio', [2, 1, 1]);
        title(sprintf('%i-%is', hf_power.sample_start_t(jj), hf_power.sample_start_t(jj) + hf_power.sample_duration));
        
    end
    
    FigureTitle(gcf, sprintf('%s', old_probe_ids{ii}))
    set(gcf, 'position', [16, 80, 1440, 890], 'papersize', [50, 40], 'renderer', 'painters');
    
    ctl.figs.save_fig(sprintf('%s_bad_chans_removed.pdf', old_probe_ids{ii}), true);
end


%% now plot averages based on the "good" batches
batches_to_use = {1:13,
                  1:17, 
                  1:14,
                  3:19,
                  1:17,
                  1:8,
                  1:20,
                  [2:4, 6:8, 9:10],
                  [3:5, 8:10, 11:15],
                  7:14,
                  7:10,
                  [2:3, 5:10],
                  [1:5, 8:13],
                  1:18,
                  1:17,
                  [1:7, 11:16],
                  1:12,
                  1:13,
                  1:11,
                  1:11};

ctl.setup_figures({'high_frequency_power', 'avg'}, true);

for ii = 1 : length(probe_ids)
    
     rec = ctl.load.spikeglx_ap_recording(probe_ids{ii});
    
    hf_power = HighFrequencyPowerProfile(rec, 0);

    probe_track     = ctl.load_track(probe_ids{ii}, 0);
    [boundaries, ~, region_str] = probe_track.region_boundaries();
    
    figure();
    from_tip = hf_power.rec.channel_from_tip_um(0:200);
    
    if ismember(probe_ids{ii}, old_probe_ids)
        idx = ismember(0:200, [hf_power.rec.reference_channel_ids; bad_channel_ids(:)]);
    else
        idx = ismember(0:200, hf_power.rec.reference_channel_ids);
    end
    
    norm_power = raw_powers{ii}(:, batches_to_use{ii});
    norm_power = bsxfun(@minus, norm_power, mean(norm_power, 1));
    norm_power = bsxfun(@rdivide, norm_power, std(norm_power, [], 1));
    
    power = mean(norm_power, 2);
    power(idx) = nan;
    
    plot(from_tip, power, 'color', 'k', 'linewidth', 0.75);
    
    for jj = 1 : length(boundaries)-1
        % boundary
        line(boundaries(jj)*[1, 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--', 'linewidth', 0.25);
        text(sum(boundaries([jj, jj+1]))/2, max(get(gca, 'ylim')), region_str{jj}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
        % text in middle of boundary
        if strcmp(region_str{jj}, 'VISp5')
            line([1, 1]*sum(boundaries([jj, jj+1]))/2, get(gca, 'ylim'), 'color', 'b');
        end
    end
    
    box off
    set(gca, 'plotboxaspectratio', [2, 1, 1]);
    
    FigureTitle(gcf, sprintf('%s', probe_ids{ii}))
    set(gcf, 'position', [16, 80, 1440, 890], 'papersize', [50, 40], 'renderer', 'painters');
    
    ctl.figs.save_fig(sprintf('%s.pdf', probe_ids{ii}), true);
end   



%% only plot along banks
ctl.setup_figures({'high_frequency_power', 'avg_banks2'}, true);

for ii = 1 : length(probe_ids)
    
    rec = ctl.load.spikeglx_ap_recording(probe_ids{ii});
    
    hf_power = HighFrequencyPowerProfile(rec, 0);
    
    probe_track     = ctl.load_track(probe_ids{ii}, 0);
    [boundaries, ~, region_str] = probe_track.region_boundaries();
    
    figure();
    from_tip = hf_power.rec.channel_from_tip_um(0:200);
    
    if ismember(probe_ids{ii}, old_probe_ids)
        idx = ismember(0:200, [hf_power.rec.reference_channel_ids; bad_channel_ids(:)]);
    else
        idx = ismember(0:200, hf_power.rec.reference_channel_ids);
    end
    
    norm_power = raw_powers{ii}(:, batches_to_use{ii});
    norm_power = bsxfun(@minus, norm_power, mean(norm_power, 1));
    norm_power = bsxfun(@rdivide, norm_power, std(norm_power, [], 1));
    
    power = mean(norm_power, 2);
    power(idx) = nan;
    
    plot(from_tip(1:4:end), power(1:4:end), 'linewidth', 0.75); hold on;
    plot(from_tip(2:4:end), power(2:4:end), 'linewidth', 0.75); hold on;
    plot(from_tip(3:4:end), power(3:4:end), 'linewidth', 0.75); hold on;
    plot(from_tip(4:4:end), power(4:4:end), 'linewidth', 0.75);
    
    for jj = 1 : length(boundaries)-1
        % boundary
        line(boundaries(jj)*[1, 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--', 'linewidth', 0.25);
        text(sum(boundaries([jj, jj+1]))/2, max(get(gca, 'ylim')), region_str{jj}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
        % text in middle of boundary
        if strcmp(region_str{jj}, 'VISp5')
            line([1, 1]*sum(boundaries([jj, jj+1]))/2, get(gca, 'ylim'), 'color', 'b');
        end
    end
    
    box off
    set(gca, 'plotboxaspectratio', [2, 1, 1]);
    
    FigureTitle(gcf, sprintf('%s', probe_ids{ii}))
    set(gcf, 'position', [16, 80, 1440, 890], 'papersize', [50, 40], 'renderer', 'painters');
    
    ctl.figs.save_fig(sprintf('%s.pdf', probe_ids{ii}), true);
end 



%% interpolate along banks and average
ctl.setup_figures({'high_frequency_power', 'avg_interp'}, true);

for ii = 1 : length(probe_ids)
    
    rec = ctl.load.spikeglx_ap_recording(probe_ids{ii});
    hf_power = HighFrequencyPowerProfile(rec, 0);

    probe_track     = ctl.load_track(probe_ids{ii}, 0);
    [boundaries, ~, region_str] = probe_track.region_boundaries();
    
    figure();
    from_tip = hf_power.rec.channel_from_tip_um(0:200);
    
    if ismember(probe_ids{ii}, old_probe_ids)
        idx = ismember(0:200, [hf_power.rec.reference_channel_ids; bad_channel_ids(:)]);
    else
        idx = ismember(0:200, hf_power.rec.reference_channel_ids);
    end
    
    norm_power = raw_powers{ii}(:, batches_to_use{ii});
    norm_power = bsxfun(@minus, norm_power, mean(norm_power, 1));
    norm_power = bsxfun(@rdivide, norm_power, std(norm_power, [], 1));
    
    power = mean(norm_power, 2);
    power(idx) = nan;
    
    tip1 = from_tip(1:4:end);
    tip2 = from_tip(2:4:end);
    tip3 = from_tip(3:4:end);
    tip4 = from_tip(4:4:end);
    
    power1 = power(1:4:end);
    power2 = power(2:4:end);
    power3 = power(3:4:end);
    power4 = power(4:4:end);
    
    % interpolate
    xnew = from_tip(1):from_tip(end);
    power1 = interp1(tip1(~isnan(power1)), power1(~isnan(power1)), xnew);
    power2 = interp1(tip2(~isnan(power2)), power2(~isnan(power2)), xnew);
    power3 = interp1(tip3(~isnan(power3)), power3(~isnan(power3)), xnew);
    power4 = interp1(tip4(~isnan(power4)), power4(~isnan(power4)), xnew);
    
    plot(xnew, power1, 'linewidth', 0.75); hold on;
    plot(xnew, power2, 'linewidth', 0.75); hold on;
    plot(xnew, power3, 'linewidth', 0.75); hold on;
    plot(xnew, power4, 'linewidth', 0.75);
    
    avg_power = (power1+power2+power3+power4)/4; 
    plot(xnew, avg_power, 'color', [0.5, 0.5, 0.5], 'linewidth', 1);
    
    avg_power_smooth = smooth(avg_power, 100);
    plot(xnew, avg_power_smooth, 'color', 'k', 'linewidth', 1);
    
    for jj = 1 : length(boundaries)-1
        % boundary
        line(boundaries(jj)*[1, 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--', 'linewidth', 0.25);
        text(sum(boundaries([jj, jj+1]))/2, max(get(gca, 'ylim')), region_str{jj}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
        % text in middle of boundary
        if strcmp(region_str{jj}, 'VISp5')
            line([1, 1]*sum(boundaries([jj, jj+1]))/2, get(gca, 'ylim'), 'color', 'b', 'linewidth', 1);
        end
    end
    
    box off
    set(gca, 'plotboxaspectratio', [2, 1, 1]);
    
    FigureTitle(gcf, sprintf('%s', probe_ids{ii}))
    set(gcf, 'position', [16, 80, 1440, 890], 'papersize', [50, 40], 'renderer', 'painters');
    
    ctl.figs.save_fig(sprintf('%s.pdf', probe_ids{ii}), true);
end






%% boundaries shifted
offset = [-16
110
-69
-65
-95
-47
-44
24
-57
-62
55
40
125
164
42
75
48
212
33
132];

ctl.setup_figures({'high_frequency_power', 'avg_shifted'}, true);

for ii = 1 : length(probe_ids)
    
    rec = ctl.load.spikeglx_ap_recording(probe_ids{ii});
    hf_power = HighFrequencyPowerProfile(rec, 0);
    
    figure();
    probe_track     = ctl.load_track(probe_ids{ii}, 0);
    [boundaries, ~, region_str] = probe_track.region_boundaries();
    
    from_tip = hf_power.rec.channel_from_tip_um(0:200);
    
    if ismember(probe_ids{ii}, old_probe_ids)
        idx = ismember(0:200, [hf_power.rec.reference_channel_ids; bad_channel_ids(:)]);
    else
        idx = ismember(0:200, hf_power.rec.reference_channel_ids);
    end
    
    norm_power = raw_powers{ii}(:, batches_to_use{ii});
%     norm_power = bsxfun(@minus, norm_power, mean(norm_power, 1));
%     norm_power = bsxfun(@rdivide, norm_power, std(norm_power, [], 1));
    
    power = mean(norm_power, 2);
    power(idx) = nan;
    
    tip1 = from_tip(1:4:end);
    tip2 = from_tip(2:4:end);
    tip3 = from_tip(3:4:end);
    tip4 = from_tip(4:4:end);
    
    power1 = power(1:4:end);
    power2 = power(2:4:end);
    power3 = power(3:4:end);
    power4 = power(4:4:end);
    
    % interpolate
    xnew = from_tip(1):from_tip(end);
    power1 = interp1(tip1(~isnan(power1)), power1(~isnan(power1)), xnew);
    power2 = interp1(tip2(~isnan(power2)), power2(~isnan(power2)), xnew);
    power3 = interp1(tip3(~isnan(power3)), power3(~isnan(power3)), xnew);
    power4 = interp1(tip4(~isnan(power4)), power4(~isnan(power4)), xnew);
    
    plot(xnew, power1, 'linewidth', 0.75); hold on;
    plot(xnew, power2, 'linewidth', 0.75); hold on;
    plot(xnew, power3, 'linewidth', 0.75); hold on;
    plot(xnew, power4, 'linewidth', 0.75);
    
    avg_power = (power1+power2+power3+power4)/4; 
    plot(xnew, avg_power, 'color', [0.5, 0.5, 0.5], 'linewidth', 1);
    
    avg_power_smooth = smooth(avg_power, 100);
    plot(xnew, avg_power_smooth, 'color', 'k', 'linewidth', 1);
    
    for jj = 1 : length(boundaries)-1
        % boundary
        line(boundaries(jj)*[1, 1]-offset(ii), get(gca, 'ylim'), 'color', 'k', 'linestyle', '--', 'linewidth', 0.25);
        text(sum(boundaries([jj, jj+1]))/2-offset(ii), max(get(gca, 'ylim')), region_str{jj}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
        % text in middle of boundary
        if strcmp(region_str{jj}, 'VISp5')
            line([1, 1]*sum(boundaries([jj, jj+1]))/2-offset(ii), get(gca, 'ylim'), 'color', 'b', 'linewidth', 1);
        end
    end
    
    box off
    set(gca, 'plotboxaspectratio', [2, 1, 1]);
%     title(gca, sprintf('%s', probe_ids{ii}));
    FigureTitle(gcf, sprintf('%s', probe_ids{ii}))
    set(gcf, 'position', [16, 80, 1440, 890], 'papersize', [50, 40], 'renderer', 'painters');
    ctl.figs.save_fig('all_probes_shifted', true);
    
end


