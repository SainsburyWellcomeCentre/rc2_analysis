probe_fname         = {'CA_176_1_rec1_rec2_rec3', ...
    'CA_176_3_rec1_rec2_rec3', ...
    'CAA-1110264_rec1_rec2', ...
    'CAA-1110264_rec1_rec2', ...
    'CAA-1110264_rec1_rec2', ...
    'CAA-1110264_rec1_rec2', ...
    'CA_176_3_rec1_rec2_rec3', ...
    'CA_176_3_rec1_rec2_rec3'};

cluster_id          = [186, ...
    275, ...
    209, ...
    209, ...
    209, ...
    209, ...
    274, ...
    274];

trial_type          = {'StageOnly', ...
    'StageOnly', ...
    'ReplayOnly', ...
    'StageOnly', ...
    'EncoderOnly', ...
    'Coupled', ...
    'EncoderOnly', ...
    'Coupled'};


cl = {[0, 30], [0, 60], [0, 70], [0, 70], [0, 70], [0, 70], [0, 70], [0, 70]};

title_code          = {'T', 'T', 'V', 'VT', 'MV', 'MVT', 'M', 'MT'};

prepad              = -0.2;
postpad             = 1;

trace_type          = 'heatmap';  % 'heatmap', 'average'

config              = RC2AnalysisConfig();
figs                = RC2Figures(config);
figs.save_on        = true;
figs.set_figure_subdir('examples');
options             = default_options();

for i = 1 : length(probe_fname)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    if ismember(i, [1, 2, 7, 8])
        exp_obj             = DarknessExperiment(data, config);
    else
        exp_obj             = VisualFlowExperiment(data, config);
    end
    cluster             = get_cluster_by_id(data.selected_clusters, cluster_id(i));
    trials              = exp_obj.trials_of_type(trial_type{i});
    bouts               = [];
    
    for trial_i = 1 : length(trials)
        
        replayed_trial = exp_obj.get_replayed_trial(trials(trial_i));
        
        if ~isempty(replayed_trial)
            
            offset = exp_obj.get_offset(replayed_trial, trials(trial_i));
            bouts_tt = replayed_trial.motion_bouts(true);
            
            clear bouts_t
            
            for bout_i = 1 : length(bouts_tt)
                s = bouts_tt(bout_i).start_idx + offset;
                e = bouts_tt(bout_i).end_idx + offset;
                bouts_t(bout_i) = MotionBout(s, e, trials(trial_i));
            end
        else
            
            bouts_t = trials(trial_i).motion_bouts(true);
        end
        
        bouts = [bouts, bouts_t];
    end
    
    
    %%
    fs = data.sample_rate;
    
    % traces will be on a common timebase
    common_t = linspace(prepad, postpad, round((postpad - prepad) * fs));
    
    %
    start_t = [bouts(:).start_time];
    
    % use bout timings to get raster data
    [spike_rates, spike_times] = ...
        get_raster_data(cluster, start_t, common_t, fs, options.spiking);
    
    % get velocity traces around those start times
    M_traces = nan(length(prepad*fs : postpad*fs-1), length(bouts));
    V_traces = nan(length(prepad*fs : postpad*fs-1), length(bouts));
    T_traces = nan(length(prepad*fs : postpad*fs-1), length(bouts));
    
    for bout_i = 1 : length(bouts)
        
        idx = bouts(bout_i).start_idx + (prepad*fs : postpad*fs-1);
        
        M_traces(:, bout_i) = bouts(bout_i).trial.treadmill_speed(idx);
        if isempty(bouts(bout_i).trial.multiplexer_output)
            V_traces(:, bout_i) = zeros(size(V_traces, 1), 1);
        else
            V_traces(:, bout_i) = bouts(bout_i).trial.multiplexer_output(idx);
        end
        T_traces(:, bout_i) = bouts(bout_i).trial.stage(idx);
    end
    
    
    all_traces = {M_traces, V_traces, T_traces};
    label = {'M', 'V', 'T'};
    
    
    %%
    
    h_fig = figs.a4figure();
    
    h_ax = [];
    M = -inf;
    m = inf;
    
    for trace_i = 1 : 3
        
        h_ax(trace_i) = subplot(4, 1, trace_i);
        hold on;
        
        vel_mean = mean(all_traces{trace_i}, 2);
        vel_sem = std(all_traces{trace_i}, [], 2) / ...
            sqrt(size(all_traces{trace_i}, 2));
        
        upper = vel_mean + vel_sem;
        lower = vel_mean - vel_sem;
        
        plot(common_t, upper, 'color', [0.8, 0.5, 0.5], 'linewidth', 1);
        plot(common_t, lower, 'color', [0.8, 0.5, 0.5], 'linewidth', 1);
        plot(common_t, vel_mean, 'r', 'linewidth', 2);
        
        box off
        set(gca, 'plotboxaspectratio', [3, 1, 1], 'clipping', 'off', ...
            'xlim', [prepad, postpad]);
        if trace_i == 1
            title(sprintf('%s, Cluster %i, %s', probe_fname{i}, cluster_id(i), title_code{i}), ...
                'interpreter', 'none');
        end
        
        ylabel('Avg. speed (cm/s)');
        
        m = min(m, min(get(gca, 'ylim')));
        M = max(M, max(get(gca, 'ylim')));
        
    end
    
    for trace_i = 1 : 3
        set(gcf, 'currentaxes', h_ax(trace_i));
        set(h_ax(trace_i), 'ylim', [m, M]);
        text(0.2, ...
            max(get(h_ax(trace_i), 'ylim')), label{trace_i}, ...
            'horizontalalignment', 'center', 'verticalalignment', 'top', ...
            'fontsize', 24);
        line([0, 0], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--');
    end
    
    
    subplot(4, 1, 4);
    hold on;
    
    
    if strcmp(trace_type, 'heatmap')
        
        load('r2bcmap.mat', 'map');
        
        h_im = imagesc(spike_rates');
        h_ax = gca;
        colormap(map);
        if isempty(cl{i})
            tcl = get(h_ax, 'clim');
        else
            tcl = cl{i};
        end
        set(h_im, 'xdata', common_t);
        set(h_ax, 'plotboxaspectratio', [3, 1, 1], ...
            'xlim', [prepad, postpad], 'clim', max(tcl)*[-1, 1]);
        box off;
        
        pos = get(h_ax, 'position');
        ax = axes('position', pos+[0, 0, 0.1, 0]);
        set(ax, 'plotboxaspectratio', [3, 1, 1], ...
            'xlim', [prepad, postpad], 'clim', max(tcl)*[-1, 1]);
        box off;
        h = colorbar;
        set(h, 'ytick', [-tcl(2), 0, tcl(2)]);
        set(get(h, 'label'), 'string', 'FR (Hz)');
        set(ax, 'visible', 'off');
        set(gcf, 'currentaxes', h_ax);
        ylabel('Bout #');
        
    elseif strcmp(trace_type, 'average')
        
        fr_mean = mean(spike_rates, 2);
        fr_sem = std(spike_rates, [], 2)/sqrt(size(spike_rates, 2));
        
        upper = fr_mean + fr_sem;
        lower = fr_mean - fr_sem;
        
        plot(common_t, upper, 'color', [0.7, 0.7, 0.7], 'linewidth', 1);
        plot(common_t, lower, 'color', [0.7, 0.7, 0.7], 'linewidth', 1);
        plot(common_t, fr_mean, 'k', 'linewidth', 2);
        
        set(gca, 'plotboxaspectratio', [3, 1, 1], 'clipping', 'off', ...
            'xlim', [prepad, postpad]);
        ylabel('Avg. firing rate (Hz)');
    end
    
    line([0, 0],  get(gca, 'ylim'), 'color', 'k', 'linestyle', '--');
    
    
    xlabel('Time (s)');
    figs.save_fig(sprintf('%s_cluster_%03i_%s.pdf', probe_fname{i}, cluster_id(i), title_code{i}))
end



