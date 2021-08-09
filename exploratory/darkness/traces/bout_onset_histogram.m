experiment = 'visual_flow';

config = RC2AnalysisConfig();

figs = RC2Figures(config);
figs.save_on = true;
figs.set_figure_subdir(experiment, 'averaged_responses_and_heatmaps');

% get details of the experiment
probe_fnames = experiment_details(experiment, 'protocols');

% number of seconds before and after the bout onset to take
prepad              = -1;
postpad             = 1;

save_fname = sprintf( ...
    'darkness_average_traces_motion_onset_%is_baseline_%is_motion.pdf', ...
    -prepad, postpad);

% default options
options = default_options();

% fundamental trial types
if strcmp(experiment, 'visual_flow')
    trial_types = ...
                    {'Coupled', 'EncoderOnly', 'StageOnly', ...
                     'StageOnly', 'ReplayOnly', 'ReplayOnly'};
    replay_types =  {'', '', 'Coupled', ...
                     'EncoderOnly', 'Coupled', 'EncoderOnly'};
    protocol_ids = 1 : 6;
    title_str =     {'MVT', 'MV', 'VT (MVT)', ...
                     'VT (MV)', 'V (MVT)', 'V (MV)'};
else
    trial_types = {'Coupled', 'EncoderOnly', 'StageOnly'};
    replay_types = {'', '', ''};
    protocol_ids = 1:3;
    title_str = {'MT', 'M', 'T'};
end

spike_rates     = cell(1, length(trial_types));
psth            = cell(1, length(trial_types));
velocity_traces = cell(1, length(trial_types));
response_type   = cell(1, length(trial_types));

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.VISp_clusters;
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(data, config);
    end
    
    % sampling rate of the recordings
    fs = data.sample_rate;
    
    % traces will be on a common timebase
    common_t = linspace(prepad, postpad, round((postpad-prepad)*fs));
    
    for type_i = 1 : length(trial_types)
        
        % get trials for this type
        trials = exp_obj.trials_of_type_replay_of_type(trial_types{type_i}, replay_types{type_i});
        
        bouts = [];
        
        for trial_i = 1 : length(trials)
            
            atrial = exp_obj.to_aligned(trials(trial_i));
            bouts = [bouts, atrial.motion_bouts(true)];
        end
        
        % start times of the bouts
        start_t = [bouts(:).start_time];
        
        % use bout timings to get raster data
        % n_sample_points x 
        spike_rates{type_i}{probe_i} = get_raster_data(clusters, start_t, common_t, fs, options.spiking);
        
        for cluster_i = 1 : length(clusters)
            
            fr = FiringRate(clusters(cluster_i).spike_times);
            [psth{type_i}{probe_i}{cluster_i}, ...
                bin_t] = fr.psth(start_t, 0.01, [-prepad, postpad]);
        end
        
        % get velocity traces around those start times
        velocity_traces{type_i}{probe_i} = nan(length(prepad*fs : postpad*fs-1), length(bouts));
        
        for bout_i = 1 : length(bouts)
            
            idx = bouts(bout_i).start_idx + (prepad*fs : postpad*fs-1);
            idx(idx > length(bouts(bout_i).trial.velocity)) = [];
            velocity_traces{type_i}{probe_i}(1:length(idx), bout_i) = bouts(bout_i).trial.velocity(idx);
        end
        
        response_type{type_i}{probe_i} = [];
        
        % get significant responders
        for cluster_i = 1 : length(clusters)
            
            prot_i = protocol_ids(type_i);
            
            x = exp_obj.trial_stationary_fr(clusters(cluster_i).id, prot_i);
            y = exp_obj.trial_motion_fr(clusters(cluster_i).id, prot_i);
            
            if signrank(x, y) < 0.05 && nanmedian(x) <= nanmedian(y)
                response_type{type_i}{probe_i}(end+1) = 1;
            elseif signrank(x, y) < 0.05 && nanmedian(x) > nanmedian(y)
                response_type{type_i}{probe_i}(end+1) = -1;
            else
                response_type{type_i}{probe_i}(end+1) = 0;
            end
        end        
    end
end


velocity_traces_cat = cell(1, length(trial_types));
for type_i = 1 : length(trial_types)
    velocity_traces_cat{type_i} = cat(2, velocity_traces{type_i});
end




%%
h_fig = figs.a4figure();
h_ax = [];
yL = [-2, -inf];
resp_type = [1, 0, -1];
resp_cols = {[0, 128, 0]/255, [0, 0, 0], [0, 0, 128]/255};

for prot_i = 1 : length(trial_types)
    
    vel_traces = cat(2, velocity_traces_cat{prot_i}{:});
    
%     mean_fr_each_cluster = [];
    cat_psth = [];
    cat_response_type = [];
    
    for probe_i = 1 : 4
        
        cat_psth = cat(1, cat_psth, cat(1, psth{prot_i}{probe_i}{:}));
        
%         mean_fr_each_cluster = ...
%             [mean_fr_each_cluster, ...
%             squeeze(mean(spike_rates{prot_i}{probe_i}, 2))];
        cat_response_type = cat(2, cat_response_type, response_type{type_i}{probe_i});
    end
    
    % subtract baseline
    baseline_idx = common_t > -1 & common_t < 0;
    
    mean_fr = {};
    
    for resp_i = 1 : 3
        mean_fr{resp_i} = sum(cat_psth(cat_response_type == resp_type(resp_i), :), 1);
    end
    
    mean_vel = nanmean(vel_traces, 2);
    sem_vel = nanstd(vel_traces, [], 2)/sqrt(size(vel_traces, 2));
    
    
    h_ax(prot_i, 1) = subplot(3, 2, prot_i);
    hold on;
    axis normal;
    
    set(gca, 'activepositionproperty', 'position');
    
    common_t_sub = common_t(1:4:end);
    
    for resp_i = 1
        for bar_i = 1 : length(bin_t)-1
            patch(bin_t([bar_i, bar_i, bar_i+1, bar_i+1]), ...
                  [0, mean_fr{resp_i}([bar_i, bar_i]), 0], 'b');
        end
    end
    
    yL(1) = min(yL(1), min(get(gca, 'ylim')));
    yL(2) = max(yL(2), max(get(gca, 'ylim')));
    box off;
    set(h_ax(prot_i, 1), 'plotboxaspectratio', [3, 1, 1]);
    title(title_str{prot_i});
    if prot_i == 1
        ylabel('\Delta FR (Hz)');
    end
    xlabel('Time (s)');
    
    
    pos = get(h_ax(prot_i, 1), 'position');
    h_ax(prot_i, 2) = axes('position', pos);
    hold on;
    axis normal;
    
    set(gca, 'activepositionproperty', 'position');
    set(h_ax(prot_i, 2), 'color', 'none');
    
    plot(common_t, mean_vel, 'color', 'r', 'linewidth', 2);
    plot(common_t, mean_vel + sem_vel, 'color', [1, 0.7, 0.7]);
    plot(common_t, mean_vel - sem_vel, 'color', [1, 0.7, 0.7]);
    
    box off
    set(gca, 'plotboxaspectratio', [3, 1, 1]);
    set(gca, 'yaxislocation', 'right', 'ycolor', 'r');
    if prot_i == 1
        ylabel('cm/s')
    end
    
    
end

set(gcf, 'renderer', 'painters')


for prot_i = 1 : length(trial_types)
    
    set(gcf, 'currentaxes', h_ax(prot_i, 2));
    
    set(h_ax(prot_i, 2), 'ylim', [-2, 30]);
    set(h_ax(prot_i, 1), 'ylim', yL);
    
    set(h_ax(prot_i, 1), 'xlim', [-1, 4]);
    set(h_ax(prot_i, 2), 'xlim', [-1, 4]); % [prepad, postpad]
    
    line([0, 0], yL, 'color', 'k', 'linestyle', '--');
end


% figs.save_fig(save_fname);