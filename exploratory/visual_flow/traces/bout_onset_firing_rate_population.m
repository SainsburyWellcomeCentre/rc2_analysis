% plot rasters for an animal
clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'visual_flow';  % 'darkness', 'visual_flow', 'head_tilt', 'mismatch_nov20'
combination = 'protocols';

config = RC2AnalysisConfig();
figs = RC2Figures(config);

% where to save figure
figs.save_on = true;
figs.set_figure_subdir('visual_flow', 'rasters');

% get details of the experiment
probe_fnames = experiment_details(experiment, combination);

% number of seconds before and after the bout onset to take
prepad              = -1;
postpad             = 1;

% minimum duration of a bout
min_bout_duration   = 2;

% default options
options = default_options();

% fundamental trial types
trial_types = {'Coupled', 'EncoderOnly'};
replay_types = {'StageOnly', 'ReplayOnly'};
protocols = 1:6;
title_str = {'MVT', 'MV', 'VT (MVT)', 'VT (MV)', 'V (MVT)', 'V (MV)'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% store ids
probe_id_store = [];
cluster_id_store = [];

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    
    % sampling rate of the recordings
    fs = data.sample_rate;
    
    % traces will be on a common timebase
    common_t = linspace(prepad, postpad, round((postpad-prepad)*fs));
    
    % object to handle visual flow
    vf = VisualFlowExperiment(data, config);
     
    % filter clusters according to the cluster ID list file
    clusters = data.VISp_clusters;
    
    bouts = [];
    bout_protocol{probe_i} = [];
    
    for type_i = 1 : length(trial_types)
        
        % get trials for this type
        these_trials = vf.trials_of_type(trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            
            bouts_t = get_motion_bouts_by_trial(these_trials(trial_i), options.stationary, min_bout_duration);
            
            % append bouts
            bouts = [bouts, bouts_t];
            if strcmp(trial_types{type_i}, 'Coupled')
                bout_protocol{probe_i} = [bout_protocol{probe_i}; ones(length(bouts_t), 1)];
            elseif strcmp(trial_types{type_i}, 'EncoderOnly')
                bout_protocol{probe_i} = [bout_protocol{probe_i}; 2*ones(length(bouts_t), 1)];
            end
            
            
            for rep_i = 1 : length(replay_types)
                
                replay_trial = vf.get_replay_of(these_trials(trial_i), replay_types{rep_i});
                assert(length(replay_trial) == 1, 'More replays than expected');
                
                offset = vf.get_offset(these_trials(trial_i), replay_trial);
                
                clear bouts_rep_t
                
                for bout_i = 1 : length(bouts_t)
                    
                    s = bouts_t(bout_i).start_idx + offset;
                    e = bouts_t(bout_i).end_idx + offset;
                    bouts_rep_t(bout_i) = MotionBout(s, e, replay_trial);
                    
                end
                
                % append
                bouts = [bouts, bouts_rep_t];
                
                if strcmp(trial_types{type_i}, 'Coupled') && strcmp(replay_types{rep_i}, 'StageOnly')
                    bout_protocol{probe_i} = [bout_protocol{probe_i}; 3*ones(length(bouts_t), 1)];
                elseif strcmp(trial_types{type_i}, 'EncoderOnly') && strcmp(replay_types{rep_i}, 'StageOnly')
                    bout_protocol{probe_i} = [bout_protocol{probe_i}; 4*ones(length(bouts_t), 1)];
                elseif strcmp(trial_types{type_i}, 'Coupled') && strcmp(replay_types{rep_i}, 'ReplayOnly')
                    bout_protocol{probe_i} = [bout_protocol{probe_i}; 5*ones(length(bouts_t), 1)];
                elseif strcmp(trial_types{type_i}, 'EncoderOnly') && strcmp(replay_types{rep_i}, 'ReplayOnly')
                    bout_protocol{probe_i} = [bout_protocol{probe_i}; 6*ones(length(bouts_t), 1)];
                end
                
            end
            
        end
        
    end
    
    % start times of the bouts
    start_t = [bouts(:).start_time];
    
    % use bout timings to get raster data
    [spike_rates{probe_i}, spike_times] = ...
        get_raster_data(clusters, start_t, common_t, fs, options.spiking);
    
    % get velocity traces around those start times
    velocity_traces{probe_i} = nan(length(prepad*fs : postpad*fs-1), length(bouts));
    
    for bout_i = 1 : length(bouts)
        
        idx = bouts(bout_i).start_idx + (prepad*fs : postpad*fs-1);
        velocity_traces{probe_i}(:, bout_i) = bouts(bout_i).trial.velocity(idx);
        
    end
    
    probe_id_store = [probe_id_store, probe_i * ones(1, length(clusters))];
    cluster_id_store = [cluster_id_store, [clusters(:).id]];
end


velocity_traces_cat = cat(2, velocity_traces);


%%
figure

h_ax = [];

yL = [-2, -inf];

for prot_i = 1 : 6
    
    mean_fr_each_cluster = [];

    vel_traces = [];
    
    for probe_i = 1 : 4
        
        idx = bout_protocol{probe_i} == prot_i;
        
        mean_fr_each_cluster = ...
            [mean_fr_each_cluster, ...
            squeeze(mean(spike_rates{probe_i}(:, idx, :), 2))];
        
        vel_traces = [vel_traces, velocity_traces{probe_i}(:, idx)];
        
    end
    
    % subtract baseline
    baseline_idx = common_t > -1 & common_t < 0;
    
    mean_fr_each_cluster = ...
        bsxfun(@minus, mean_fr_each_cluster, ...
                       mean(mean_fr_each_cluster(baseline_idx, :), 1));
                   
    mean_fr = mean(mean_fr_each_cluster, 2);
    
    sem_fr = std(mean_fr_each_cluster, [], 2)/sqrt(size(mean_fr_each_cluster, 2));
    
    mean_vel = mean(vel_traces, 2);
    
    sem_vel = std(vel_traces, [], 2)/sqrt(size(vel_traces, 2));
    
    h_ax(prot_i, 1) = subplot(3, 2, prot_i);
    
    hold on;
    
    axis normal;
    
    set(gca, 'activepositionproperty', 'position');
    
    plot(common_t, mean_vel, 'color', 'r', 'linewidth', 2);
    
    plot(common_t, mean_vel + sem_vel, 'color', [1, 0.7, 0.7]);
    
    plot(common_t, mean_vel - sem_vel, 'color', [1, 0.7, 0.7]);
    
    box off
    
    set(gca, 'plotboxaspectratio', [3, 1, 1], 'clipping', 'off');
    
    set(gca, 'yaxislocation', 'right', 'ycolor', 'r');
    
    ylabel('cm/s')
    
    pos = get(h_ax(prot_i, 1), 'position');
    
    h_ax(prot_i, 2) = axes('position', pos);
    
    hold on;
    
    axis normal;
    
    set(gca, 'activepositionproperty', 'position');
    
    set(h_ax(prot_i, 2), 'color', 'none');
    
    plot(common_t, mean_fr, 'color', 'k', 'linewidth', 2);
    
    plot(common_t, mean_fr - sem_fr, 'color', [0.7, 0.7, 0.7]);
    
    plot(common_t, mean_fr + sem_fr, 'color', [0.7, 0.7, 0.7]);
    
    yL(1) = min(yL(1), min(get(gca, 'ylim')));
    yL(2) = max(yL(2), max(get(gca, 'ylim')));
    
    box off
    
    set(h_ax(prot_i, 2), 'plotboxaspectratio', [3, 1, 1], 'clipping', 'off');
    
    title(title_str{prot_i});
    
    ylabel('\Delta FR (Hz)');
    
    xlabel('Time (s)');
    
end

set(gcf, 'renderer', 'painters')


for prot_i = 1 : 6
    
    set(gcf, 'currentaxes', h_ax(prot_i, 2));
    
    set(h_ax(prot_i, 1), 'ylim', [-2, 10]);
    set(h_ax(prot_i, 2), 'ylim', yL);
    
    line([0, 0], yL, 'color', 'k', 'linestyle', '--');

end



%% Heatmaps
load('r2bcmap', 'map');

cmap = map;

figure;

h_ax = [];

cl = [-8, 8];

for prot_i = 1 : 6
    
    mean_fr_each_cluster = [];

    vel_traces = [];
    
    for probe_i = 1 : 4
        
        idx = bout_protocol{probe_i} == prot_i;
        
        mean_fr_each_cluster = ...
            [mean_fr_each_cluster, ...
            squeeze(mean(spike_rates{probe_i}(:, idx, :), 2))];
        
    end
    
    % subtract baseline
    baseline_idx = common_t > -1 & common_t < 0;
    
    mean_fr_each_cluster = ...
        bsxfun(@minus, mean_fr_each_cluster, ...
                       mean(mean_fr_each_cluster(baseline_idx, :), 1));
    
    if prot_i == 1
        
        baseline_idx = common_t > -1 & common_t <= 0;
        response_idx = common_t > 0 & common_t < 1;
        
        resp = mean(mean_fr_each_cluster(response_idx, :), 1) - ...
               mean(mean_fr_each_cluster(baseline_idx, :), 1);
        
        [~, order_idx] = sort(resp, 'ascend');
    end
                   
    h_ax(prot_i) = subplot(3, 2, prot_i);
    
    hold on;
    
    axis normal;
    
    set(gca, 'activepositionproperty', 'position');
    
    h_im = imagesc(mean_fr_each_cluster(:, order_idx)');
    
    set(h_im, 'xdata', common_t, 'ydata', size(mean_fr_each_cluster, 2):-1:1);
    
    set(gca, 'clim', cl, 'ydir', 'reverse');
    
    set(gca, 'plotboxaspectratio', [1, 2, 1], 'clipping', 'off');
    
    box off
    
    title(title_str{prot_i});
    
    ylabel('Cluster #');
    
    colormap(cmap);
    
    if prot_i == 6
        xlabel('Time (s)');
        pos = get(h_ax(prot_i), 'position');
        ax = axes('position', pos);
        set(ax, 'clim', cl);
        h = colorbar;
        set(get(h, 'label'), 'string', '\DeltaFR (Hz)');
        set(ax, 'visible', 'off');
    end
    
end

set(gcf, 'renderer', 'painters')


for prot_i = 1 : 6
    
    set(gcf, 'currentaxes', h_ax(prot_i));
    
    line([0, 0], [0, 37], 'color', 'k', 'linestyle', '--');

end

