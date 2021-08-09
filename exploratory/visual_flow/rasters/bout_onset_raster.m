% plot rasters for an animal
clear all
close all

import helper.*

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

for probe_i = 4%1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    
    % sampling rate of the recordings
    fs = data.sample_rate;
    
    % traces will be on a common timebase
    common_t = linspace(prepad, postpad, round((postpad-prepad)*fs));
    
    % object to handle visual flow
    vf = VisualFlowExperiment(data, config);
     
    % filter clusters according to the cluster ID list file
    clusters = data.selected_clusters;
    
    bouts = [];
    bout_protocol = [];
    
    for type_i = 1 : length(trial_types)
        
        % get trials for this type
        these_trials = vf.trials_of_type(trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            
            bouts_t = get_motion_bouts_by_trial(these_trials(trial_i), options.stationary, min_bout_duration);
            
            % append bouts
            bouts = [bouts, bouts_t];
            if strcmp(trial_types{type_i}, 'Coupled')
                bout_protocol = [bout_protocol; ones(length(bouts_t), 1)];
            elseif strcmp(trial_types{type_i}, 'EncoderOnly')
                bout_protocol = [bout_protocol; 2*ones(length(bouts_t), 1)];
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
                    bout_protocol = [bout_protocol; 3*ones(length(bouts_t), 1)];
                elseif strcmp(trial_types{type_i}, 'EncoderOnly') && strcmp(replay_types{rep_i}, 'StageOnly')
                    bout_protocol = [bout_protocol; 4*ones(length(bouts_t), 1)];
                elseif strcmp(trial_types{type_i}, 'Coupled') && strcmp(replay_types{rep_i}, 'ReplayOnly')
                    bout_protocol = [bout_protocol; 5*ones(length(bouts_t), 1)];
                elseif strcmp(trial_types{type_i}, 'EncoderOnly') && strcmp(replay_types{rep_i}, 'ReplayOnly')
                    bout_protocol = [bout_protocol; 6*ones(length(bouts_t), 1)];
                end
                
            end
            
        end
        
    end
    
    
    % start times of the bouts
    start_t = [bouts(:).start_time];
    
    % use bout timings to get raster data
    [spike_rates, spike_times] = ...
        get_raster_data(clusters, start_t, common_t, fs, options.spiking);
    
    % get velocity traces around those start times
    velocity_traces = nan(length(prepad*fs : postpad*fs-1), length(bouts));
    
    for bout_i = 1 : length(bouts)
        
        idx = bouts(bout_i).start_idx + (prepad*fs : postpad*fs-1);
        velocity_traces(:, bout_i) = bouts(bout_i).trial.velocity(idx);
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLOT RASTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %for each cluster
    for clust_i = 20%1 : length(clusters)
        
        % create a raster display object
        r = RasterDisplayFigure(length(protocols));
        
        % get the largest spike rate across the protocols
        M = -inf;
        
        % for each protocol
        for prot_i = 1 : length(protocols)
            
            % find the bouts belonging to the protocol
            idx = bout_protocol == prot_i;
            
            % get the required velocity traces and spike rates
            v =  velocity_traces(:, idx);
            sr = spike_rates(:, idx, clust_i);
            
            % make sure they are the same size
            assert(isequal(size(v), size(sr)));
            
            % update the maximum spike rate seen
            M = max([M, max(nanmean(sr, 2))]);
            
            % which subsection to fill
            x_i = mod(prot_i-1, 2) + 1;
            y_i = ceil(prot_i/2);
            
            % fill the section with the required data
            r.fill_data(x_i, y_i, spike_times{clust_i}(idx), v, sr, common_t, title_str{prot_i})
            
            % add grey bars to the end of the rasters
            %r.h_raster{x_i, y_i}.add_bars(end_t(idx) - start_t(idx), common_t(length(common_t)))
        end
        
        % set the limits and synchronize the axes of all the sections
        r.x_lim([prepad, postpad]);
        r.y_lim([0, 1.1*M], 3);
        r.sync_sections();
        
        % give the page a title
        FigureTitle(r.h_fig, sprintf('%s, Cluster %i, %s', probe_fnames{probe_i}, clusters(clust_i).id, clusters(clust_i).region_str));
        
        figs.save_fig_to_join(); % true, 400
        
    end
    
    figs.join_figs(sprintf('%s_rasters.pdf', probe_fnames{probe_i}));
    figs.clear_figs();
    
end
