experiment = 'passive';

config = RC2AnalysisConfig();

figs = RC2Figures(config);
figs.save_on = true;
figs.set_figure_subdir(experiment, 'rasters');

% get details of the experiment
probe_fnames = experiment_details(experiment);

% number of seconds before and after the bout onset to take
prepad              = -1;
postpad             = 1;

% default options
options = default_options();

% fundamental trial types
protocol_ids    = 1 : 3;
% trial_types     = PassiveExperiment.protocol_type;
% visual_types    = PassiveExperiment.protocol_vis_stim;
title_str       = {'VT', 'V', 'T'};




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    
    exp_obj = PassiveExperiment(data, config);
    
    % filter clusters according to the cluster ID list file
    clusters = data.selected_clusters();
    
    % sampling rate of the recordings
    fs = data.sample_rate;
    
    % traces will be on a common timebase
    common_t = linspace(prepad, postpad, round((postpad-prepad)*fs));
    
    bouts           = [];
    bout_protocol   = [];
    
    for prot_i = protocol_ids
        
        trials = exp_obj.trials_of_type(prot_i);
        
        for trial_i = 1 : length(trials)
            
            bouts_t = trials(trial_i).motion_bouts(true);
            bouts = [bouts, bouts_t];
            bout_protocol = [bout_protocol; protocol_ids(prot_i)*ones(length(bouts_t), 1)];
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
    for clust_i = 1 : length(clusters)
        
        % create a raster display object
        r = RasterDisplayFigure(length(protocol_ids));
        
        % get the largest spike rate across the protocols
        M = -inf;
        
        % for each protocol
        for prot_i = protocol_ids
            
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
        
        figs.save_fig_to_join(true, 400); % 
        
    end
    
    figs.join_figs(sprintf('%s_rasters.pdf', probe_fnames{probe_i}));
    figs.clear_figs();
    
end
