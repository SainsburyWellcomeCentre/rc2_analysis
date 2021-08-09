% plot rasters for an animal
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experiment              = 'mismatch_nov20';  % 'darkness', 'visual_flow', 'head_tilt', 'mismatch_nov20'

config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, 'rasters');

% get details of the experiment
probe_fnames            = experiment_details(experiment);

% number of seconds before and after the bout onset to take
prepad                  = -1;
postpad                 = 1;


if strcmp(experiment, 'visual_flow')
    protocol_ids    = 1 : 6;
elseif strcmp(experiment, 'mismatch_nov20')
    protocol_ids    = [2, 4];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.selected_clusters();
    
    % sampling rate of the recordings
    fs                  = data.sample_rate;
    
    % traces will be on a common timebase
    common_t = linspace(prepad, postpad, round((postpad-prepad)*fs));
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(data, config);
    elseif strcmp(experiment, 'mismatch_nov20')
        exp_obj         = MismatchExperiment(data, config);
    end
    
    % store spike times, spike rates and velocity traces for each protocol
    spike_rates = {};
    spike_times = {};
    velocity_traces = {};
    
    for prot_i = 1 : length(protocol_ids)
        
        all_bouts = exp_obj.motion_bouts_by_protocol(protocol_ids(prot_i), true, true, 1);
        
        % move bouts below 2s
        all_bouts = all_bouts([all_bouts(:).duration] > 2);
        
        % start time of the bouts
        start_t = [all_bouts(:).start_time];
        start_idx = [all_bouts(:).start_idx];
        
        % use bout timings to get raster data
        for clust_i = 1 : length(clusters)
            fr = FiringRate(clusters(clust_i).spike_times);
            for bout_i = 1 : length(all_bouts)
                spike_rates{prot_i}{clust_i}(:, bout_i) = fr.get_convolution(start_t(bout_i) + common_t);
                spike_times{prot_i}{clust_i}{bout_i} = ...
                    fr.restrict_times(start_t(bout_i) + common_t([1, end])) - start_t(bout_i);
            end
        end
        
        velocity_traces{prot_i} = [];
        
        for bout_i = 1 : length(start_idx)
            idx = start_idx(bout_i) + (prepad*fs : postpad*fs-1);
            velocity_traces{prot_i}(:, bout_i) = all_bouts(bout_i).trial.velocity(idx);
        end
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
        for prot_i = 1 : length(protocol_ids)
            
            % get the required velocity traces and spike rates
            tm = spike_times{prot_i}{clust_i};
            v =  velocity_traces{prot_i};
            sr = spike_rates{prot_i}{clust_i};
            
            % make sure they are the same size
            assert(isequal(size(v), size(sr)));
            
            % update the maximum spike rate seen
%             M = max([M, max(nanmean(sr, 2))]);
            
            % which subsection to fill
            x_i = mod(prot_i-1, 2) + 1;
            y_i = ceil(prot_i/2);
            
            % fill the section with the required data
            r.fill_data(x_i, y_i, tm, v, sr, common_t, exp_obj.get_protocol_label(protocol_ids(prot_i)))
            
            
            M = max([M, max(r.h_rates{x_i, y_i}.h_ax.YLim)]);
            % add grey bars to the end of the rasters
            %r.h_raster{x_i, y_i}.add_bars(end_t(idx) - start_t(idx), common_t(length(common_t)))
        end
        
        % set the limits and synchronize the axes of all the sections
        r.x_lim([prepad, postpad]);
        r.y_lim([0, 1.1*M], 3);
        r.sync_sections();
        
        % give the page a title
        c = Cluster(clusters(clust_i));
        
        FigureTitle(r.h_fig, sprintf('%s, Cluster %i, %s, %s', probe_fnames{probe_i}, clusters(clust_i).id, clusters(clust_i).region_str, c.spiking_class));
        
        figs.save_fig_to_join(true, 400); %
    end
    
    figs.join_figs(sprintf('%s_rasters.pdf', probe_fnames{probe_i}));
    figs.clear_figs();
end
