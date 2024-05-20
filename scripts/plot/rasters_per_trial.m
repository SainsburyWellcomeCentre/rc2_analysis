close all;

% modified from rasters.m - this code is not finished and might never be
%
%   Specify options:
%
%       experiment_groups:      Will generate raster plots for all clusters
%                               and all probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21'
%                               Should be a cell array of strings with each
%                               entry an experiment group
%
%       trial_group_labels:     Will generate a raster for all trials specified.
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%                               e.g. {'R', {'T_bank', 'T_RT', 'T_R'}, 'RT'}
%                               will create three heatmaps, the first
%                               for all 'R' (running) trials, the
%                               second for all trials of any of 
%                               'T_bank', 'T_RT' or 'T_R' type, and the
%                               third for all 'RT'
%                               (running+translation) trials.
%
%       raster_trigger:         For each of the entries in
%                               'trial_group_labels', a string specifying which
%                               aspect of the trials to trigger on. Should
%                               be either 'motion' (onset of motion) or
%                               'mismatch' (onset of a mismatch event).
%                               Should be a cell array of the same  
%                               length as `trial_group_labels`.
%
%       limits:                 time in seconds around the event to display
%                               for the raster. e.g. [-1, 1] will display
%                               the rasters from 1 second before to 1
%                               second after the event.
%
%       common_fs:              sampling frequency to compute the
%                               firing rate convolutions for each cluster
%                               (e.g. 60Hz) 
%
%       save_figs:              true or false, whether to save the figures to pdf
%
%       overwrite:              true or false. If figure pdf's already exist,
%                               whether to overwrite 
%       
%       figure_dir:             cell array of strings specifying which
%                               directory to save pdf's. The directory will
%                               be relative to the directory specified by
%                               path_config.figure_dir (in
%                               `path_config.m`), so that {'one', 'two',
%                               'three'} will save .pdfs to:
%                               <path_config.figure_dir>\one\two\three\
%
%   If any of the elemenets in `raster_trigger` are 'motion', the
%   following options should also be specified:
%
%       min_bout_duration:      the minimum duration in s for a motion bout to be included
%
%       include_200ms:          whether or not to include the first 200ms
%                               after the solenoid goes low to look for
%                               motion onset  
%
%
% If `save_figs` is true, one pdf will be created for each probe recording,
% and contain a A4 page for each cluster, containing the raster data for
% several conditions.


%%
experiment_groups       = {'passive_same_luminance'};
trial_group_labels      = {'T_Vstatic', 'V', 'VT'};
raster_trigger          = {'motion', 'motion', 'motion'};
limits                  = [-5, 10];
common_fs               = 60;
save_figs               = true;
overwrite               = true;
figure_dir              = {'rasters', 'passive_same_luminance'};
min_bout_duration       = 2.95;
include_200ms           = true;


%%
% experiment_groups       = {'mismatch_jul21'};
% trial_group_labels      = {'R', 'T', 'RVT_gain_up', 'RV_gain_up'};
% raster_trigger          = {'motion', 'motion', 'mismatch', 'mismatch'};
% limits                  = [-1, 1];
% common_fs               = 60;
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'rasters', 'mismatch_jul21', 'motion'};
% min_bout_duration       = 2;
% include_200ms           = true;


%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
rd                      = cell(size(trial_group_labels));

ctl.setup_figures(figure_dir, save_figs);



for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.selected_clusters();
    
    traces      = {};
    
     for kk = 1 : length(trial_group_labels)
        
        if strcmp(raster_trigger{kk}, 'motion')
            options.min_bout_duration = min_bout_duration;
            options.include_200ms = include_200ms;
            [traces{kk}, t, times] = data.get_traces_around_motion_onset(trial_group_labels{kk}, limits, common_fs, options);
        elseif strcmp(raster_trigger{kk}, 'mismatch')
            [traces{kk}, t, times] = data.get_traces_around_mismatch_onset(trial_group_labels{kk}, limits, common_fs);
        elseif strcmp(raster_trigger{kk}, 'solenoid')
            [traces{kk}, t, times] = data.get_traces_around_solenoid_up(trial_group_labels{kk}, limits, common_fs);
        end
        
        for jj = 1 : length(clusters)
            
            rd{jj}{kk} = RasterData(clusters(jj), limits);
            rd{jj}{kk}.trigger_times = times;
        end
     end
    
    
    % for each cluster
    jj = 1
     for jj = 1 : length(clusters)
        %
        % create a raster display object
        r = RasterDisplayFigure(length(trial_group_labels));
        
        % for each protocol
         
         for kk = 1 : length(trial_group_labels)
            
            % get the required velocity traces and spike rates
            v =  traces{kk};
            sr = rd{jj}{kk}.spike_convolutions();
            
            % make sure they are the same size
            assert(isequal(size(v), size(sr)));
            
            % which subsection to fill
            x_i = mod(kk-1, 2) + 1;
            y_i = ceil(kk/2);
            
            % fill the section with the required data
            r.fill_data(x_i, y_i, rd{jj}{kk}.spike_array(), v, sr, t, trial_group_labels{kk})
         end
        
        % set the limits and synchronize the axes of all the sections
        r.x_lim(limits);
%         r.y_lim([0, 1.1*r.max_spike_rate], 3);
        r.sync_sections();
        
        % give the page a title
        FigureTitle(r.h_fig, sprintf('%s, Cluster %i, %s', probe_ids{ii}, clusters(jj).id, clusters(jj).region_str));
        
        ctl.figs.save_fig_to_join(true, 400); % 
        
     end
    
    ctl.figs.join_figs(sprintf('%s.pdf', probe_ids{ii}), overwrite);
    ctl.figs.clear_figs();
end
