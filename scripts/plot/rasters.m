%% settings for mismatch experiments
% list the experimental groups to plot, the label of the trial groups and
% the kind of raster you want plotted ('motion', 'mismatch', 'solenoid')
%   limits = time limits in s around the trigger point
%   common_fs = frequency of sampling for rasters
%   figure_dir = cell array containing location relative to main figure
%                   directory
experiment_groups       = {'mismatch_jul21'};
trial_group_labels      = {'R', 'T', 'RVT_gain_up', 'RV_gain_up'};
raster_trigger          = {'motion', 'motion', 'mismatch', 'mismatch'};  % {'solenoid', 'solenoid', 'solenoid', 'solenoid'};  % 
limits                  = [-1, 1];
common_fs               = 60;
save_figs               = true;
overwrite               = true;
figure_dir              = {'rasters', 'mismatch_jul21', 'motion'};



%% settings for visual flow experiments
% plot.experiment_groups       = {'visual_flow'};
% plot.trial_group_labels      = {'RVT', 'RV', 'RVT_gain_up', 'RT_gain_up'};
% plot.raster_trigger          = {'motion', 'motion', 'mismatch', 'mismatch'};
% plot.limits                  = [-1, 1];
% plot.common_fs               = 60;
% plot.save_figs               = true;
% plot.figure_dir              = {'rasters', 'visual_flow'};



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
            [traces{kk}, t, times] = data.get_traces_around_motion_onset(trial_group_labels{kk}, limits, common_fs);
        elseif strcmp(raster_trigger{kk}, 'mismatch')
            [traces{kk}, t, times] = data.get_traces_around_mismatch_onset(trial_group_labels{kk}, limits, common_fs);
        elseif strcmp(raster_trigger{kk}, 'solenoid')
            [traces{kk}, t, times] = data.get_traces_around_solenoid_up(trial_group_labels{kk}, limits, common_fs);
        end
        
        for jj = 1 : length(clusters)
            
            rd{jj}{kk} = RasterData(clusters(jj));
            rd{jj}{kk}.trigger_times = times;
        end
    end
    
    
    % for each cluster
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
