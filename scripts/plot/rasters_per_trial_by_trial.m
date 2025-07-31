close all;

% Modified version of rasters_per_trial.m that groups by trial instead of by cluster
%
%   This script creates raster plots where each page shows one trial,
%   with all clusters displayed together in the same raster plot.
%
%   Specify options:
%
%       experiment_groups:      Will generate raster plots for all trials
%                               and all probe recordings 
%                               in the specified experiment group.
%
%       trial_group_labels:     Will generate a raster for all trials specified.
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%
%       raster_trigger:         For each of the entries in
%                               'trial_group_labels', a string specifying which
%                               aspect of the trials to trigger on. Should
%                               be either 'motion' (onset of motion) or
%                               'mismatch' (onset of a mismatch event).
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
%                               directory to save pdf's.
%
%       min_bout_duration:      the minimum duration in s for a motion bout to be included
%
%       include_200ms:          whether or not to include the first 200ms
%                               after the solenoid goes low to look for
%                               motion onset  

%%
experiment_groups       = {'training_running'};
trial_group_labels      = {'RT'};
raster_trigger          = {'motion'};
limits                  = [-5, 10];
common_fs               = 60;
save_figs               = true;
overwrite               = true;
figure_dir              = {'rasters_by_trial', 'training_running'};
min_bout_duration       = 2;
include_200ms           = true;

%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.selected_clusters();
    
    % Get all trials for the specified trial group
    all_trials = {};
    for kk = 1 : length(trial_group_labels)
        trials = data.get_trials_with_trial_group_label(trial_group_labels{kk});
        all_trials = [all_trials; trials];
    end
    
    % For each trial
    for trial_idx = 1 : length(all_trials)
        
        trial = all_trials{trial_idx};
        
        % Get motion onset times for this trial
        if strcmp(raster_trigger{1}, 'motion')
            options.min_bout_duration = min_bout_duration;
            options.include_200ms = include_200ms;
            
            % Get motion bouts for this specific trial
            bouts = trial.motion_bouts(options.include_200ms, true);
            motion_times = cellfun(@(x)(x.start_time), bouts);
            
            % Filter bouts by minimum duration
            bout_durations = cellfun(@(x)(x.duration), bouts);
            valid_bouts = bout_durations >= options.min_bout_duration;
            motion_times = motion_times(valid_bouts);
            
        elseif strcmp(raster_trigger{1}, 'mismatch')
            motion_times = trial.mismatch_onset_t;
        end
        
        if isempty(motion_times)
            continue; % Skip trials with no events
        end
        
        % Create one raster display object for this trial (with 1 section)
        r = RasterDisplayFigure(1);
        
        % Combine all clusters' spike data for this trial
        all_spike_times = {};
        all_spike_rates = [];
        
        for cluster_idx = 1 : length(clusters)
            
            cluster = clusters(cluster_idx);
            
            % Create RasterData for this cluster and trial
            rd = RasterData(cluster, limits);
            rd.trigger_times = motion_times;
            
            % Get spike data for this cluster around motion times
            spike_times = rd.spike_array();
            spike_rates = rd.spike_convolutions();
            
            % Add to combined data
            all_spike_times = [all_spike_times; spike_times];
            all_spike_rates = [all_spike_rates, spike_rates];
        end
        
        % Get velocity traces for this trial around motion times
        % This would need to be implemented based on your data structure
        % For now, we'll use placeholder data
        velocity_traces = nan(length(rd.common_t), length(motion_times));
        
        % fill the section with the combined data
        r.fill_data(1, 1, all_spike_times, velocity_traces, all_spike_rates, rd.common_t, sprintf('Trial %i - All Clusters', trial.trial_id));
        
        % set the limits and synchronize the axes of all the sections
        r.x_lim(limits);
        r.sync_sections();
        
        % give the page a title
        FigureTitle(r.h_fig, sprintf('%s, Trial %i (%s)', probe_ids{ii}, trial.trial_id, trial.trial_group_label));
        
        ctl.figs.save_fig_to_join(true, 400);
    end
    
    ctl.figs.join_figs(sprintf('%s.pdf', probe_ids{ii}), overwrite);
    ctl.figs.clear_figs();
end 