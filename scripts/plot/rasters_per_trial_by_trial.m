close all;

% Modified version of rasters_per_trial.m that groups by trial instead of by cluster
%
%   This script creates raster plots where each page shows one trial,
%   with all clusters displayed together in the same raster plot.

%   The time axis is aligned to solenoid timing:
%       - t=0 is when solenoid goes DOWN
%       - Display range: from -5s (before solenoid down) to +10s after solenoid goes UP
%       - This means the displayed time range varies per trial
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

%%
experiment_groups       = {'ambient_light'};
trial_group_labels      = {'RT'};
common_fs               = 60;
save_figs               = true;
overwrite               = true;
figure_dir              = {'rasters_by_trial', 'ambient_light'};

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
        
        % Find solenoid timing for this trial
        % Solenoid goes down (high to low): diff(solenoid > 2.5) == -1
        % Solenoid goes up (low to high): diff(solenoid > 2.5) == 1
        sol_down_idx = find(diff(trial.solenoid > 2.5) == -1, 1);
        sol_up_idx = find(diff(trial.solenoid > 2.5) == 1, 1);
        
        if isempty(sol_down_idx) || isempty(sol_up_idx)
            warning('Trial %d: Could not find solenoid transitions, skipping...', trial.trial_id);
            continue;
        end
        
        % Get absolute times for solenoid events
        sol_down_time = trial.probe_t(sol_down_idx);
        sol_up_time = trial.probe_t(sol_up_idx);
        
        % Use the entire trial time range, but make it relative to solenoid down
        % So t=0 is at solenoid down
        trial_start_time = trial.probe_t(1);
        trial_end_time = trial.probe_t(end);
        
        % Set limits relative to solenoid down (t=0)
        limits = [trial_start_time - sol_down_time, trial_end_time - sol_down_time];
        
        % Create one raster display object for this trial (with 1 section)
        r = RasterDisplayFigure(1);
        
        % Combine all clusters' spike data for this trial
        all_spike_times = {};
        all_spike_rates = [];
        
        for cluster_idx = 1 : length(clusters)
            
            cluster = clusters(cluster_idx);
            
            % Create RasterData for this cluster and trial
            rd = RasterData(cluster, limits);
            rd.trigger_times = sol_down_time;  % Set t=0 at solenoid down
            
            % Get spike data for this cluster around solenoid down time
            spike_times = rd.spike_array();
            spike_rates = rd.spike_convolutions();
            
            % Add to combined data
            all_spike_times = [all_spike_times; spike_times];
            all_spike_rates = [all_spike_rates, spike_rates];
        end
        
        % Get velocity trace for this trial
        % Extract velocity from trial and interpolate to common time points
        trial_velocity = trial.velocity;
        trial_times = trial.probe_t;
        
        velocity_trace = nan(length(rd.common_t), 1);  % Single column vector
        
        % Interpolate entire trial velocity to common time points (relative to solenoid down)
        common_times_absolute = rd.common_t + sol_down_time;
        velocity_trace(:, 1) = interp1(trial_times, trial_velocity, common_times_absolute, 'linear', nan);
        
        % fill the section with the combined data
        sol_duration = sol_up_time - sol_down_time;
        title_str = sprintf('Trial %i - All Clusters (Solenoid down at t=0, up at t=%.2fs)', trial.trial_id, sol_duration);
        r.fill_data(1, 1, all_spike_times, velocity_trace, all_spike_rates, rd.common_t, title_str);
        
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