% Plot heatmaps and population FR average, triggered on some feature of a trial
%
%   Specify options:
%
%       experiment_groups:      Will generate plots combining all clusters
%                               for all probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21'
%                               Should be a cell array of strings with each
%                               entry an experiment group
%
%       trial_group_labels:     Will generate a heatmap by averaging firing
%                               rates across all trials specified.
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%                               e.g. {'R', {'T_bank', 'T_RT', 'T_R'}, 'RT'}
%                               will create three heatmaps, the first
%                               averaging across 'R' (running) trials, the
%                               second averaging across any trial of
%                               'T_bank', 'T_RT' or 'T_R' type, and the
%                               third averaging across 'RT'
%                               (running+translation) trials.
%
%       heatmap_trigger:        For each of the entries in
%                               'trial_group_labels', a string specifying which
%                               aspect of the trials to trigger on. Should
%                               be either 'motion' (onset of motion) or
%                               'mismatch' (onset of a mismatch event).
%
%       limits:                 time in seconds around the event to display
%                               for the heatmap. e.g. [-1, 1] will display
%                               the heatmap from 1 second before to 1
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
%   If any of the elemenets in `heatmap_trigger` are 'motion', the
%   following options should also be specified:
%
%       min_bout_duration:      the minimum duration in s for a motion bout to be included
%
%       include_200ms:          whether or not to include the first 200ms
%                               after the solenoid goes low to look for
%                               motion onset  
%
%
% If `save_figs` is true, one pdf will be created for the heatmpas and one
% for the population FR average.


%%
% experiment_groups       = {'darkness'};
% trial_group_labels      = {'R', {'T_bank', 'T_RT', 'T_R'}, 'RT'};
% heatmap_trigger         = {'motion', 'motion', 'motion'};
% limits                  = [-1, 1];
% common_fs               = 60;
% save_figs               = false;
% overwrite               = true;
% figure_dir              = {'heatmaps', 'darkness'};
% min_bout_duration       = 2;
% include_200ms           = true;


%%
% experiment_groups       = {'visual_flow'};
% trial_group_labels      = {'RVT', 'RV', 'VT_RVT', 'VT_RV', 'V_RVT', 'V_RV'};
% heatmap_trigger         = {'motion', 'motion', 'motion', 'motion', 'motion', 'motion'};
% limits                  = [-1, 1];
% common_fs               = 60;
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'heatmaps', 'visual_flow'};
% min_bout_duration       = 2;
% include_200ms           = true;


%%
% experiment_groups       = {'mismatch_jul21'};
% trial_group_labels      = {'R', 'T', 'RVT_gain_up', 'RV_gain_up'};
% heatmap_trigger         = {'motion', 'motion', 'mismatch', 'mismatch'};  % {'solenoid', 'solenoid', 'solenoid', 'solenoid'};  % 
% limits                  = [-1, 1];
% common_fs               = 60;
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'heatmaps', 'mismatch_jul21'};
% min_bout_duration       = 2;
% include_200ms           = true;


%%
% experiment_groups       = {'mismatch_darkness_oct21'};
% trial_group_labels      = {'R', 'T', 'RT_gain_up'};
% heatmap_trigger         = {'motion', 'motion', 'mismatch'};
% limits                  = [-1, 1];
% common_fs               = 60;
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'heatmaps', 'mismatch_darkness_oct21'};
% min_bout_duration       = 2;
% include_200ms           = true;


%%
experiment_groups       = {'visual_flow', 'mismatch_nov20', 'mismatch_jul21'};
trial_group_labels      = {{'RVT', 'RVT_gain_up'}, {'RV', 'RV_gain_up'}};
heatmap_trigger         = {'motion', 'motion'};
limits                  = [-1, 1];
common_fs               = 60;
save_figs               = true;
overwrite               = true;
figure_dir              = {'heatmaps', 'visual_flow+mismatch_nov20+mismatch_jul21'};
min_bout_duration       = 2;
include_200ms           = true;


%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
rd                      = cell(1, length(trial_group_labels));
mm_response             = cell(1, length(trial_group_labels));


ctl.setup_figures(figure_dir, save_figs);


for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.VISp_clusters();
    
    for kk = 1 : length(trial_group_labels)
        
        if strcmp(heatmap_trigger{kk}, 'motion')
            options.min_bout_duration = min_bout_duration;
            options.include_200ms = include_200ms;
            [~, t, times] = data.get_traces_around_motion_onset(trial_group_labels{kk}, limits, common_fs, options);
        elseif strcmp(heatmap_trigger{kk}, 'mismatch')
            [~, t, times] = data.get_traces_around_mismatch_onset(trial_group_labels{kk}, limits, common_fs);
        elseif strcmp(heatmap_trigger{kk}, 'solenoid')
            [~, t, times] = data.get_traces_around_solenoid_up(trial_group_labels{kk}, limits, common_fs);
        end
        
        for jj = 1 : length(clusters)
            
            rd{kk}{end+1} = RasterData(clusters(jj));
            rd{kk}{end}.trigger_times = times;
            rd{kk}{end}.fs = 10e3;
            
            if strcmp(heatmap_trigger{kk}, 'mismatch')
                mm_response{kk}(end+1) = data.get_mismatch_response(clusters(jj).id, trial_group_labels{kk});
            end 
        end
    end
end



%% heatmaps
ctl.setup_figures(figure_dir, save_figs);
ctl.figs.a4figure();
plot_array              = PlotArray(3, 2);
heatmap                 = cell(1, length(trial_group_labels));

for ii = 1 : length(trial_group_labels)
    
    hd = HeatmapData(rd{ii});
    
    if strcmp(heatmap_trigger{ii}, 'mismatch')
        [~, idx] = sort(mm_response{ii});
    else
        idx = hd.heatmap_response_order();
    end
    
    heatmap{ii} = hd.sort_heatmap(idx);
    
    pos         = plot_array.get_position(ii);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    hp = HeatmapPlot(h_ax);
    hp.plot(hd.common_t, heatmap{ii});
    
    hp.title(trial_group_labels{ii});
end

ctl.figs.save_fig('heatmaps', overwrite);
    



%% population average trace
h_fig                   = ctl.figs.a4figure();
plot_array              = PlotArray(3, 2);
pp                      = {};

for ii = 1 : length(trial_group_labels)
    
    population_average = mean(heatmap{ii}, 1);
    population_sem = std(heatmap{ii}, [], 1) / sqrt(size(heatmap{ii}, 1));
    
    pos         = plot_array.get_position(ii);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    pp{ii} = PopulationFRPlot(h_ax);
    pp{ii}.plot(hd.common_t, population_average, population_sem);
    pp{ii}.title(trial_group_labels{ii});
end

%sync the yaxes
m = min(cellfun(@(x)(x.ymin), pp));
M = max(cellfun(@(x)(x.ymax), pp));

for ii = 1 : length(trial_group_labels)
    pp{ii}.ylim([m, M]);
end

ctl.figs.save_fig('population_average_fr', overwrite);
