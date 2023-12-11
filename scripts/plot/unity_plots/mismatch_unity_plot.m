% Plot mismatch unity plots
%
%   Specify options:
%
%       experiment_groups:      Will generate unity plots including all
%                               clusters in the probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21'
%                               Should be a cell array of strings with each
%                               entry an experiment group
%
%       trial_group_labels:     Will generate a unity plot for each of the
%                               trials specified. 
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%                               e.g. {'RVT_gain_up', 'RV_gain_up'}
%                               will create two unity plots, the first
%                               using data from 'RVT_gain_up' trials and
%                               the second using data from 'RV_gain_up'
%                               trials. 
%
%       marker_type:            the type of marker to use on the plot (any
%                               of those accepted by matlab) e.g. 'o' or
%                               'v'. Should be a cell array of the same
%                               length as `trial_group_labels`.
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
% If `save_figs` is true, one pdf with a unity plot for each of the trial
% groups specified


%%
% experiment_groups       = {'darkness','mismatch_darkness_oct21'};
% trial_group_labels      = {{'T_bank', 'T_RT', 'T_R', 'T'}};
experiment_groups       = {'mismatch_nov20'};
trial_group_labels      = {'RVT_gain_down'};
marker_style            = {'o'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'mismatch_unity_plot'};


%%
ctl                     = RC2Analysis();
ctl.setup_figures(figure_dir, save_figs);

probe_ids               = ctl.get_probe_ids(experiment_groups{:});
mm                      = MismatchAnalysis();

x_all                   = cell(1, length(trial_group_labels));
y_all                   = cell(1, length(trial_group_labels));
p_val                   = cell(1, length(trial_group_labels));
direction               = cell(1, length(trial_group_labels));

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.VISp_clusters;
    
    for jj = 1 : length(trial_group_labels)
        
        trials      = data.get_trials_with_trial_group_label(trial_group_labels{jj});
        
        for kk = 1 : length(clusters)
            
            x_all{jj}(end+1) = mm.get_avg_baseline_fr(clusters(kk), trials);
            y_all{jj}(end+1) = mm.get_avg_response_fr(clusters(kk), trials);
            [~, p_val{jj}(end+1), direction{jj}(end+1)] = mm.is_response_significant(clusters(kk), trials);
        end
    end
end



%%
h_fig                   = a4figure();
plot_array              = PlotArray(3, 2);
u                       = UnityPlotPopulation.empty();


for ii = 1 : length(trial_group_labels)
    
    pos         = plot_array.get_position(ii);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    u(end+1)   = UnityPlotPopulation(x_all{ii}, ...
        y_all{ii}, ...
        p_val{ii}, ...
        direction{ii}, ...
        h_ax);
    
    u(end).marker_style = marker_style{ii};
    
    u(end).plot();
    
    u(end).xlabel('Baseline FR (Hz)');
    u(end).ylabel('Response FR (Hz)');
    
    u(end).title(trial_group_labels{ii});
    u(end).add_histogram(1);
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for ii = 1 : length(u)
    u(ii).xlim([m, M]);
end


ctl.figs.save_fig('unity_plot', overwrite);
