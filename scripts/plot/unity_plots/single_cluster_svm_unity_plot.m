% Plot stationary vs. motion unity plot for an experiment group and trial
% group for each cluster individually. i.e. show the individual trials on
% the unity plot (c.f. population_svm_unity_and_mi_plot).
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
%                               entry an experiment group.
%
%       trial_group_labels:     Will generate a unity plot for each of the
%                               trial groups specified. 
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%                               e.g. {'R', {'T_bank', 'T_RT', 'T_R'}, 'RT'}
%                               will create three unity plots, the first
%                               combining across 'R' (running) trials, the
%                               second combining across any trial of
%                               'T_bank', 'T_RT' or 'T_R' type, and the
%                               third combining across 'RT'
%                               (running+translation) trials.
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
%
% If `save_figs` is true, one pdf is saved for each probe recording with a
% page for each cluster showing the unity plots for each of the trial
% groups specified with `trial_group_labels`.


%%
% experiment_groups       = {'visual_flow'};
% trial_group_labels      = {'RVT', 'RV', 'VT_RVT', 'VT_RV', 'V_RVT', 'V_RV'};
% marker_style            = {'o', 'o', 'o', 'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'visual_flow', 'stationary_vs_motion', 'single_cluster'};


%%
% experiment_groups       = {'darkness'};
% trial_group_labels      = {'RT', 'R', 'T_bank', 'T_RT', 'T_R'};
% marker_style            = {'o', 'o', 'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'darkness', 'stationary_vs_motion', 'single_cluster'};


%%
experiment_groups       = {'mismatch_jul21'};
trial_group_labels      = {'RVT_gain_up', 'RV_gain_up', 'R', 'T'};
marker_style            = {'o', 'o', 'o', 'o'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'unity_plots', 'mismatch_jul21', 'stationary_vs_motion', 'single_cluster'};



%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
mm                      = MismatchAnalysis();

x_all                   = cell(1, length(probe_ids));
y_all                   = cell(1, length(probe_ids));
p_val                   = cell(1, length(probe_ids));
direction               = cell(1, length(probe_ids));
cluster_ids             = cell(1, length(probe_ids));
region_str              = cell(1, length(probe_ids));


for ii = 1 : length(probe_ids)
    
    data                = ctl.load_formatted_data(probe_ids{ii});
    clusters            = data.VISp_clusters;
    
    % get list of cluster IDs and the re
    cluster_ids{ii}     = data.VISp_cluster_ids;
    region_str{ii}      = {clusters(:).region_str};
    
    for jj = 1 : length(trial_group_labels)
        
        % skip if the trial group label is not in the experiment
        if ~data.check_trial_group(trial_group_labels{jj})
            continue
        end
        
        for kk = 1 : length(cluster_ids{ii})
            
            x_all{ii}{kk}{jj} = data.stationary_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            y_all{ii}{kk}{jj} = data.motion_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            [~, p_val{ii}{kk}(jj), direction{ii}{kk}(jj)] = ...
                data.is_stationary_vs_motion_significant(cluster_ids{ii}(kk), trial_group_labels{jj});
        end
    end
end




%%
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    for jj = 1 : length(cluster_ids)
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        u                       = UnityPlotSingleCluster.empty();

        for kk = 1 : length(trial_group_labels)
            
            pos         = plot_array.get_position(kk);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            u(end+1)   = UnityPlotSingleCluster(x_all{ii}{jj}{kk}, ...
                                            y_all{ii}{jj}{kk}, ...
                                            p_val{ii}{jj}(kk), ...
                                            direction{ii}{jj}(kk), ...
                                            h_ax);
    
            u(end).marker_style = marker_style{kk};

            u(end).plot();

            u(end).xlabel('Baseline FR (Hz)');
            u(end).ylabel('Response FR (Hz)');

            u(end).title(trial_group_labels{kk});
%             u(end).add_histogram(1);
        end
        
        % sync across all axes
        m               = min([u(:).min]);
        M               = max([u(:).max]);
        
        for kk = 1 : length(u)
            u(kk).xlim([m, M]);
        end
        
         % give the page a title
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', probe_ids{ii}, cluster_ids{ii}(jj), region_str{ii}{jj}));
        
        ctl.figs.save_fig_to_join(); %
    end
    
    ctl.figs.join_figs(sprintf('%s.pdf', probe_ids{ii}), overwrite);
    ctl.figs.clear_figs();
end


