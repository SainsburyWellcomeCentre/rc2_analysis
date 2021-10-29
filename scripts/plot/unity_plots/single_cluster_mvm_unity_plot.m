%%script for producing a unity plot for stationary vs. motion for each
%%cluster and trial group

experiment_groups       = {'visual_flow'};
trial_group_labels      = {'RVT', 'RV', 'VT_RVT', 'VT_RV', 'V_RVT', 'V_RV'};
marker_style            = {'o', 'o', 'o', 'o', 'o', 'o'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'unity_plots', 'visual_flow', 'motion_vs_motion', 'single_cluster'};



%%
% experiment_groups       = {'darkness'};
% trial_group_labels      = {'RT', 'R', 'T_bank', 'T_RT', 'T_R'};
% marker_style            = {'o', 'o', 'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'darkness', 'motion_vs_motion', 'single_cluster'};



%%
% experiment_groups       = {'mismatch_jul21'};
% trial_group_labels      = {'RVT_gain_up', 'RV_gain_up', 'R', 'T'};
% marker_style            = {'o', 'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'mismatch_jul21', 'motion_vs_motion', 'single_cluster'};



%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
mm                      = MismatchAnalysis();

x_all                   = cell(length(probe_ids));
y_all                   = cell(length(probe_ids));
p_val                   = cell(length(probe_ids));
direction               = cell(length(probe_ids));
cluster_ids             = cell(1, length(probe_ids));
region_str              = cell(1, length(probe_ids));


for ii = 1 : length(probe_ids)
    
    data                = ctl.load_formatted_data(probe_ids{ii});
    clusters            = data.VISp_clusters;
    
    % get list of cluster IDs and the re
    cluster_ids{ii}     = data.VISp_cluster_ids;
    region_str{ii}      = {clusters(:).region_str};
    
    for jj = 1 : length(trial_group_labels)-1
        for kk = jj+1 : length(trial_group_labels)
            
            % skip if the trial group label is not in the experiment
            if ~data.check_trial_group(trial_group_labels{jj})
                continue
            end
            
            for zz = 1 : length(cluster_ids{ii})
                
                x_fr = data.motion_fr_for_trial_group(cluster_ids{ii}(zz), trial_group_labels{kk});
                y_fr = data.motion_fr_for_trial_group(cluster_ids{ii}(zz), trial_group_labels{jj});
                
                M = max(length(x_fr), length(y_fr));
                
                x_all{ii}{zz}{jj, kk} = x_fr(1:M);
                y_all{ii}{zz}{jj, kk} = y_fr(1:M);
                
                [~, ~, p_val{ii}{zz}(jj, kk), direction{ii}{zz}(jj, kk)] = ...
                    compare_groups_with_signrank(x_all{ii}{zz}{jj, kk}, y_all{ii}{zz}{jj, kk});
            end
        end
    end
end




%%
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    for jj = 1 : length(cluster_ids)
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(6, 6);
        u                       = UnityPlotSingleCluster.empty();

        for kk = 1 : length(trial_group_labels)-1
            for zz = kk+1 : length(trial_group_labels)

                sp_idx      = (kk-1)*length(trial_group_labels) + zz;
                pos         = plot_array.get_position(sp_idx);
                h_ax        = axes('units', 'centimeters', 'position', pos);

                u(end+1)   = UnityPlotSingleCluster(x_all{ii}{jj}{kk, zz}, ...
                                                y_all{ii}{jj}{kk, zz}, ...
                                                p_val{ii}{jj}(kk, zz), ...
                                                direction{ii}{jj}(kk, zz), ...
                                                h_ax);

                u(end).marker_style = 'o';

                u(end).plot();

                u(end).xlabel(sprintf('%s (Hz)', trial_group_labels{zz}));
                u(end).ylabel(sprintf('%s (Hz)', trial_group_labels{kk}));
            end
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


