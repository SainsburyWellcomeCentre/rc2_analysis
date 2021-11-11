% plot speed tuning curve for an experimental group
%  tuning cuves must have been generated and saved using
%  RC2Analysis.create_tuning_curves

% experiment_groups       = {'mismatch_darkness_oct21'};
% trial_group_labels      = {'R', 'T', 'RT_gain_up'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'tuning', 'mismatch_darkness_oct21'};


%%
experiment_groups       = {'darkness'};
trial_group_labels      = {'RT', 'R', {'T_bank', 'T_RT', 'T_R'}};
save_figs               = true;
overwrite               = false;
figure_dir              = {'tuning_curves', 'darkness'};



%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

tuning          = {};
p_svm           = [];
direction       = [];
probe_id        = {};
cluster_id      = [];
cluster_region  = {};

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.VISp_clusters();
    
    for jj = 1 : length(clusters)
        
        for kk = 1 : length(trial_group_labels)
            
            [~, p_svm(ii, jj, kk), direction(ii, jj, kk)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{kk});
            tuning{ii}{jj}{kk} = data.load_tuning_curves(clusters(jj).id, trial_group_labels{kk});
        end
        
        % extra info to plot on the figure
        probe_id{ii, jj} = probe_ids{ii};
        cluster_id(ii, jj) = clusters(jj).id;
        cluster_region{ii, jj} = clusters(jj).region_str;
    end
end
        

%% plot

for ii = 1 : length(probe_ids)  
    
    for jj = 1 : length(tuning{ii})
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        
        tuning_curve_plot       = {};
        shuff_histogram         = {};
        
        for kk = 1 : length(trial_group_labels)
            
            pos         = plot_array.get_position(2*kk-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            tuning_curve_plot{kk} = TuningCurvePlot(h_ax);
            tuning_curve_plot{kk}.plot(tuning{ii}{jj}{kk}, p_svm(ii, jj, kk), direction(ii, jj, kk));
            title(gca, trial_group_labels{kk});
            if kk == length(trial_group_labels)
                tuning_curve_plot{kk}.xlabel('Speed (cm/s)');
                tuning_curve_plot{kk}.ylabel('Firing rate (Hz)');
            end
            
            pos         = plot_array.get_position(2*kk);
            h_ax2        = axes('units', 'centimeters', 'position', pos);
            
            shuff_histogram{kk} = TuningCurveHistogram(h_ax2);
            shuff_histogram{kk}.plot(tuning{ii}{jj}{kk});
        end
        
        mx              = min(cellfun(@(x)(x.xmin), tuning_curve_plot));
        Mx              = max(cellfun(@(x)(x.xmax), tuning_curve_plot));
        My              = max(cellfun(@(x)(x.ymax), tuning_curve_plot));
        
        for kk = 1 : length(tuning_curve_plot)
            tuning_curve_plot{kk}.xlim([mx, Mx]);
            tuning_curve_plot{kk}.ylim([0, My]);
        end
        
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                probe_id{ii, jj}, ...
                cluster_id(ii, jj), ...
                cluster_region{ii, jj}));
        
        ctl.figs.save_fig_to_join();
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end

