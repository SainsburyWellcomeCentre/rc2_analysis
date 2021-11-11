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
        
            c = c + 1;
            [~, p_svm(ii, jj, kk), direction(ii, jj, kk)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{kk});
            tuning{ii, jj, kk} = data.load_tuning_curves(clusters(jj).id, trial_group_labels{kk});
        end
        
        % extra info to plot on the figure
        probe_id{ii, jj} = probe_ids{ii};
        cluster_id(ii, jj) = clusters(jj).id;
        cluster_region{ii, jj} = clusters(jj).region_str;
    end
end
        

%% plot

for ii = 1 : length(probe_ids)  
    
    for jj = 1 : length(clusters)
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        
        tuning_curve_plot       = {};
        shuff_histogram         = {};
        
        for kk = 1 : length(trial_group_labels)
            
            pos         = plot_array.get_position(2*kk-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            tuning_curve_plot{kk} = TuningCurvePlot(h_ax);
            tuning_curve_plot{kk}.plot(tuning{ii, jj, kk}, p_svm(ii, jj, kk), direction(ii, jj, kk));
            title(gca, trial_group_labels{kk});
            if kk == length(trial_group_labels)
                tuning_curve_plot{kk}.xlabel('Speed (cm/s)');
                tuning_curve_plot{kk}.ylabel('Firing rate (Hz)');
            end
            
            pos         = plot_array.get_position(2*kk);
            h_ax2        = axes('units', 'centimeters', 'position', pos);
            
            shuff_histogram{kk} = TuningCurveHistogram(h_ax2);
            shuff_histogram{kk}.plot(tuning{ii, jj, kk});
        end
        
        mx              = cellfun(@(x)(x.xmin), tuning_curve_plot);
        Mx              = cellfun(@(x)(x.xmax), tuning_curve_plot);
        My              = cellfun(@(x)(x.ymax), tuning_curve_plot);
        
        for kk = 1 : length(tuning_curve_plot)
            tuning_curve_plot{kk}.xlim([mx, Mx]);
            tuning_curve_plot{kk}.ylim([0, My]);
        end
        
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                probe_id{ii, jj}, ...
                cluster_id(ii, jj), ...
                cluster_region{ii, jj}));
        
        figs.save_fig_to_join();
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end


        
        
        for kk = 1 : length(trial_group_labels)
            
            pos         = plot_array.get_position(2*kk-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            % get info about the tuning curve for this cluster
            [fr, sd, n, x, shuff, stat_fr, stat_sd, stat_n] = ...
                data.tuning_curve_info(clusters(jj).id, trial_group_labels{kk});
            
            [~, p_signrank] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{kk});
            
            tuning_curve(kk)   = TuningCurvePlot(fr, sd, n, x, shuff, stat_fr, stat_sd, stat_n, p_signrank, h_ax);
            
            pos         = plot_array.get_position(2*prot_i);
            h_ax2        = axes('units', 'centimeters', 'position', pos);
            
            shuff_histogram(prot_i)   = TuningCurveHistogram(shuff, h_ax2);
            
            if kk == length(trial_group_labels)
                shuff_histogram(kk).xlabel('Speed (cm/s)');
                shuff_histogram(kk).ylabel('Firing rate (Hz)');
            end
            
            shuff_histogram(kk).title(trial_group_labels{kk});
            
%             recording_name{end+1, 1}   = probe_ids{ii};
%             cluster_id(end+1, 1)       = clusters(jj).id;
%             protocol{end+1, 1}         = exp_obj.protocol_label{prot_i};
%             r(end+1, 1)                = shuff.r;
%             R_sq(end+1, 1)             = shuff.rsq;
%             fit_slope(end+1, 1)        = shuff.beta(1);
%             p_val_shuffled(end+1, 1)   = shuff.p;
%             p_val_linear_fit(end+1, 1) = shuff.p_lm;
%             is_significant_shuffled(end+1, 1)   = shuff.p < 0.05;
%             is_significant_linear_fit(end+1, 1)   = shuff.p_lm < 0.05;
        end
        
        mx               = min([tuning_curve(:).xmin]);
        Mx               = max([tuning_curve(:).xmax]);
        my               = min([tuning_curve(:).ymin]);
        My               = max([tuning_curve(:).ymax]);
        
        for kk = 1 : length(tuning_curve)
            tuning_curve(kk).xlim([mx, Mx]);
            tuning_curve(kk).ylim([my, My]);
        end
        
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                probe_ids{ii}, ...
                clusters(jj).id, ...
                clusters(jj).region_str));
        
        figs.save_fig_to_join();
        
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end
