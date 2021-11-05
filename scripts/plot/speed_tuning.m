% plot speed tuning curve for an experimental group


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
figure_dir              = {'tuning', 'darkness'};



%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.VISp_clusters();
    
    for jj = 1 : length(clusters)
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        
        tuning_curve            = TuningCurvePlot.empty();
        shuff_histogram         = TuningCurveHistogram.empty();
        
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
