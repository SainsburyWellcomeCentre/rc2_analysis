experiment              = 'darkness';

config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, 'tuning', 'all_conditions', datestr(now, 'yyyymmdd'));

probe_fnames            = experiment_details(experiment, 'protocols');

plot_array              = PlotArray(3, 2);

recording_name = {};
cluster_id = [];
protocol = {};
r = [];
R_sq = [];
fit_slope = [];
p_val_shuffled = [];
p_val_linear_fit = [];
is_significant_shuffled = [];
is_significant_linear_fit = [];

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.selected_clusters;
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(data, config);
    end
    
    for cluster_i = 1 : length(clusters)
        %%
        h_fig           = figs.a4figure();
        
        u               = TuningCurvePlot.empty();
        u2              = TuningCurveHistogram.empty();
        
        if clusters(cluster_i).id == 257
            disp('');
        end
        
        for prot_i = 1 : length(exp_obj.protocol_ids)
            
            pos         = plot_array.get_position(2*prot_i-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            [fr, sd, n, x, shuff, stat_fr, stat_sd, stat_n] = ...
                exp_obj.tuning_curve(clusters(cluster_i).id, prot_i);
            
            stat_rate = exp_obj.trial_stationary_fr(clusters(cluster_i).id, prot_i);
            mot_rate = exp_obj.trial_motion_fr(clusters(cluster_i).id, prot_i);
            
            p_signrank = signrank(stat_rate, mot_rate);
            
            u(prot_i)   = TuningCurvePlot(fr, sd, n, x, shuff, stat_fr, stat_sd, stat_n, p_signrank, h_ax);
            
            pos         = plot_array.get_position(2*prot_i);
            h_ax2        = axes('units', 'centimeters', 'position', pos);
            
            u2(prot_i)   = TuningCurveHistogram(shuff, h_ax2);
            
            if prot_i == length(exp_obj.protocol_ids)
                u(prot_i).xlabel('Speed (cm/s)');
                u(prot_i).ylabel('Firing rate (Hz)');
            end
            
            u(prot_i).title(exp_obj.protocol_label{prot_i});
            
            recording_name{end+1, 1}   = probe_fnames{probe_i};
            cluster_id(end+1, 1)       = clusters(cluster_i).id;
            protocol{end+1, 1}         = exp_obj.protocol_label{prot_i};
            r(end+1, 1)                = shuff.r;
            R_sq(end+1, 1)             = shuff.rsq;
            fit_slope(end+1, 1)        = shuff.beta(1);
            p_val_shuffled(end+1, 1)   = shuff.p;
            p_val_linear_fit(end+1, 1) = shuff.p_lm;
            is_significant_shuffled(end+1, 1)   = shuff.p < 0.05;
            is_significant_linear_fit(end+1, 1)   = shuff.p_lm < 0.05;
        end
        
        mx               = min([u(:).xmin]);
        Mx               = max([u(:).xmax]);
        my               = min([u(:).ymin]);
        My               = max([u(:).ymax]);
        
        for prot_i = 1 : length(u)
            u(prot_i).xlim([mx, Mx]);
            u(prot_i).ylim([my, My]);
        end
        
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                probe_fnames{probe_i}, ...
                clusters(cluster_i).id, ...
                clusters(cluster_i).region_str));
        
        figs.save_fig_to_join();
        
    end
    
    figs.join_figs(sprintf('%s.pdf', probe_fnames{probe_i}));
    figs.clear_figs();
end


results_table = table(recording_name, ...
    cluster_id, ...
    protocol, ...
    r, ...
    R_sq, ...
    fit_slope, ...
    p_val_shuffled, ...
    p_val_linear_fit, ...
    is_significant_shuffled, ...
    is_significant_linear_fit);

writetable(results_table, fullfile(figs.curr_dir, 'summary.csv'))
