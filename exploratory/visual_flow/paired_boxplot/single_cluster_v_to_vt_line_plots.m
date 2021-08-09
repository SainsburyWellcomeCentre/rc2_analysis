config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('visual_flow', datestr(now, 'yyyymmdd'), ...
                       'single_cluster_V_to_VT_line_plot');

probe_fnames            = experiment_details('visual_flow', 'protocols');

plot_array             = PlotArray(1, 1);
plot_array.nx_total = 1;
plot_array.ny_total = 1;
plot_array.ax_size_cm = 10;

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.selected_clusters;
    
    exp_obj             = VisualFlowExperiment(data, config);
    
    V_trials            = exp_obj.trials_of_type('ReplayOnly');
    VT_trials           = exp_obj.trials_of_type('StageOnly');
    
    clear V_replayed_trials
    for trial_i = 1 : length(V_trials)    
        V_replayed_trials(trial_i) = exp_obj.get_replayed_trial(V_trials(trial_i));
    end
    
    clear VT_replayed_trials
    for trial_i = 1 : length(VT_trials)    
        VT_replayed_trials(trial_i) = exp_obj.get_replayed_trial(VT_trials(trial_i));
    end
    
    V_stationary_fr = nan(length(V_trials), length(clusters));
    V_motion_fr = nan(length(V_trials), length(clusters));
    VT_stationary_fr = nan(length(V_trials), length(clusters));
    VT_motion_fr = nan(length(V_trials), length(clusters));
    
    
    MV_trials            = exp_obj.trials_of_type('EncoderOnly');
    MVT_trials           = exp_obj.trials_of_type('Coupled');
    
    MV_stationary_fr    = nan(length(MV_trials), length(clusters));
    MV_motion_fr        = nan(length(MV_trials), length(clusters));
    MVT_stationary_fr   = nan(length(MV_trials), length(clusters));
    MVT_motion_fr       = nan(length(MV_trials), length(clusters));
    
    
    
    for cluster_i = 1 : length(clusters)
        
        for trial_i = 1 : length(V_trials)
            
            idx = exp_obj.svm_table.cluster_id == clusters(cluster_i).id & ...
                exp_obj.svm_table.trial_id == V_trials(trial_i).id;
            
            assert(sum(idx) == 1);
            
            V_stationary_fr(trial_i, cluster_i) = exp_obj.svm_table.stationary_firing_rate(idx);
            V_motion_fr(trial_i, cluster_i) = exp_obj.svm_table.motion_firing_rate(idx);
            
            VT_trial_idx = find([VT_replayed_trials(:).id] == V_replayed_trials(trial_i).id);
            
            idx = exp_obj.svm_table.cluster_id == clusters(cluster_i).id & ...
                exp_obj.svm_table.trial_id == VT_trials(VT_trial_idx).id;
            
            assert(sum(idx) == 1);
            
            VT_stationary_fr(trial_i, cluster_i) = exp_obj.svm_table.stationary_firing_rate(idx);
            VT_motion_fr(trial_i, cluster_i) = exp_obj.svm_table.motion_firing_rate(idx);
        end
        
        
        for trial_i = 1 : length(MV_trials)
            
            idx = exp_obj.svm_table.cluster_id == clusters(cluster_i).id & ...
                exp_obj.svm_table.trial_id == MV_trials(trial_i).id;
            
            MV_stationary_fr(trial_i, cluster_i) = exp_obj.svm_table.stationary_firing_rate(idx);
            
            idx = exp_obj.svm_table.cluster_id == clusters(cluster_i).id & ...
                exp_obj.svm_table.trial_id == MVT_trials(trial_i).id;
            
            MV_stationary_fr(trial_i, cluster_i) = exp_obj.svm_table.stationary_firing_rate(idx);
        end
    end
    
    %%
    for cluster_i = 1 : length(clusters)
        
        X = [V_stationary_fr(:, cluster_i), V_motion_fr(:, cluster_i), ...
            VT_stationary_fr(:, cluster_i), VT_motion_fr(:, cluster_i)];
        
        h_fig           = figs.a4figure();
        pos         = plot_array.get_position(1);
        h_ax        = axes('units', 'centimeters', 'position', pos);
        
        u = ExtendedPairedPlot(X, h_ax);
        u.ylabel('FR (Hz)');
        str = {'Stationary', 'Motion', 'Stationary', 'Motion'; 'V', 'V', 'VT', 'VT'};
        str = strtrim(sprintf('%s\\newline%s\n', str{:}));
        set(u.h_ax, 'xticklabel', str);
        
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                probe_fnames{probe_i}, ...
                clusters(cluster_i).id, ...
                clusters(cluster_i).region_str));
            
        figs.save_fig_to_join();
        
    end
    
    figs.join_figs(sprintf('%s.pdf', probe_fnames{probe_i}));
    figs.clear_figs();
    
end
