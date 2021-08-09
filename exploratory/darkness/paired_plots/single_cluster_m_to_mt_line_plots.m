config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('darkness', datestr(now, 'yyyymmdd'), ...
                       'single_cluster_M_to_MT_line_plot');

probe_fnames            = experiment_details('darkness', 'protocols');

plot_array             = PlotArray(1, 1);
plot_array.nx_total = 1;
plot_array.ny_total = 1;
plot_array.ax_size_cm = 10;

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.selected_clusters;
    
    exp_obj             = DarknessExperiment(data, config);
    
    M_trials            = exp_obj.trials_of_type('EncoderOnly');
    MT_trials           = exp_obj.trials_of_type('Coupled');
    
    assert(length(M_trials) == length(MT_trials));
    
    M_stationary_fr    = nan(length(M_trials), length(clusters));
    M_motion_fr        = nan(length(M_trials), length(clusters));
    MT_stationary_fr   = nan(length(MT_trials), length(clusters));
    MT_motion_fr       = nan(length(MT_trials), length(clusters));
    
    for cluster_i = 1 : length(clusters)
        
        for trial_i = 1 : length(M_trials)
            
            idx = exp_obj.svm_table.cluster_id == clusters(cluster_i).id & ...
                exp_obj.svm_table.trial_id == M_trials(trial_i).id;
            
            M_stationary_fr(trial_i, cluster_i) = exp_obj.svm_table.stationary_firing_rate(idx);
            M_motion_fr(trial_i, cluster_i) = exp_obj.svm_table.motion_firing_rate(idx);
            
            idx = exp_obj.svm_table.cluster_id == clusters(cluster_i).id & ...
                exp_obj.svm_table.trial_id == MT_trials(trial_i).id;
            
            MT_stationary_fr(trial_i, cluster_i) = exp_obj.svm_table.stationary_firing_rate(idx);
            MT_motion_fr(trial_i, cluster_i) = exp_obj.svm_table.motion_firing_rate(idx);
        end
    end
    
    %%
    for cluster_i = 1 : length(clusters)
        
        X = [M_stationary_fr(:, cluster_i),  M_motion_fr(:, cluster_i), ...
             MT_stationary_fr(:, cluster_i), MT_motion_fr(:, cluster_i)];
        
        h_fig       = figs.a4figure();
        pos         = plot_array.get_position(1);
        h_ax        = axes('units', 'centimeters', 'position', pos);
        
        u = ExtendedPairedPlot(X, h_ax);
        u.ylabel('FR (Hz)');
        str = {'Stationary', 'Motion', 'Stationary', 'Motion'; 'M', 'M', 'MT', 'MT'};
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
