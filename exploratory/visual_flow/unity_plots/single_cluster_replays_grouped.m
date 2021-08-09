% general config info
config                  = RC2AnalysisConfig();

% where to save figure
figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('visual_flow', 'single_cluster_unity_plots', 'replay_conditions_grouped');

% get details of the experiment
probe_fnames            = experiment_details('visual_flow', 'protocols');

plot_array             = PlotArray(1, 3);

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    vf                  = VisualFlowExperiment(data, config);
    clusters            = data.selected_clusters;
    
    for cluster_i = 1 : length(clusters)
        
        h_fig           = figs.a4figure();
        u               = UnityPlotSingleCluster.empty();
        
        for prot_i = 1 : 2
            
            pos         = plot_array.get_position(prot_i);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            if prot_i == 1
                
                x           = [vf.trial_stationary_fr(clusters(cluster_i).id, 3); ...
                               vf.trial_stationary_fr(clusters(cluster_i).id, 4)];
                y           = [vf.trial_motion_fr(clusters(cluster_i).id, 3); ...
                               vf.trial_motion_fr(clusters(cluster_i).id, 4)];
                
            elseif prot_i == 2
            
                x           = [vf.trial_stationary_fr(clusters(cluster_i).id, 5); ...
                               vf.trial_stationary_fr(clusters(cluster_i).id, 6)];
                y           = [vf.trial_motion_fr(clusters(cluster_i).id, 5); ...
                               vf.trial_motion_fr(clusters(cluster_i).id, 6)];
                
            elseif prot_i == 3
                
                x           = [vf.trial_motion_fr(clusters(cluster_i).id, 5); ...
                               vf.trial_motion_fr(clusters(cluster_i).id, 6)];
                y           = [vf.trial_motion_fr(clusters(cluster_i).id, 3); ...
                               vf.trial_motion_fr(clusters(cluster_i).id, 4)];
                
            end
            
            u(prot_i)   = UnityPlotSingleCluster(x, y, h_ax);
            
            if prot_i == 1
                u(prot_i).xlabel('Stationary (Hz)');
                u(prot_i).ylabel('Motion (Hz)');
            elseif prot_i == 2
                u(prot_i).xlabel('Stationary (Hz)');
                u(prot_i).ylabel('Motion (Hz)');
            elseif prot_i == 3
                u(prot_i).xlabel('V (MVT & MV)');
                u(prot_i).ylabel('VT (MVT & MV)');
            end
            
            if prot_i == 1
                u(prot_i).title('VT (MVT & MV)');
            elseif prot_i == 2
                u(prot_i).title('V (MVT & MV)');
            end
            
        end
        
        m               = min([u(:).min]);
        M               = max([u(:).max]);
        
        for prot_i = 1 : length(u)
            u(prot_i).xlim([m, M])
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
