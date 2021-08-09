experiment              = 'passive';  % 'darkness' or 'visual_flow'
plot_type               = 'box_plots';  % 'box_plots' or 'unity_plots'

config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, sprintf('single_cluster_%s', plot_type), ...
                       'all_conditions');

probe_fnames            = experiment_details(experiment, 'protocols');

plot_array             = PlotArray(3, 2);

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.selected_clusters;
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(data, config);
    elseif strcmp(experiment, 'passive')
        exp_obj         = PassiveExperiment(data, config);
    end
    
    for cluster_i = 1 : length(clusters)
        
        h_fig           = figs.a4figure();
        
        switch plot_type
            case 'unity_plots'
                u = UnityPlotSingleCluster.empty();
            case 'box_plots'
                u = PairedBoxPlot.empty();
        end
        
        for prot_i = 1 : length(exp_obj.protocol_ids)
            
            pos         = plot_array.get_position(prot_i);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            x           = exp_obj.trial_stationary_fr(clusters(cluster_i).id, prot_i);
            y           = exp_obj.trial_motion_fr(clusters(cluster_i).id, prot_i);
            
            switch plot_type
                case 'unity_plots'
                    
                    u(prot_i)   = UnityPlotSingleCluster(x, y, h_ax);
                    
                    if prot_i == length(exp_obj.protocol_ids)
                        u(prot_i).xlabel('Stationary (Hz)');
                        u(prot_i).ylabel('Motion (Hz)');
                    end
                    
                case 'box_plots'
                    
                    u(prot_i)   = PairedBoxPlot(x, y, h_ax);
                    
                    u(prot_i).xticklabel('Stationary', 'Motion');
                    u(prot_i).ylabel('(Hz)');        
            end
            
            u(prot_i).title(exp_obj.protocol_label{prot_i});
        end
        
        m               = min([u(:).min]);
        M               = max([u(:).max]);
        
        for prot_i = 1 : length(u)
            u(prot_i).ylim([m, M]);
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
