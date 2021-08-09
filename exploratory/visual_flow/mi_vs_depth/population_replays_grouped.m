config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('visual_flow', 'population_mi_vs_depth');

probe_fnames            = experiment_details('visual_flow', 'protocols');

protocols               = VisualFlowExperiment.protocol_ids;

x_all                   = cell(length(protocols), 1);
y_all                   = cell(length(protocols), 1);
p_all                   = cell(length(protocols), 1);
relative_depth          = [];
layer                   = {};
anatomies               = Anatomy.empty();

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    vf                  = VisualFlowExperiment(data, config);
    clusters            = data.VISp_clusters;
    anatomies(probe_i)  = Anatomy(data);
    
    for cluster_i = 1 : length(clusters)
        
        [relative_depth(end+1), layer{end+1}] = ...
                anatomies(probe_i).VISp_layer_relative_depth(clusters(cluster_i).distance_from_probe_tip);
        assert(strcmp(layer{end}, clusters(cluster_i).region_str));
        
        for prot_i = 1 : 3
            
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
            
            x_all{prot_i}(end+1) = nanmedian(x);
            y_all{prot_i}(end+1) = nanmedian(y);
            p_all{prot_i}(end+1) = signrank(x, y);
            
        end
    end
end

% average the anatomy
avg_anatomy = AverageAnatomy(anatomies);
[boundaries, cluster_positions] = ...
    avg_anatomy.mi_vs_depth_positions(relative_depth, layer);

%%
h_fig                   = figs.a4figure();
array                   = PlotArray(1, 3);
mi                      = MIDepthPlot.empty();

for prot_i = 1 : 3
    
    pos         = array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'innerposition', pos);
    
    mi(end+1)   = MIDepthPlot(x_all{prot_i}, ...
                            y_all{prot_i}, ...
                            p_all{prot_i}, ...
                            cluster_positions, boundaries, ...
                            avg_anatomy.VISp_layers, h_ax);
    
    if prot_i == 1
        mi(end).title('VT (MVT & MV)');
    elseif prot_i == 2
        mi(end).title('V (MVT & MV)');
    elseif prot_i == 3
        mi(end).print_layers();
        mi(end).title({'VT (MVT & MV)/', 'V (MVT & MV)'});
        mi(end).xlabel('MI');
    end
    
end

FigureTitle(h_fig, 'population, replays grouped');
figs.save_fig('population_replays_grouped.pdf');