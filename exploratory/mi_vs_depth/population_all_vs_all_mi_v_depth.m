function population_all_vs_all_mi_v_depth(experiment, spiking_class)

% experiment              = 'visual_flow';
% spiking_class            = 'FS';   % 'any', 'RS', 'FS'


%%
config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, 'population_mi_vs_depth', spiking_class);

if strcmp(experiment, 'mismatch_nov20+visual_flow')
    probe_fnames            = experiment_details('visual_flow', 'protocols');
    probe_fnames            = [probe_fnames, experiment_details('mismatch_nov20', 'protocols')];
else
    probe_fnames            = experiment_details(experiment, 'protocols');
end

if strcmp(experiment, 'visual_flow')
    protocols           = VisualFlowExperiment.protocol_ids;
    protocol_labels     = VisualFlowExperiment.protocol_label;
elseif strcmp(experiment, 'darkness')
    protocols           = DarknessExperiment.protocol_ids(1:2);
    protocol_labels     = DarknessExperiment.protocol_label(1:2);
elseif strcmp(experiment, 'mismatch_nov20+visual_flow')
    protocol_labels     = {'MVT', 'MV'};
end

x_all                   = cell(length(protocol_labels));
y_all                   = cell(length(protocol_labels));
is_increase             = cell(length(protocol_labels));
p_all                   = cell(length(protocol_labels));
relative_depth          = [];
layer                   = {};
anatomies               = Anatomy.empty();

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.VISp_clusters([], spiking_class);
    anatomies(probe_i)  = Anatomy(data);
    
    exp_obj = get_experiment(data, config);
    
    if strcmp(experiment, 'mismatch_nov20+visual_flow')
        if strcmp(data.experiment_group, 'visual_flow')
            protocols = [1, 2];
        else
            protocols = [2, 4];
        end
    end
    
    for cluster_i = 1 : length(clusters)
        
        [relative_depth(end+1), layer{end+1}] = ...
                anatomies(probe_i).VISp_layer_relative_depth(clusters(cluster_i).distance_from_probe_tip);
        assert(strcmp(layer{end}, clusters(cluster_i).region_str));
        
        for prot_y = 1 : length(protocols)-1
            for prot_x = prot_y+1 : length(protocols)
            
                x           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_x));
                y           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_y));
                
                if strcmp(data.experiment_group, 'mismatch_nov20')
                    x = x(1:10);
                    y = y(1:10);
                end
                
                [x_all{prot_y, prot_x}(end+1, 1), ...
                    y_all{prot_y, prot_x}(end+1, 1), ...
                    p_all{prot_y, prot_x}(end+1, 1), ...
                    is_increase{prot_y, prot_x}(end+1, 1)] = compare_groups_with_signrank(x, y);
                
                fprintf('%s, nnan x: %i\n', probe_fnames{probe_i}, sum(isnan(x)));
                fprintf('%s, nnan y: %i\n', probe_fnames{probe_i}, sum(isnan(y)));
                
            end
        end
    end
end

% average the anatomy
avg_anatomy = AverageAnatomy(anatomies);
[boundaries, cluster_positions] = ...
    avg_anatomy.mi_vs_depth_positions(relative_depth, layer);


%%
h_fig                   = figs.a4figure();
array                   = PlotArray(6, 6);
mi                      = MIDepthPlot.empty();

for prot_y = 1 : length(protocols)-1
    for prot_x = prot_y+1 : length(protocols)
        
        sp_idx      = (prot_y-1)*length(protocols) + prot_x;
        pos         = array.get_position(sp_idx);
        h_ax        = axes('units', 'centimeters', 'position', pos);
        
        mi(end+1)   = MIDepthPlot(x_all{prot_y, prot_x}, ...
                            y_all{prot_y, prot_x}, ...
                            p_all{prot_y, prot_x}, ...
                            is_increase{prot_y, prot_x}, ...
                            cluster_positions, boundaries, ...
                            avg_anatomy.VISp_layers, h_ax);
        
        mi(end).plot();
                        
        if prot_y == 1 && prot_x == 2 
            mi(end).print_layers('left');
            mi(end).xlabel('MI');
        end
        
        str = [protocol_labels{prot_y}, '/', protocol_labels{prot_x}];
        mi(end).title(str);
    end
end


FigureTitle(h_fig, 'population, all vs. all');
figs.save_fig('population_all_vs_all.pdf');