function population_all_vs_all_paired(data, experiment, spiking_class)

config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, 'population_unity_plots', spiking_class);

csvs                    = CSVManager(config);
csvs.save_on            = false;
csvs.set_csv_fulldir(figs.curr_dir);


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
    protocols = 1:3;
elseif strcmp(experiment, 'mismatch_nov20')
    protocols           = MismatchExperiment.protocol_ids;
    protocol_labels     = MismatchExperiment.protocol_label;
elseif strcmp(experiment, 'mismatch_nov20+visual_flow')
     protocol_labels     = {'MVT', 'MV'};
end


cluster_id              = cell(length(protocol_labels));
probe_name              = cell(length(protocol_labels));
protocol_x_id           = cell(length(protocol_labels));
protocol_y_id           = cell(length(protocol_labels));

x_all                   = cell(length(protocol_labels));
y_all                   = cell(length(protocol_labels));
is_increase             = cell(length(protocol_labels));
p_all                   = cell(length(protocol_labels));
store_spiking_class     = cell(length(protocol_labels));

for probe_i = 1 : length(probe_fnames)
    
%     data                = load_formatted_data(probe_fnames{probe_i}, config);
    this_data           = get_data_for_recording_id(data, probe_fnames{probe_i});
    clusters            = this_data.VISp_clusters([], spiking_class);
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(this_data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(this_data, config);
        if any(strcmp(probe_fnames{probe_i}, {'CA_176_1_rec1_rec2_rec3', 'CA_176_3_rec1_rec2_rec3'}))
            protocols           = DarknessExperiment.protocol_ids(1:3);
            protocol_labels     = DarknessExperiment.protocol_label(1:3);
        else
            protocols           = DarknessExperiment.protocol_ids([1:2, 4]);
            protocol_labels     = DarknessExperiment.protocol_label([1:2, 4]);
        end
    elseif strcmp(experiment, 'mismatch_nov20')
        exp_obj         = MismatchExperiment(this_data, config);
    elseif strcmp(experiment, 'mismatch_nov20+visual_flow')
        exp_obj         = get_experiment(this_data, config);
        if strcmp(this_data.experiment_group, 'visual_flow')
            protocols = [1, 2];
        else
            protocols = [2, 4];
        end
    end
    
    
    for prot_y = 1 : length(protocols)-1
        for prot_x = prot_y+1 : length(protocols)
            
            for cluster_i = 1 : length(clusters)
                
                x           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_x));
                y           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_y));
                
                if strcmp(this_data.experiment_group, 'mismatch_nov20')
                    x = x(1:10);
                    y = y(1:10);
                end
                
                % store the cluster
                probe_name{prot_y, prot_x}{end+1, 1} = probe_fnames{probe_i};
                protocol_x_id{prot_y, prot_x}(end+1, 1) = protocols(prot_x);
                protocol_y_id{prot_y, prot_x}(end+1, 1) = protocols(prot_y);
                cluster_id{prot_y, prot_x}(end+1, 1) = clusters(cluster_i).id;
                store_spiking_class{prot_y, prot_x}(end+1, 1) = clusters(cluster_i).duration < 0.45;
                
                [x_all{prot_y, prot_x}(end+1, 1), ...
                    y_all{prot_y, prot_x}(end+1, 1), ...
                    p_all{prot_y, prot_x}(end+1, 1), ...
                    is_increase{prot_y, prot_x}(end+1, 1)] = compare_groups_with_signrank(x, y);
            end
        end
    end
end


csvs.create_table(  'probe_name',       cat(1, probe_name{:}), ...
                    'cluster_id',       cat(1, cluster_id{:}), ...
                    'protocol_x_id',    cat(1, protocol_x_id{:}), ...
                    'protocol_y_id',    cat(1, protocol_y_id{:}), ...
                    'x_fr',             cat(1, x_all{:}), ...
                    'y_fr',             cat(1, y_all{:}), ...
                    'p_val_signrank',   cat(1, p_all{:}), ...
                    'is_increase',      cat(1, is_increase{:}))
csvs.save('population_all_vs_all');



if strcmp(experiment, 'darkness')
    % there are two different protocols for the translation 'darkness' experiments
    %   we just give them the label T
    protocol_labels{3} = 'T';
end



%%
n_protocols             = length(protocols);
h_fig                   = figs.a4figure();
plot_array              = PlotArray(n_protocols, n_protocols);
u                       = UnityPlotPopulation.empty();

for prot_y = 1 : n_protocols-1
    for prot_x = prot_y+1 : n_protocols
        
        sp_idx      = (prot_y-1)*n_protocols + prot_x;
        pos         = plot_array.get_position(sp_idx);
        h_ax        = axes('units', 'centimeters', 'position', pos);
        
        u(end+1)   = UnityPlotPopulation(x_all{prot_y, prot_x}, ...
            y_all{prot_y, prot_x}, ...
            p_all{prot_y, prot_x}, ...
            is_increase{prot_y, prot_x}, ...
            h_ax);
        
        u(end).plot();
        
        u(end).xlabel(protocol_labels{prot_x});
        u(end).ylabel(protocol_labels{prot_y});
        u(end).add_histogram(1);
    end
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for prot_i = 1 : length(u)
    u(prot_i).xlim([m, M])
end

FigureTitle(h_fig, 'population, all vs. all');
figs.save_fig('population_all_vs_all.pdf');


%%
n_protocols             = length(protocols);
h_fig                   = figs.a4figure();
plot_array              = PlotArray(n_protocols, n_protocols);
u                       = UnityPlotPopulation.empty();

for prot_y = 1 : n_protocols-1
    for prot_x = prot_y+1 : n_protocols
        
        sp_idx      = (prot_y-1)*n_protocols + prot_x;
        pos         = plot_array.get_position(sp_idx);
        h_ax        = axes('units', 'centimeters', 'position', pos);
        
        narrow_idx = strcmp(store_spiking_class{prot_y, prot_x}, 'FS');
        wide_idx = ~narrow_idx;
        
        fprintf('%i narrow spiking cells, %i wide spiking cells\n', ...
            sum(narrow_idx), ...
            sum(wide_idx));
        
        scatter(h_ax, x_all{prot_y, prot_x}(wide_idx), y_all{prot_y, prot_x}(wide_idx), scatterball_size(1), 'k', 'fill');
        scatter(h_ax, x_all{prot_y, prot_x}(narrow_idx), y_all{prot_y, prot_x}(narrow_idx), scatterball_size(1), [0.5, 0.5, 0.5], 'fill');
        
        line([m, M], [m, M], 'color', 'k', 'linestyle', '--');
        
        set(h_ax, 'xlim', [m, M], 'ylim', [m, M]);
        set(h_ax, 'xtick', 0:20:60, 'ytick', 0:20:60, 'fontsize', 8);
        
        xlabel(protocol_labels{prot_x});
        ylabel(protocol_labels{prot_y});
    end
end

FigureTitle(h_fig, 'population, all vs. all, spiking class');
figs.save_fig('population_all_vs_all_spiking_class.pdf');