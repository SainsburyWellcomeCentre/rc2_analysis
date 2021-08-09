function population_paired(data, experiment, spiking_class)
%%

config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = false;
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
    protocols           = DarknessExperiment.protocol_ids;
    protocol_labels     = DarknessExperiment.protocol_label;
elseif strcmp(experiment, 'passive')
    protocols           = PassiveExperiment.protocol_ids;
    protocol_labels     = PassiveExperiment.protocol_label;
elseif strcmp(experiment, 'head_tilt')
    protocols           = HeadTiltExperiment.protocol_ids;
    protocol_labels     = HeadTiltExperiment.protocol_label;
elseif strcmp(experiment, 'mismatch_nov20')
    protocols           = MismatchExperiment.protocol_ids;
    protocol_labels     = MismatchExperiment.protocol_label;
elseif strcmp(experiment, 'mismatch_nov20+visual_flow')
     protocol_labels     = {'MVT', 'MV'};
end

cluster_id              = cell(length(protocol_labels), 1);
probe_name              = cell(length(protocol_labels), 1);
protocol_id             = cell(length(protocol_labels), 1);
store_spiking_class     = cell(length(protocol_labels), 1);


x_all                   = cell(length(protocol_labels), 1);
y_all                   = cell(length(protocol_labels), 1);
is_increase             = cell(length(protocol_labels), 1);
p_all                   = cell(length(protocol_labels), 1);


for probe_i = 1 : length(probe_fnames)
    
%     data                = load_formatted_data(probe_fnames{probe_i}, config);
    this_data           = get_data_for_recording_id(data, probe_fnames{probe_i});
    clusters            = this_data.VISp_clusters([], spiking_class);
    
    if isempty(clusters)
        continue
    end
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(this_data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(this_data, config);
    elseif strcmp(experiment, 'passive')
        exp_obj         = PassiveExperiment(this_data, config);
    elseif strcmp(experiment, 'head_tilt')
        exp_obj         = HeadTiltExperiment(this_data, config);
    elseif strcmp(experiment, 'mismatch_nov20')
        exp_obj         = MismatchExperiment(this_data, config);
    elseif strcmp(experiment, 'mismatch_nov20+visual_flow')
        exp_obj         = get_experiment(this_data, config);
        if strcmp(this_data.experiment_type, 'visual_flow')
            protocols = [1, 2];
        else
            protocols = [2, 4];
        end
    end
    
    for prot_i = 1 : length(protocols)
        
        for cluster_i = 1 : length(clusters)
            
            x           = exp_obj.trial_stationary_fr(clusters(cluster_i).id, protocols(prot_i));
            y           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_i));
            
            fprintf('%s: length, %i; n_nan: %i\n', probe_fnames{probe_i}, length(x), sum(isnan(x)));
            
            % store the cluster
            probe_name{prot_i}{end+1, 1} = probe_fnames{probe_i};
            protocol_id{prot_i}(end+1, 1) = protocols(prot_i);
            cluster_id{prot_i}(end+1, 1) = clusters(cluster_i).id;
            
            clust_obj = Cluster(clusters(cluster_i));
            store_spiking_class{prot_i}{end+1} = clust_obj.spiking_class;
            
            [x_all{prot_i}(end+1, 1), ...
             y_all{prot_i}(end+1, 1), ...
             p_all{prot_i}(end+1, 1), ...
             is_increase{prot_i}(end+1, 1)] = compare_groups_with_signrank(x, y);
            
        end
    end
end


csvs.create_table(  'probe_name',       cat(1, probe_name{:}), ...
                    'cluster_id',       cat(1, cluster_id{:}), ...
                    'protocol_id',      cat(1, protocol_id{:}), ...
                    'stationary_fr',    cat(1, x_all{:}), ...
                    'motion_fr',        cat(1, y_all{:}), ...
                    'p_val_signrank',   cat(1, p_all{:}), ...
                    'is_increase',      cat(1, is_increase{:}))
csvs.save('population_motion_vs_stationary');


%%
h_fig                   = figs.a4figure();
plot_array              = PlotArray(4, 2);
u                       = UnityPlotPopulation.empty();

for prot_i = 1 : length(protocols)
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    u(end+1)   = UnityPlotPopulation(x_all{prot_i}, ...
        y_all{prot_i}, ...
        p_all{prot_i}, ...
        is_increase{prot_i}, ...
        h_ax);
    
    u(end).plot();
    
    u(end).xlabel('Stationary (Hz)');
    u(end).ylabel('Motion (Hz)');
    
    u(end).title(protocol_labels{prot_i});
    u(end).add_histogram(1);
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for prot_i = 1 : length(u)
    u(prot_i).xlim([m, M])
end

FigureTitle(h_fig, 'population, motion vs. stationary');
figs.save_fig('population_motion_vs_stationary.pdf');



%% color by cell class
h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);
u                       = UnityPlotPopulation.empty();

for prot_i = 1 : length(protocols)
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    hold on;
    
    narrow_idx = strcmp(store_spiking_class{prot_i}, 'FS');
    wide_idx = ~narrow_idx;
    
    fprintf('%i narrow spiking cells, %i wide spiking cells\n', ...
        sum(narrow_idx), ...
        sum(wide_idx));
    
    scatter(h_ax, x_all{prot_i}(wide_idx), y_all{prot_i}(wide_idx), scatterball_size(1), 'k', 'fill');
    scatter(h_ax, x_all{prot_i}(narrow_idx), y_all{prot_i}(narrow_idx), scatterball_size(1), [0.5, 0.5, 0.5], 'fill');
    
    line([m, M], [m, M], 'color', 'k', 'linestyle', '--');
    
    set(h_ax, 'xlim', [m, M], 'ylim', [m, M]);
    set(h_ax, 'xtick', 0:20:60, 'ytick', 0:20:60, 'fontsize', 8);
    
    xlabel('Stationary (Hz)');
    ylabel('Motion (Hz)');
    
    title(protocol_labels{prot_i});
end

FigureTitle(h_fig, 'population, motion vs. stationary, spiking class');
figs.save_fig('population_motion_vs_stationary_spiking_class.pdf');