spiking_class           = 'any';

config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('visual_flow', 'population_unity_plots', spiking_class);

probe_fnames            = experiment_details('visual_flow', 'protocols');

protocols               = VisualFlowExperiment.protocol_ids;

x_all                   = cell(length(protocols), 1);
y_all                   = cell(length(protocols), 1);
is_increase             = cell(length(protocols), 1);
p_all                   = cell(length(protocols), 1);
store_spiking_class     = cell(length(protocols), 1);

for probe_i = 1 : length(probe_fnames)
    
    this_data           = get_data_for_recording_id(data, probe_fnames{probe_i});
    clusters            = this_data.VISp_clusters([], spiking_class);
    
    vf                  = VisualFlowExperiment(this_data, config);
    
    for prot_i = 1 : 3
        
        for cluster_i = 1 : length(clusters)
            
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
            
            store_spiking_class{prot_i}(end+1, 1) = clusters(cluster_i).duration < 0.45;
            
            [x_all{prot_i}(end+1), ...
             y_all{prot_i}(end+1), ...
             p_all{prot_i}(end+1), ...
             is_increase{prot_i}(end+1)] = compare_groups_with_signrank(x, y);
        end
    end
end


%%
h_fig                   = figs.a4figure();
plot_array             = PlotArray(1, 3);

u                       = UnityPlotPopulation.empty();

for prot_i = 1 : 3
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'innerposition', pos);
    
    u(end+1)   = UnityPlotPopulation(x_all{prot_i}, y_all{prot_i}, p_all{prot_i}, is_increase{prot_i}, h_ax);
    
    u(end).plot();
    
    if prot_i == 1
        u(end).xlabel('Stationary (Hz)');
        u(end).ylabel('Motion (Hz)');
        u(end).title('VT (MVT & MV)');
    elseif prot_i == 2
        u(end).xlabel('Stationary (Hz)');
        u(end).ylabel('Motion (Hz)');
        u(end).title('V (MVT & MV)');
    elseif prot_i == 3
        u(end).xlabel('V (MVT & MV)');
        u(end).ylabel('VT (MVT & MV)');
        u(end).title({'VT (MVT & MV)', 'vs. V (MVT & MV)'});
    end
    u(end).add_histogram(1);
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for prot_i = 1 : length(u)
    u(prot_i).xlim([m, M])
end

FigureTitle(h_fig, 'population, replays grouped');
figs.save_fig('population_replays_grouped.pdf');



%%
h_fig                   = figs.a4figure();
plot_array             = PlotArray(1, 3);

u                       = UnityPlotPopulation.empty();

for prot_i = 1 : 3
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'innerposition', pos);
    hold on;
    
    narrow_idx = store_spiking_class{prot_i};
    wide_idx = ~narrow_idx;
    
    fprintf('%i narrow spiking cells, %i wide spiking cells\n', ...
        sum(narrow_idx), ...
        sum(wide_idx));
    
    scatter(h_ax, x_all{prot_i}(wide_idx), y_all{prot_i}(wide_idx), scatterball_size(1), 'k', 'fill');
    scatter(h_ax, x_all{prot_i}(narrow_idx), y_all{prot_i}(narrow_idx), scatterball_size(1), [0.5, 0.5, 0.5], 'fill');
    
    line([m, M], [m, M], 'color', 'k', 'linestyle', '--');
    
    set(h_ax, 'xlim', [m, M], 'ylim', [m, M]);
    set(h_ax, 'xtick', 0:20:60, 'ytick', 0:20:60, 'fontsize', 8);
    
    if prot_i == 1
        xlabel('Stationary (Hz)');
        ylabel('Motion (Hz)');
        title('VT (MVT & MV)');
    elseif prot_i == 2
        xlabel('Stationary (Hz)');
        ylabel('Motion (Hz)');
        title('V (MVT & MV)');
    elseif prot_i == 3
        xlabel('V (MVT & MV)');
        ylabel('VT (MVT & MV)');
        title({'VT (MVT & MV)', 'vs. V (MVT & MV)'});
    end
end

FigureTitle(h_fig, 'population, replays grouped, spiking class');
figs.save_fig('population_replays_grouped_spiking_class.pdf');