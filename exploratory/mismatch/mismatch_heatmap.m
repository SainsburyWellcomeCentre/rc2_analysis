% heatmap showing firing rate of cells in response to mismatch onset
clearvars -except data

spiking_class   = 'any';

baseline_t      = [-1, 0];
heatmap_padding = [-1, 1];


%%
config                  = config_rc2_analysis();

figs                    = RC2Figures(config);
figs.save_on            = false;
figs.set_figure_subdir('mismatch_nov20', spiking_class);

probe_ids           = experiment_details('mismatch_nov20');

protocol_ids            = 1:4;
protocol_labels         = MismatchExperiment.protocol_label(protocol_ids);

all_cluster_fr          = cell(length(protocol_labels), 1);
response_magnitude      = cell(length(protocol_labels), 1);
store_spiking_class     = cell(length(protocol_labels), 1);

for rec_i = 1 : length(probe_ids)
    
    this_data           = get_data_for_recording_id(data, probe_ids{rec_i});
    clusters            = this_data.VISp_clusters([], spiking_class);
    
    if isempty(clusters)
        continue
    end
    
    exp_obj             = get_experiment(this_data, config);
    
    for prot_i = 1 : length(protocol_ids)
        
        for cluster_i = 1 : length(clusters)
            
            [spike_rate, t] = exp_obj.average_firing_around_mismatch_by_protocol(clusters(cluster_i), protocol_ids(prot_i), heatmap_padding);
            response = exp_obj.average_mismatch_response_by_protocol(clusters(cluster_i), protocol_ids(prot_i));
            
            all_cluster_fr{prot_i}(:, end+1) = spike_rate;
            response_magnitude{prot_i}(end+1, 1) = response;
            store_spiking_class{prot_i}(end+1, 1) = clusters(cluster_i).duration < 0.45;
        end
        
    end
end

% reorder the heatmaps
[~, order_idx] = sort(response_magnitude{2});
heatmaps = cellfun(@(x)(x(:, order_idx)'), all_cluster_fr, 'uniformoutput', false);

% subtract baseline
baseline_idx = t >= baseline_t(1) & t < baseline_t(2);
heatmaps = cellfun(@(x)(bsxfun(@minus, x, mean(x(:, baseline_idx), 2))), heatmaps, 'uniformoutput', false);



%% heatmap
load('r2bcmap.mat', 'map')
figs.a4figure();

clims                   = [-8, 8];
plot_array              = PlotArray(3, 2);
n_clusters              = size(all_cluster_fr{1}, 2);

for prot_i = 1 : length(protocol_ids)
    
    % create an axis in a particular position
    pos                 = plot_array.get_position(prot_i);
    h_ax                = axes('units', 'centimeters', 'position', pos); hold on;
    
    % plot the heatmap
    h_im                = imagesc(h_ax, heatmaps{prot_i});
    
    % set the x-data to time
    set(h_im, 'xdata', t);
    % give the data the red/blue colormap
    colormap(map)
    
    % set the color limits to 8Hz above and below
    set(gca, 'clim', clims);
    
    % line at 0 time
    line([0, 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'k');
    line([0.2, 0.2], get(gca, 'ylim'), 'linestyle', '--', 'color', 'k');

    set(h_ax, 'ytick', [1, n_clusters], 'xlim', [-1, 1]);
    
    title(protocol_labels{prot_i});
    
    if prot_i == length(protocol_ids)
        pos = get(h_ax, 'position');
        new_ax = axes('units', 'centimeters', 'position', pos, 'positionconstraint', 'innerposition');
        set(new_ax, 'clim', clims);
        h = colorbar;
        set(get(h, 'label'), 'string', '\DeltaFR (Hz)');
        set(new_ax, 'position', pos);
        set(new_ax, 'visible', 'off');
        set(h, 'ytick', [-8, 0, 8], 'tickdir', 'out');
        xlabel(h_ax, 'time (ms)');
    end
    
    set(h_ax, 'layer', 'top');
end


figs.save_fig('mismatch_heatmap.pdf');



%% population average

h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);

for prot_i = 1 : length(protocol_ids)
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    hold on
    
    population_average = mean(heatmaps{prot_i}, 1);
    population_sem = std(heatmaps{prot_i}, [], 1) ./ ...
                     sqrt(sum(~isnan(heatmaps{prot_i}), 1));
    
    h_ax = plot(common_t, population_average, 'k');
    plot(common_t, population_average + population_sem, 'color', [0.6, 0.6, 0.6]);
    plot(common_t, population_average - population_sem, 'color', [0.6, 0.6, 0.6]);
    
    ylabel('\Delta FR (Hz)')
    xlabel('Time (s)');
    
    ylim([-2, 6]);
    line([0, 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'k');
    title(protocol_labels{prot_i});
end

figs.save_fig('mismatch_population_average.pdf');




%% population average by cell class

h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);


for prot_i = 1 : length(protocol_ids)
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    hold on
    
    narrow_idx = store_spiking_class{prot_i};
    wide_idx = ~narrow_idx;
    
    fprintf('%i narrow spiking cells, %i wide spiking cells\n', ...
        sum(narrow_idx), ...
        sum(wide_idx));
    
    narrow_population_average = mean(heatmaps{prot_i}(narrow_idx, :), 1);
    wide_population_average = mean(heatmaps{prot_i}(wide_idx, :), 1);
    
    narrow_population_sem = std(heatmaps{prot_i}(narrow_idx, :), [], 1) ./ ...
                     sqrt(sum(~isnan(heatmaps{prot_i}(narrow_idx, :)), 1));
    wide_population_sem = std(heatmaps{prot_i}(wide_idx, :), [], 1) ./ ...
                     sqrt(sum(~isnan(heatmaps{prot_i}(wide_idx, :)), 1));
    
    
    plot(common_t, wide_population_average + wide_population_sem, 'color', 'k');
    plot(common_t, wide_population_average - wide_population_sem, 'color', 'k');
    plot(common_t, narrow_population_average + narrow_population_sem, 'color', [0.5, 0.5, 0.5]);
    plot(common_t, narrow_population_average - narrow_population_sem, 'color', [0.5, 0.5, 0.5]);
    
    plot(common_t, wide_population_average, 'k');
    plot(common_t, narrow_population_average, 'color', [0.5, 0.5, 0.5]);
    
    ylabel('\Delta FR (Hz)')
    xlabel('Time (s)');
    
    ylim([-3, 10]);
    line([0, 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'k');
    title(protocol_labels{prot_i});
end

figs.save_fig('mismatch_population_average_cell_class.pdf');
