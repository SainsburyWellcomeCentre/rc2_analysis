clearvars -except data

% heatmaps showing transition from baseline to motion for all conditions
experiment      = 'mismatch_nov20';
spiking_class   = 'any';

baseline_t      = [-0.4, 0];
response_t      = [0, 0.4];
heatmap_padding = [-1, 1];


%%
config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, 'motion_onset', spiking_class);

if strcmp(experiment, 'mismatch_nov20+visual_flow')
    probe_fnames            = experiment_details('visual_flow', 'protocols');
    probe_fnames            = [probe_fnames, experiment_details('mismatch_nov20', 'protocols')];
else
    probe_fnames            = experiment_details(experiment, 'protocols');
end

if strcmp(experiment, 'visual_flow')
    protocol_ids        = [1, 2, 7, 8];
    protocol_labels     = VisualFlowExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'darkness')
    protocol_ids        = 1:3;
    protocol_labels     = DarknessExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'mismatch_nov20')
    protocol_ids        = [2, 4];
    protocol_labels     = MismatchExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'passive')
    protocol_ids        = 1:3;
    protocol_labels     = PassiveExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'head_tilt')
    protocol_ids        = 1:3;
    protocol_labels     = HeadTiltExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'mismatch_nov20+visual_flow')
    protocol_labels     = {'MVT', 'MV'};
end


% common timebase on which to compute the firing rates
n_samples = ceil(range(heatmap_padding)*10000);
common_t = linspace(heatmap_padding(1), heatmap_padding(2), n_samples)';
baseline_idx = common_t < baseline_t(2) & common_t >= baseline_t(1);


all_cluster_fr  = cell(length(protocol_labels), 1);
p_all           = cell(length(protocol_labels), 1);
is_increase     = cell(length(protocol_labels), 1);
store_spiking_class   = cell(length(protocol_labels), 1);

for probe_i = 1 : length(probe_fnames)
    
%     data                = load_formatted_data(probe_fnames{probe_i}, config);
    this_data           = get_data_for_recording_id(data, probe_fnames{probe_i});
    clusters            = this_data.VISp_clusters([], spiking_class);
    
    if isempty(clusters)
        continue
    end
    
    exp_obj = get_experiment(this_data, config);
    
    if strcmp(experiment, 'mismatch_nov20+visual_flow')
        if strcmp(this_data.experiment_type, 'visual_flow')
            protocol_ids = [1, 2];
        else
            protocol_ids = [2, 4];
        end
    end
    
    
    for prot_i = 1 : length(protocol_ids)
        
        % get time of motion bout onsets for this protocol
        motion_bouts = exp_obj.motion_bouts_by_protocol(protocol_ids(prot_i), true, true);
        
        %
        motion_bouts = motion_bouts([motion_bouts(:).duration] > 2);
        
        if isempty(motion_bouts)
            continue
        end
        
        start_t = [motion_bouts(:).start_time];
        
        for cluster_i = 1 : length(clusters)
            
            clust_obj = Cluster(clusters(cluster_i));
            fr = FiringRate(clusters(cluster_i).spike_times);
            cluster_fr = nan(n_samples, length(motion_bouts));
            
            for bout_i = 1 : length(motion_bouts)
                
                this_t = start_t(bout_i) + common_t; 
                cluster_fr(:, bout_i) = fr.get_convolution(this_t);
            end
            
            avg_fr = mean(cluster_fr, 2);
            
            % subtract baseline
            all_cluster_fr{prot_i}(:, end+1) = avg_fr - mean(avg_fr(baseline_idx));
            store_spiking_class{prot_i}{end+1} = clust_obj.spiking_class;
            
            % get significant clusters
            x           = exp_obj.trial_stationary_fr(clusters(cluster_i).id, protocol_ids(prot_i));
            y           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocol_ids(prot_i));
            
            [~, ~, p_all{prot_i}(end+1, 1), is_increase{prot_i}(end+1, 1)] = compare_groups_with_signrank(x, y);
        end
    end
end


%% heatmap
load('r2bcmap.mat', 'map')
figs.a4figure();

plot_array              = PlotArray(3, 2);
n_clusters              = length(p_all{1});

% rearrange heatmaps
sig_increase                = find(p_all{1} < 0.05 & is_increase{1});
sig_decrease                = find(p_all{1} < 0.05 & ~is_increase{1});
no_change                   = find(p_all{1} >= 0.05);


% time index of response period
response_idx            = common_t >= response_t(1) & common_t < response_t(2);


for prot_i = 1 : length(protocol_ids)
    
    % the delta FR response for each cluster to this protocol
    delta_response      = mean(all_cluster_fr{prot_i}(response_idx, :), 1);
    
    % index of clusters, sorted by magnitude of response
    [~, sorted_idx]     = sort(delta_response, 'ascend');
    
    % create an axis in a particular position
    pos                 = plot_array.get_position(prot_i);
    h_ax                = axes('units', 'centimeters', 'position', pos); hold on;
    
    % get the heatmap
    heatmap             = all_cluster_fr{prot_i}(:, sorted_idx)';
    
    % plot the heatmap
    h_im                = imagesc(h_ax, heatmap);
    
    % set the x-data to time
    set(h_im, 'xdata', common_t);
    
    % give the data the red/blue colormap
    colormap(map)
    
    % set the color limits to 8Hz above and below
    set(gca, 'clim', [-8, 8]);
    
    % line at 0 time
    line([0, 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'k');

    % loop through clusters and print: 1. whether significant response, 2.
    % whether RS or FS
    for clust_i = 1 : n_clusters
        
        cluster_idx = sorted_idx(clust_i);
        
        p = p_all{prot_i}(cluster_idx);
        inc = is_increase{prot_i}(cluster_idx);
        
        if p < 0.05 && inc
            h = fill([-1.2, -1.1, -1.1, -1.2], ...
                     [clust_i-0.5, clust_i-0.5, clust_i+0.5, clust_i+0.5], ...
                     'r', 'edgecolor', 'none');
        elseif p < 0.05 && ~inc
            h = fill([-1.2, -1.1, -1.1, -1.2], ...
                     [clust_i-0.5, clust_i-0.5, clust_i+0.5, clust_i+0.5], ...
                     [30, 144, 255] / 255, 'edgecolor', 'none');
        elseif p >= 0.05
            
        else
            error('wtf?');
        end
        
        
        if strcmp(store_spiking_class{prot_i}(cluster_idx), 'RS')
            h = fill([-1.1, -1, -1, -1.1], ...
                     [clust_i-0.5, clust_i-0.5, clust_i+0.5, clust_i+0.5], ...
                     'k', 'edgecolor', 'none');
        elseif strcmp(store_spiking_class{prot_i}(cluster_idx), 'FS')
            h = fill([-1.1, -1, -1, -1.1], ...
                     [clust_i-0.5, clust_i-0.5, clust_i+0.5, clust_i+0.5], ...
                     [0.6, 0.6, 0.6], 'edgecolor', 'none');
        else
            error('wtf2?');
        end
    end
    
    set(h_ax, 'ytick', [1, n_clusters], 'xlim', [-1.2, 1]);
    
    title(protocol_labels{prot_i});
end

pos         = plot_array.get_position(5);
h_ax        = axes('units', 'centimeters', 'position', pos);
colormap(map)
set(h_ax, 'clim', [-8, 8])
colorbar
axis off

figs.save_fig('motion_onset_heatmap.pdf');



%% population average

h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);

n_clusters              = length(p_all{1});

for prot_i = 1 : length(protocol_ids)
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    hold on
    
    population_average = mean(all_cluster_fr{prot_i}, 2);
    population_sem = std(all_cluster_fr{prot_i}, [], 2) ./ ...
                     sqrt(sum(~isnan(all_cluster_fr{prot_i}), 2));
    
    h_ax = plot(common_t, population_average, 'k');
    plot(common_t, population_average + population_sem, 'color', [0.6, 0.6, 0.6]);
    plot(common_t, population_average - population_sem, 'color', [0.6, 0.6, 0.6]);
    
    line([0, 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'k');
    
    ylabel('\Delta FR (Hz)')
    xlabel('Time (s)');
    
    ylim([-2, 6]);
    title(protocol_labels{prot_i});
end

figs.save_fig('motion_onset_population_average.pdf');




%% population average by cell class

h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);

n_clusters              = length(p_all{1});

for prot_i = 1 : length(protocol_ids)
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    hold on
    
    narrow_idx = strcmp(store_spiking_class{prot_i}, 'FS');
    wide_idx = ~narrow_idx;
    
    fprintf('%i narrow spiking cells, %i wide spiking cells\n', ...
        sum(narrow_idx), ...
        sum(wide_idx));
    
    narrow_population_average = mean(all_cluster_fr{prot_i}(:, narrow_idx), 2);
    wide_population_average = mean(all_cluster_fr{prot_i}(:, wide_idx), 2);
    
    narrow_population_sem = std(all_cluster_fr{prot_i}(:, narrow_idx), [], 2) ./ ...
                     sqrt(sum(~isnan(all_cluster_fr{prot_i}(:, narrow_idx)), 2));
    wide_population_sem = std(all_cluster_fr{prot_i}(:, wide_idx), [], 2) ./ ...
                     sqrt(sum(~isnan(all_cluster_fr{prot_i}(:, wide_idx)), 2));
    
    
    plot(common_t, wide_population_average + wide_population_sem, 'color', 'k');
    plot(common_t, wide_population_average - wide_population_sem, 'color', 'k');
    plot(common_t, narrow_population_average + narrow_population_sem, 'color', [0.5, 0.5, 0.5]);
    plot(common_t, narrow_population_average - narrow_population_sem, 'color', [0.5, 0.5, 0.5]);
    
    plot(common_t, wide_population_average, 'k');
    plot(common_t, narrow_population_average, 'color', [0.5, 0.5, 0.5]);
    
    line([0, 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'k');
    
    ylabel('\Delta FR (Hz)')
    xlabel('Time (s)');
    
    ylim([-3, 10]);
    title(protocol_labels{prot_i});
end

figs.save_fig('motion_onset_population_average_cell_class.pdf');
