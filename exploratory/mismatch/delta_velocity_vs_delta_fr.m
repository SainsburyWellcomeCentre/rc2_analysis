clear all
close all

window_t                = 0.4;
n_sds                   = 1;
protocol                = 2;
display_window          = [-1, 1];


%%
config                  = RC2AnalysisConfig();

protocols               = MismatchExperiment.protocol_ids;
protocol_labels         = MismatchExperiment.protocol_label;

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('mismatch', 'running_change_at_mm', sprintf('%ims_%isd_offset_protocol_%i', 1e3*window_t, n_sds, protocol));

probe_fnames            = experiment_details('mismatch_nov20', 'protocol');

delta_speed             = cell(1, length(probe_fnames));
delta_fr                = cell(1, length(probe_fnames));

store_population_spiking = [];

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.VISp_clusters;
    
    exp_obj             = MismatchExperiment(data, config);
    
    
    for prot_i = protocol  % 1 : length(protocols)  % 
        
        [running, t] = exp_obj.running_around_mismatch_by_protocol(prot_i, display_window);
        
        spike_rate = cell(1, length(clusters));
        for cluster_i = 1 : length(clusters)
            spike_rate{cluster_i} = exp_obj.firing_around_mismatch_by_protocol(clusters(cluster_i), prot_i, display_window);
        end
        
        trials = exp_obj.trials_of_type(prot_i);
        
        assert(length(trials) == size(running, 2));
        
        delta_speed{probe_i} = nan(length(trials), 1);
        delta_fr{probe_i} = nan(length(trials), 1);
        
        for trial_i = 1 : length(trials)
            
            baseline_idx = t > -window_t & t < 0;
            response_idx = t >= 0 & t < window_t;
            
            baseline_speed = mean(running(baseline_idx, trial_i));
            response_speed = mean(running(response_idx, trial_i));
            
            delta_speed{probe_i}(trial_i) = response_speed - baseline_speed;
            
            population_spiking = cellfun(@(x)(x(:, trial_i)), spike_rate, 'uniformoutput', false);
            
            population_spiking = mean([population_spiking{:}], 2);
            
            baseline_fr = mean(population_spiking(baseline_idx));
            response_fr = mean(population_spiking(response_idx));
            
            delta_fr{probe_i}(trial_i) = response_fr - baseline_fr;
        end
    end
    
    % correlate change in speed and change in firing rate
    [r, p] = corr(delta_speed{probe_i}, delta_fr{probe_i});
    
    figure
    scatter(delta_speed{probe_i}, delta_fr{probe_i}, [], 'k', 'fill');
    
    title(sprintf('%s, r=%.2f, p=%.2f', probe_fnames{probe_i}, r, p), 'interpreter', 'none');
    
    xlabel('\Delta running speed (cm/s)');
    ylabel('\Delta FR (Hz)');
    
    Mx = max(abs(get(gca, 'xlim')));
    My = max(abs(get(gca, 'ylim')));
    
    set(gca, 'xlim', [-Mx, Mx], 'ylim', [-My, My], 'box', 'off');
end


%%

cols = lines(length(probe_fnames));
figure
hold on

for probe_i = 1 : length(probe_fnames)
    
    scatter(delta_speed{probe_i}, delta_fr{probe_i}, [], cols(probe_i, :), 'fill');
end

legend(probe_fnames, 'interpreter', 'none');

xlabel('\Delta running speed (cm/s)');
ylabel('\Delta FR (Hz)');

Mx = max(abs(get(gca, 'xlim')));
My = max(abs(get(gca, 'ylim')));

set(gca, 'xlim', [-Mx, Mx], 'ylim', [-My, My], 'box', 'off');
