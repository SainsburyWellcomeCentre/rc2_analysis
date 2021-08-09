experiment              = 'darkness';

config                  = RC2AnalysisConfig();

probe_fnames            = experiment_details(experiment, 'protocols');

if strcmp(experiment, 'visual_flow')
    protocols           = VisualFlowExperiment.protocol_ids;
    protocol_labels     = VisualFlowExperiment.protocol_label;
elseif strcmp(experiment, 'darkness')
    protocols           = DarknessExperiment.protocol_ids(1:3);
    protocol_labels     = DarknessExperiment.protocol_label(1:3);
end

results                       = cell(1, length(protocols));

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    clusters            = data.VISp_clusters;
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(data, config);
    end
    
    for prot_i = 1 : length(protocols)
        
        for cluster_i = 1 : length(clusters)
            
            [fr, sd, n, x, p_anova, beta, stat_fr, stat_sd, stat_n] = ...
                exp_obj.tuning_curve(clusters(cluster_i).id, prot_i);
            
            stat_rate = exp_obj.trial_stationary_fr(clusters(cluster_i).id, prot_i);
            mot_rate = exp_obj.trial_motion_fr(clusters(cluster_i).id, prot_i);
            
            p_signrank = signrank(stat_rate, mot_rate);
            
            results{prot_i}(end+1, 1) = p_signrank;
            results{prot_i}(end, 2) = p_anova;
            results{prot_i}(end, 3) = mean(fr);
            results{prot_i}(end, 4) = stat_fr;
            results{prot_i}(end, 5) = nanmedian(mot_rate);
            results{prot_i}(end, 6) = nanmedian(stat_rate);
            
            if p_signrank < 0.05 && ...
                    p_anova >= 0.05 && ...
                    mean(fr) > stat_fr
                % tonic increase
                results{prot_i}(end, 7) = 1;    
            elseif p_signrank < 0.05 && ...
                    p_anova >= 0.05 && ...
                    mean(fr) < stat_fr
                % tonic decrease
                results{prot_i}(end, 7) = 2;
            elseif p_anova < 0.05 && beta(1) >= 0 %obj.fr(end) >= obj.fr(1)
                % speed tuned +
                results{prot_i}(end, 7) = 3;
            elseif p_anova < 0.05 && beta(1) < 0 %obj.fr(end) < obj.fr(1)
                % speed tuned -
                results{prot_i}(end, 7) = 4;
            else
                results{prot_i}(end, 7) = 5;
            end
            
            results{prot_i}(end, 8) = beta(1);
        end
    end
end
