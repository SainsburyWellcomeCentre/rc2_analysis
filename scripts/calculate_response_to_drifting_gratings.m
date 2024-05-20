close all

ctl = RC2Analysis();
probe_ids  = ctl.get_probe_ids('visual_flow', ...
                               'mismatch_nov20', ...
                               'mismatch_jul21', ...
                               'mismatch_darkness_oct21', ...
                               'darkness');
mean_clusters_responses = [];
mean_clusters_baselines = [];
for probe_i = 1 : length(probe_ids)
    
    this_probe = probe_ids(probe_i)
    if strcmp(this_probe, 'CAA-1112872_rec1_rec1b_rec2_rec3') || strcmp(this_probe, 'CAA-1112874_rec1_rec2_rec3') ||  ...
       strcmp(this_probe, 'CAA-1114977_rec1_rec2_rec3') || strcmp(this_probe, 'CAA-1114978_rec1_rec2') || ...
        strcmp(this_probe, 'CAA-1114979_rec1_rec2') || strcmp(this_probe, 'CAA-1114980_rec1_rec2') || ...
        strcmp(this_probe, 'CAA-1115689_rec1_rec2') || strcmp(this_probe, 'CAA-1115691_rec1_rec2') || ...
        strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3') || strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3')
        % skip this probes as the photodiode signal is corrupted or no sf_tf
        % session present
        continue
    end

    data = ctl.load_formatted_data(this_probe{1});
    drifting_gratings_session = data.sessions{1, 3};
    
    animal_id = split(this_probe, "_");
    if strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3') || strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3')
        animal_id = join([animal_id(1), animal_id(2), animal_id(3)], '_');
    else
        animal_id = animal_id(1);
    end
    starts = pd_times(animal_id, drifting_gratings_session, "sf_tf");
    
    
    clusters = data.VISp_clusters();

    for clust_i = 1 : length(clusters)
        responses = cell(1, length(starts));
        baselines = cell(1, 6);
       
        % Identify total spikes in the time interval of interest
        baseline_counter = 1;
        baselines{1, baseline_counter} = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(1) : ...
                                                                                   drifting_gratings_session.probe_t(starts(1))));
        deltas = [];
        for i = 1 : length(starts) - 1
            delta = starts(i + 1) - starts(i);
            if delta < 50000
                deltas(end+1) = delta;
                responses{1, i} = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i)) : ...
                                                                            drifting_gratings_session.probe_t(starts(i+1))));
            else
                mean_duration = cast(mean(deltas), "uint32");
                responses{1, i} = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i)) : ...
                                                                            drifting_gratings_session.probe_t(starts(i) + mean_duration)));
                baseline_counter = baseline_counter + 1;
                baselines{1, baseline_counter} = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(i) + mean_duration) : ...
                                                                                           drifting_gratings_session.probe_t(starts(i+1))));
            end
        end
        baseline_counter = baseline_counter + 1;
        baselines{1, baseline_counter} = mean(clusters(clust_i).fr.get_convolution(drifting_gratings_session.probe_t(starts(end)) : ...
                                                                                   drifting_gratings_session.probe_t(starts(end) + 3 * mean_duration)));        
        mean_clusters_responses(end+1) = mean(responses{1});
        mean_clusters_baselines(end+1) = mean(baselines{1});
   end
end

mean(mean_clusters_responses)
std(mean_clusters_responses)

mean(mean_clusters_baselines)
std(mean_clusters_baselines)


% Lineplot
figure(1)
hold on
for clust_i = 1 : length(mean_clusters_responses)
    scatter(1, mean_clusters_baselines(clust_i), scatterball_size(1), 'blue', 'o');
    scatter(2, mean_clusters_responses(clust_i), scatterball_size(1), 'blue', 'o');
    plot([1 2], [mean_clusters_baselines(clust_i), mean_clusters_responses(clust_i)], 'blue', 'LineWidth', 0.1);
end
scatter(1, mean(mean_clusters_baselines), scatterball_size(3), 'black');
scatter(2, mean(mean_clusters_responses), scatterball_size(3), 'black');
plot([1 2], [mean(mean_clusters_baselines), mean(mean_clusters_responses)], 'black', 'LineWidth', 5);
xlim([0 3])
ylim([-2 100])

[p] = signrank(mean_clusters_baselines, mean_clusters_responses)







