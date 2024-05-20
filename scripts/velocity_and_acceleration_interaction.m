%  First, for each probe, understand if certain values of speed are
%  associated with significantly different accelerations
close all 
clear all

experiment_groups       = {'darkness'}; %{'darkness', 'mismatch_darkness_oct21'};
trial_group_labels      = {'T_bank', 'T_RT', 'T_R'}; %{{'T_bank', 'T_RT', 'T_R'}, 'T'};

ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

cluster_mean = [];
for probe_i = 1 : length(probe_ids)
    data        = ctl.load_formatted_data(probe_ids{probe_i});
    
    trials      = data.motion_trials();
    clusters    = data.VISp_clusters();
    
%     figure(probe_i)
%     hold on
    
    for cluster_i = 1: length(clusters)
        tuning_vel = data.load_tuning_curves(clusters(cluster_i).id, trial_group_labels);
        tuning_acc = data.load_tuning_curves_acceleration(clusters(cluster_i).id, trial_group_labels);
        tuning_acc = tuning_acc{1, 1}; % Both acceleration and deceleration
        
        bin_edges_vel = tuning_vel.bin_edges;
        bin_edges_acc = tuning_acc.bin_edges;
        
        fr_per_bin = zeros(length(trials), length(bin_edges_vel) - 1, length(bin_edges_acc) -1);
        
        for trial_i = 1 : length(trials)
            trial = trials{trial_i}.to_aligned;
        
        % See interaction of velocity and acceleration
%         rounded_velocity = round(trial.velocity, 2);
%         rounded_acceleration = round(trial.acceleration, 2);
%         scatter(rounded_velocity, rounded_acceleration, 0.1);
           
            mmask = trials{trial_i}.motion_mask();

            acc_bin_mask = false(length(trial.acceleration), length(bin_edges_acc) - 1);
            for bin_i = 1 : (length(bin_edges_acc) - 1)   

               mask = trial.acceleration >= bin_edges_acc(bin_i) & ...
                    trial.acceleration < bin_edges_acc(bin_i+1);

                acc_bin_mask(:, bin_i) = mask & mmask(1:length(mask));
            end

            vel_bin_mask = false(length(trial.velocity), length(bin_edges_vel) - 1);
            for bin_i = 1 : (length(bin_edges_vel) - 1)   

                mask = trial.velocity >= bin_edges_vel(bin_i) & ...
                    trial.velocity < bin_edges_vel(bin_i+1);

                vel_bin_mask(:, bin_i) = mask & mmask(1:length(mask));
            end
            
            
            convolved_fr = clusters(cluster_i).fr.get_convolution(trial.probe_t);
           
            for bin_vel = 1 : (length(bin_edges_vel) - 1)   
                for bin_acc = 1 : (length(bin_edges_acc) - 1)   
                    this_acc_this_vel_mask = vel_bin_mask(:, bin_vel) & acc_bin_mask(:, bin_acc);
                    fr_per_bin(trial_i, bin_vel, bin_acc) = nanmean(convolved_fr(this_acc_this_vel_mask));
                end
            end
         end
%         cluster_mean(end+1) = squeeze(nanmean(fr_per_bin));
        figure(cluster_i);
        
        mean_response = squeeze(nanmean(fr_per_bin));
        
        subplot(2,2,1);
        % This mean looks different from the one of the tuning because it's
        % calculated across bins and not across trials...
        area(gca, flip(nanmean(mean_response, 2)))
%         area(gca, nanmean(tuning_acc.tuning, 2))
        xlim([1 18])
        xlabel('Accelereation (m/s^2)')
        ylabel('Firing rate (Hz)')
        set(gca,'view',[-90 90])
        len = (length(bin_edges_acc) - 1);
        xlim([1 len])
        xticks([1 : len])
        xticklabels(round(tuning_acc.bin_centers, 1))
        
        subplot(2,2,2);
        heatmap(mean_response);
        
        subplot(2,2,4);
        area(gca, nanmean(mean_response, 1))
         xlabel('Velocity (m/s)')
        ylabel('Firing rate (Hz)')
        len = (length(bin_edges_vel) - 1);
        xlim([1 len])
        xticks([1 : len])
        xticklabels(round(tuning_vel.bin_centers, 1))
        
    end
end