
close all
%%
experiment_groups       = {'darkness', 'mismatch_darkness_oct21'};
trial_group_labels      = {{'T_bank', 'T_RT', 'T_R'}, 'T'};
% experiment_groups = {'passive_same_luminance'};
% trial_group_labels = {'T_Vstatic'};

modalities = ["all", "acc", "dec"];
% figure_dir  = {'tuning_curves', 'passive_same_luminance_acc'};

max_number = 82; %151;
%%
ctl                     = RC2Analysis();
% ctl.setup_figures(figure_dir, save_figs);

classification      = {};
c                   = 0;
speed_tuned         = {};
for ll = 1 : length(trial_group_labels)
    tuning_acc      = {};
    tuning_vel      = {};
    p_svm           = [];
    direction       = [];
    probe_id        = {};
    cluster_id      = [];
    cluster_region  = {};

    
    probe_ids           = ctl.get_probe_ids(experiment_groups{ll});

    for ii = 1 : length(probe_ids)

        data        = ctl.load_formatted_data(probe_ids{ii});
        clusters    = data.VISp_clusters();

        for jj = 1 : length(clusters)
          [~, p_svm(ii, jj), direction(ii, jj)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{ll});
          tuning_acc{ii}{jj} = data.load_tuning_curves_acceleration(clusters(jj).id, trial_group_labels{ll});
          tuning_vel{ii}{jj} = data.load_tuning_curves(clusters(jj).id, trial_group_labels{ll});
        end
    end

        

    %% plot

    for ii = 1 : length(probe_ids)  
        for kk = 1 : length(tuning_acc{ii})
            c = c + 1;
            
            if p_svm(ii, kk) < 0.05 && direction(ii, kk) == 1
                direction
            elseif p_svm(ii, kk) < 0.05 && direction(ii, kk) == -1
                direction
            end
                
            % Acceleration tuning
            for acc_dec_i = 2 : 3

                this_tuning = tuning_acc{ii}{kk}{acc_dec_i};
                classification{acc_dec_i}{1}{max_number} = [];
                classification{acc_dec_i}{2}{max_number} = [];
                
                % exc
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)
                    classification{acc_dec_i}{1}{c} = 2; % high
                end
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)
                    classification{acc_dec_i}{1}{c} = 3; % low
                end
                
                % supp
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr)
                    classification{acc_dec_i}{2}{c} = 2;
                end
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr)
                    classification{acc_dec_i}{2}{c} = 3;
                end
            end
            
            % Velocity tuning
            this_tuning = tuning_vel{ii}{kk};
            speed_tuned{82} = [];
            if (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)) || ...
               (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)) || ...
               (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr)) || ...
               (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr))
               
                speed_tuned{c} = 1;
            end 
        end
    end
end

%% All acceleration

not_tuned_acc_exc = sum(cellfun('isempty', classification{2}{1}))
high_acc_exc = sum(cell2mat(classification{2}{1}) == 2)
low_acc_exc = sum(cell2mat(classification{2}{1}) == 3)


not_tuned_acc_supp = sum(cellfun('isempty', classification{2}{2}))
high_acc_supp = sum(cell2mat(classification{2}{2}) == 2)
low_acc_supp = sum(cell2mat(classification{2}{2}) == 3)


not_tuned_dec_exc = sum(cellfun('isempty', classification{3}{1}))
high_dec_exc = sum(cell2mat(classification{3}{1}) == 2)
low_dec_exc = sum(cell2mat(classification{3}{1}) == 3)


not_tuned_dec_supp = sum(cellfun('isempty', classification{3}{2}))
high_dec_supp = sum(cell2mat(classification{3}{2}) == 2)
low_dec_supp = sum(cell2mat(classification{3}{2}) == 3)


%% Not speed tuned
% {acc / dec} {exc / supp} {cluster}

tuned_to_speed =  sum(~cellfun('isempty', speed_tuned))
                     
tuned_to_acc_or_dec = sum((~cellfun('isempty', classification{2}{1}) | ~cellfun('isempty', classification{2}{2}) | ...
                           ~cellfun('isempty', classification{3}{1}) | ~cellfun('isempty', classification{3}{2})) & ...
                           cellfun('isempty', speed_tuned))

                       
tuned_exc_only_acc = sum((~cellfun('isempty', classification{2}{1}) | ~cellfun('isempty', classification{3}{1})) & cellfun('isempty', speed_tuned))
tuned_supp_only_acc = sum((~cellfun('isempty', classification{2}{2}) | ~cellfun('isempty', classification{3}{2})) & cellfun('isempty', speed_tuned))

tuned_only_acc = sum((~cellfun('isempty', classification{2}{1}) | ~cellfun('isempty', classification{2}{2})) & cellfun('isempty', speed_tuned))
tuned_only_dec = sum((~cellfun('isempty', classification{3}{1}) | ~cellfun('isempty', classification{3}{2})) & cellfun('isempty', speed_tuned))
tuned_only_acc_dec_both = sum(((~cellfun('isempty', classification{2}{1}) | ~cellfun('isempty', classification{2}{2})) & ...
                           (~cellfun('isempty', classification{3}{1}) | ~cellfun('isempty', classification{3}{2}))) & ...
                           cellfun('isempty', speed_tuned))

acc_exc_not_speed_tuned = [classification{2}{1}{cellfun('isempty', speed_tuned)}];
high_acc_exc = sum(acc_exc_not_speed_tuned == 2)
low_acc_exc = sum(acc_exc_not_speed_tuned == 3)

acc_supp_not_speed_tuned = [classification{2}{2}{cellfun('isempty', speed_tuned)}];
high_acc_supp = sum(acc_supp_not_speed_tuned == 2)
low_acc_supp = sum(acc_supp_not_speed_tuned == 3)

dec_exc_not_speed_tuned = [classification{3}{1}{cellfun('isempty', speed_tuned)}];
high_dec_exc = sum(dec_exc_not_speed_tuned == 2)
low_dec_exc = sum(dec_exc_not_speed_tuned == 3)

dec_supp_not_speed_tuned = [classification{3}{2}{cellfun('isempty', speed_tuned)}];
high_dec_supp = sum(dec_supp_not_speed_tuned == 2)
low_dec_supp = sum(dec_supp_not_speed_tuned == 3)




