% Count how many clusters in the darkness datasets are tuned for speed
% and exclude speed tuned clusters
% Currently not counting here how many are tonically responsive 

% Tuning filters explanation
% tuning.shuffled.p = p-value calculated by ShuffleTuning
%   comparing the shuffled distribution to the cluster tuning
% tuning.shuffled.beta = polyfit result calculated by ShuffleTuning
%   it is the slope of the line fitted to the cluster tuning  
% filter: mean(nanmean(tuning.tuning)) > nanmean(tuning.stationary_fr)
%   compare mean response in motion periods across bins to the mean stationary period



% Initialize parameters
experiment_groups       = {'darkness', 'mismatch_darkness_oct21'};
trial_group_labels      = {{'T_bank', 'T_RT', 'T_R'}, 'T'};
modalities              = ["all", "acc", "dec"];


% Instantiate controller
ctl                 = RC2Analysis();

% initialize cell arrays
tuning_acc      = {};
tuning_vel      = {};

% initialize the table in which to save probe id, cluster id and all tuning
% iformation
columnNames = {'probe_id', 'cluster_id', 'tonic', 'speed_tuned', 'acc_exc_H', 'acc_exc_L', 'acc_supp_H', 'acc_supp_L', ...
                                                                 'dec_exc_H', 'dec_exc_L', 'dec_supp_H', 'dec_supp_L'};
variableTypes = {'string', 'int32', 'logical', 'logical', 'logical', 'logical', ...
                 'logical', 'logical', 'logical', 'logical', 'logical', 'logical'};

% Create the table with the specified sizes and types
T = table('Size', [0, length(columnNames)], 'VariableTypes', variableTypes, 'VariableNames', columnNames);

% loop through experimental groups
for ll = 1 : length(experiment_groups)

    % get probe ids
    probe_ids           = ctl.get_probe_ids(experiment_groups{ll});
    
    % loop across probe_ids
    for ii = 1 : length(probe_ids)
        % load data for this probe
        data        = ctl.load_formatted_data(probe_ids{ii});
        clusters    = data.VISp_clusters();
        
        % loop across clusters
        for jj = 1 : length(clusters)
            [~, p_svm, direction] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{ll});
            
            % load and store tuning curves    
            tuning_acc = data.load_tuning_curves_acceleration(clusters(jj).id, trial_group_labels{ll});
            tuning_vel = data.load_tuning_curves(clusters(jj).id, trial_group_labels{ll});
            
            % Tonic 
            if p_svm < 0.05 && direction ~= 0
                tonic_ = true;
            else
                tonic_ = false;
            end
            
            % Initialize all acceleration tunings to false
            acc_exc_H  = false;
            acc_exc_L  = false;
            acc_supp_H = false;
            acc_supp_L = false;
            dec_exc_H  = false;
            dec_exc_L  = false;
            dec_supp_H = false;
            dec_supp_L = false;
            
            % Filter for all acceleration tuning combinations
            for acc_dec_i = 2 : 3
                % copy tuning data to a variable to make filtering easier to read
                this_tuning = tuning_acc{acc_dec_i};
                
                % filter for excited clusters 
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)
                    % high excited (positive trend)
                    if acc_dec_i == 2 acc_exc_H = true, else dec_exc_H = true, end;
                end
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)
                    % low excited (negative trend)
                    if acc_dec_i == 2 acc_exc_L = true, else dec_exc_L = true, end;
                end
                
                % filter for suppressed clusters
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr)
                    % high suppressed (positive trend)
                    if acc_dec_i == 2 acc_supp_H = true, else dec_supp_H = true, end;
                end
                if this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr)
                    % low suppressed (negative trend)
                    if acc_dec_i == 2 acc_supp_L = true, else dec_supp_L = true, end;
                end
            end
            
            % Velocity tuning
            % copy tuning data to a variable to make filtering easier to read
            this_tuning = tuning_vel;
            
            % Save tuning for velocity in any combination
            if (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)) || ...
               (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) > nanmean(this_tuning.stationary_fr)) || ...
               (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) >= 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr)) || ...
               (this_tuning.shuffled.p < 0.05 && this_tuning.shuffled.beta(1) < 0 && mean(nanmean(this_tuning.tuning)) < nanmean(this_tuning.stationary_fr))
               
                speed_tuned_ = true;
            else 
                speed_tuned_ = false;
            end 
            
            % Save tuning information in a new row and add it to the table
            newRow = {probe_ids{ii}, clusters(jj).id, tonic_, speed_tuned_, ...
                                     acc_exc_H, acc_exc_L, acc_supp_H, acc_supp_L, ...
                                     dec_exc_H, dec_exc_L, dec_supp_H, dec_supp_L};
                                 
            T = [T; cell2table(newRow, 'VariableNames', columnNames)];
        end
    end
end


% ACCELERATION TUNED

% Extract the boolean columns
booleanColumns = T{:, 3:end};

% Sum the true values (ones) for each boolean column
trueCounts = sum(booleanColumns);

% Display the results
for i = 1:length(trueCounts)
    fprintf('Column %s has %d true values.\n', columnNames{i+2}, trueCounts(i));
end


% ACCELERATION TUNED BUT NOT SPEED TUNED
% Filter out rows where 'speed_tuned' is true
filteredT = T(~T.speed_tuned, :);

% Define the relevant column names
relevantColumns = {'acc_exc_H', 'acc_exc_L', 'acc_supp_H', 'acc_supp_L', ...
                   'dec_exc_H', 'dec_exc_L', 'dec_supp_H', 'dec_supp_L'};

% Iterate through each relevant column and count the true values
for i = 1:length(relevantColumns)
    columnName = relevantColumns{i};
    trueCount = sum(filteredT{:, columnName});
    fprintf('Column %s has %d true values (excluding speed_tuned).\n', columnName, trueCount);
end




