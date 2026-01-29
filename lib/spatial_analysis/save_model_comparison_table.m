function save_model_comparison_table(tuning_curves, accel_tuning_curves, cluster_ids, output_dir, probe_id)
% SAVE_MODEL_COMPARISON_TABLE Save BIC values for all models to CSV files
%
% SYNTAX:
%   save_model_comparison_table(tuning_curves, accel_tuning_curves, cluster_ids, output_dir, probe_id)
%
% INPUTS:
%   tuning_curves       - Cell array of velocity tuning curves (one per cluster)
%   accel_tuning_curves - Cell array of acceleration tuning curves (one per cluster)
%   cluster_ids         - Array of cluster IDs
%   output_dir          - Directory where CSV files will be saved
%   probe_id            - String identifier for the probe
%
% OUTPUTS:
%   Creates two CSV files:
%     - <probe_id>_velocity_model_comparison.csv
%     - <probe_id>_acceleration_model_comparison.csv
%
% DESCRIPTION:
%   Extracts BIC values for all fitted models from tuning curve data and saves
%   them to CSV files for easy inspection. Each row corresponds to a cluster,
%   and columns show BIC values for each model type, plus the best model name.

    % Ensure output directory exists
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Save velocity model comparison
    if ~isempty(tuning_curves)
        fprintf('  Saving velocity model comparison to CSV...\n');
        velocity_csv = fullfile(output_dir, sprintf('%s_velocity_model_comparison.csv', probe_id));
        save_model_table(tuning_curves, cluster_ids, velocity_csv, 'Velocity');
    end
    
    % Save acceleration model comparison
    if ~isempty(accel_tuning_curves)
        fprintf('  Saving acceleration model comparison to CSV...\n');
        accel_csv = fullfile(output_dir, sprintf('%s_acceleration_model_comparison.csv', probe_id));
        save_model_table(accel_tuning_curves, cluster_ids, accel_csv, 'Acceleration');
    end
end


function save_model_table(tuning_curves, cluster_ids, filepath, curve_type)
% Helper function to save model comparison table for one type of tuning
    
    n_clusters = length(cluster_ids);
    
    % Initialize data structure
    model_names = {'linear', 'quadratic', 'cubic', 'gaussian', 'relu', 'sigmoid'};
    n_models = length(model_names);
    
    % Preallocate arrays
    bic_matrix = nan(n_clusters, n_models);
    loglik_matrix = nan(n_clusters, n_models);
    best_models = cell(n_clusters, 1);
    p_values = nan(n_clusters, 1);
    
    % Extract data from each cluster
    for c = 1:n_clusters
        if ~isempty(tuning_curves{c}) && isstruct(tuning_curves{c}) && ...
           isfield(tuning_curves{c}, 'shuffled') && isstruct(tuning_curves{c}.shuffled)
            
            shuffled = tuning_curves{c}.shuffled;
            
            % Get best model info
            if isfield(shuffled, 'best_model')
                best_models{c} = shuffled.best_model;
            end
            
            % Get p-value
            if isfield(shuffled, 'p')
                p_values(c) = shuffled.p;
            end
            
            % Extract BIC and log-likelihood for each model
            if isfield(shuffled, 'all_models') && isstruct(shuffled.all_models)
                all_models = shuffled.all_models;
                
                for m = 1:n_models
                    model = model_names{m};
                    if isfield(all_models, model) && isstruct(all_models.(model))
                        if isfield(all_models.(model), 'bic')
                            bic_matrix(c, m) = all_models.(model).bic;
                        end
                        if isfield(all_models.(model), 'log_likelihood')
                            loglik_matrix(c, m) = all_models.(model).log_likelihood;
                        end
                    end
                end
            end
        else
            best_models{c} = 'N/A';
        end
    end
    
    % Create table
    T = table(cluster_ids(:), ...
              bic_matrix(:,1), bic_matrix(:,2), bic_matrix(:,3), ...
              bic_matrix(:,4), bic_matrix(:,5), bic_matrix(:,6), ...
              loglik_matrix(:,1), loglik_matrix(:,2), loglik_matrix(:,3), ...
              loglik_matrix(:,4), loglik_matrix(:,5), loglik_matrix(:,6), ...
              best_models, p_values, ...
              'VariableNames', {'ClusterID', ...
                                'BIC_Linear', 'BIC_Quadratic', 'BIC_Cubic', ...
                                'BIC_Gaussian', 'BIC_ReLU', 'BIC_Sigmoid', ...
                                'LogLik_Linear', 'LogLik_Quadratic', 'LogLik_Cubic', ...
                                'LogLik_Gaussian', 'LogLik_ReLU', 'LogLik_Sigmoid', ...
                                'BestModel', 'PValue'});
    
    % Write to CSV
    writetable(T, filepath);
    fprintf('    Saved %s model comparison for %d clusters to: %s\n', curve_type, n_clusters, filepath);
end
