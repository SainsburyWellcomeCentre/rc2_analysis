function stats = compute_tuning_statistics(data, clusters, trial_group_label, p_thresh)
% COMPUTE_TUNING_STATISTICS Compute statistics for velocity and acceleration tuning
%
% Inputs:
%   data              - Data object with load_tuning_curves methods
%   clusters          - Array of cluster structures with .id field
%   trial_group_label - Label for trial group (e.g., 'RT')
%   p_thresh          - Significance threshold (default: 0.05)
%
% Outputs:
%   stats - Structure with fields:
%           .n_clusters         - Total number of clusters
%           .n_vel_tuned        - Number of velocity-tuned clusters
%           .n_acc_tuned        - Number of acceleration-tuned clusters
%           .n_vel_by_model     - Struct with counts per model type (flat, linear, etc.)
%           .n_acc_by_model     - Struct with counts per model type
%
% Note: Updated to work with ModelSelectionTuning which selects best model using BIC

if nargin < 4
    p_thresh = 0.05;
end

n_clusters = length(clusters);
stats.n_clusters = n_clusters;
stats.n_vel_tuned = 0;
stats.n_acc_tuned = 0;

% Initialize model type counters
model_types = {'flat', 'linear', 'quadratic', 'cubic', 'gaussian', 'relu', 'sigmoid'};
stats.n_vel_by_model = struct();
stats.n_acc_by_model = struct();
for i = 1:length(model_types)
    stats.n_vel_by_model.(model_types{i}) = 0;
    stats.n_acc_by_model.(model_types{i}) = 0;
end

% Legacy counts for backward compatibility
stats.n_vel_linear = 0;
stats.n_vel_quadratic = 0;
stats.n_acc_linear = 0;
stats.n_acc_quadratic = 0;

for c = 1:n_clusters
    cluster_id = clusters(c).id;
    
    % Velocity tuning
    try
        tuning_vel = data.load_tuning_curves(cluster_id, trial_group_label);
        if isstruct(tuning_vel) && isfield(tuning_vel, 'shuffled') && ...
           isstruct(tuning_vel.shuffled) && isfield(tuning_vel.shuffled, 'p')
            if tuning_vel.shuffled.p < p_thresh
                stats.n_vel_tuned = stats.n_vel_tuned + 1;
                
                % Get model type from ModelSelectionTuning
                if isfield(tuning_vel.shuffled, 'best_model')
                    model_name = tuning_vel.shuffled.best_model;
                    if isfield(stats.n_vel_by_model, model_name)
                        stats.n_vel_by_model.(model_name) = stats.n_vel_by_model.(model_name) + 1;
                    end
                    
                    % Update legacy counters
                    if strcmp(model_name, 'linear')
                        stats.n_vel_linear = stats.n_vel_linear + 1;
                    elseif strcmp(model_name, 'quadratic')
                        stats.n_vel_quadratic = stats.n_vel_quadratic + 1;
                    end
                end
            end
        end
    catch
        % Silently continue on errors
    end
    
    % Acceleration tuning
    try
        tuning_acc = data.load_tuning_curves_acceleration(cluster_id, trial_group_label);
        if iscell(tuning_acc) && length(tuning_acc) >= 1
            if isstruct(tuning_acc{1}) && isfield(tuning_acc{1}, 'shuffled') && ...
               isstruct(tuning_acc{1}.shuffled) && isfield(tuning_acc{1}.shuffled, 'p')
                if tuning_acc{1}.shuffled.p < p_thresh
                    stats.n_acc_tuned = stats.n_acc_tuned + 1;
                    
                    % Get model type from ModelSelectionTuning
                    if isfield(tuning_acc{1}.shuffled, 'best_model')
                        model_name = tuning_acc{1}.shuffled.best_model;
                        if isfield(stats.n_acc_by_model, model_name)
                            stats.n_acc_by_model.(model_name) = stats.n_acc_by_model.(model_name) + 1;
                        end
                        
                        % Update legacy counters
                        if strcmp(model_name, 'linear')
                            stats.n_acc_linear = stats.n_acc_linear + 1;
                        elseif strcmp(model_name, 'quadratic')
                            stats.n_acc_quadratic = stats.n_acc_quadratic + 1;
                        end
                    end
                end
            end
        end
    catch
        % Silently continue on errors
    end
end

end
