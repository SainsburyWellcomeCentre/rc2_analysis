function [beta_matrix, sig_vector, model_type_vector] = load_and_classify_tuning(data, clusters, trial_group_label, tuning_type, p_thresh)
% LOAD_AND_CLASSIFY_TUNING Load tuning data and classify by best model type
%
% Inputs:
%   data              - Data object with load_tuning_curves methods
%   clusters          - Array of cluster structures with .id field
%   trial_group_label - Label for trial group (e.g., 'RT')
%   tuning_type       - Type of tuning: 'velocity' or 'acceleration'
%   p_thresh          - Significance threshold (default: 0.05)
%
% Outputs:
%   beta_matrix        - n_clusters x max_params matrix of model parameters
%   sig_vector         - n_clusters logical vector indicating significance
%   model_type_vector  - n_clusters cell array with best model names
%
% Note: This function now uses ModelSelectionTuning which selects the best
% model (flat, linear, quadratic, cubic, gaussian, relu, sigmoid) based on BIC.

if nargin < 5
    p_thresh = 0.05;
end

n_clusters = length(clusters);
max_params = 4;  % Maximum number of parameters across all models
beta_matrix = nan(n_clusters, max_params);
sig_vector = false(n_clusters, 1);
model_type_vector = cell(n_clusters, 1);

for c = 1:n_clusters
    cluster_id = clusters(c).id;
    
    try
        % Load tuning data based on type
        if strcmpi(tuning_type, 'velocity')
            tuning = data.load_tuning_curves(cluster_id, trial_group_label);
        elseif strcmpi(tuning_type, 'acceleration')
            tuning_cell = data.load_tuning_curves_acceleration(cluster_id, trial_group_label);
            if iscell(tuning_cell) && length(tuning_cell) >= 1
                tuning = tuning_cell{1};  % Use "all" acceleration data
            else
                continue;
            end
        else
            error('Invalid tuning_type. Must be "velocity" or "acceleration".');
        end
        
        % Check if tuning data has new ModelSelectionTuning structure
        if ~isstruct(tuning) || ~isfield(tuning, 'shuffled') || ...
           ~isstruct(tuning.shuffled)
            continue;
        end
        
        % Extract model selection results
        if isfield(tuning.shuffled, 'best_model') && isfield(tuning.shuffled, 'best_model_info')
            model_name = tuning.shuffled.best_model;
            model_info = tuning.shuffled.best_model_info;
            
            % Store model type
            model_type_vector{c} = model_name;
            
            % Store beta coefficients (pad with NaN if fewer parameters)
            if isfield(model_info, 'beta')
                n_beta = length(model_info.beta);
                beta_matrix(c, 1:n_beta) = model_info.beta;
            end
            
            % Check significance
            if isfield(tuning.shuffled, 'p')
                sig_vector(c) = tuning.shuffled.p < p_thresh;
            end
        else
            % Fallback for old ShuffleTuning format (backward compatibility)
            if isfield(tuning.shuffled, 'beta') && length(tuning.shuffled.beta) >= 3
                beta_matrix(c, 1:3) = tuning.shuffled.beta;
                model_type_vector{c} = 'quadratic';  % Assume old format was quadratic
                
                if isfield(tuning.shuffled, 'p')
                    sig_vector(c) = tuning.shuffled.p < p_thresh;
                end
            end
        end
        
    catch ME
        % Silently continue on errors
        continue;
    end
end

end
