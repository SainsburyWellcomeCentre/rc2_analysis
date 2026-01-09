function [beta_matrix, sig_vector, is_linear_vector] = load_and_classify_tuning(data, clusters, trial_group_label, tuning_type, p_thresh)
% LOAD_AND_CLASSIFY_TUNING Load tuning data and classify as linear or quadratic
%
% Inputs:
%   data              - Data object with load_tuning_curves methods
%   clusters          - Array of cluster structures with .id field
%   trial_group_label - Label for trial group (e.g., 'RT')
%   tuning_type       - Type of tuning: 'velocity' or 'acceleration'
%   p_thresh          - Significance threshold (default: 0.05)
%
% Outputs:
%   beta_matrix       - n_clusters x 3 matrix of beta coefficients [x², x¹, x⁰]
%   sig_vector        - n_clusters logical vector indicating significance
%   is_linear_vector  - n_clusters logical vector indicating linear fit is sufficient

if nargin < 5
    p_thresh = 0.05;
end

n_clusters = length(clusters);
beta_matrix = nan(n_clusters, 3);
sig_vector = false(n_clusters, 1);
is_linear_vector = false(n_clusters, 1);

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
        
        % Check if tuning data is valid
        if ~isstruct(tuning) || ~isfield(tuning, 'shuffled') || ...
           ~isstruct(tuning.shuffled) || ~isfield(tuning.shuffled, 'beta') || ...
           length(tuning.shuffled.beta) < 3
            continue;
        end
        
        % Store beta coefficients
        beta_matrix(c, :) = tuning.shuffled.beta;
        
        % Check significance
        if isfield(tuning.shuffled, 'p')
            sig_vector(c) = tuning.shuffled.p < p_thresh;
        end
        
        % Classify as linear or quadratic using F-test
        if sig_vector(c) && isfield(tuning, 'tuning') && isfield(tuning, 'bin_centers')
            x_data = tuning.bin_centers(:);
            y_data = nanmean(tuning.tuning, 2);  % Average firing rate per bin
            
            is_linear_vector(c) = test_linear_vs_quadratic_fit(x_data, y_data, tuning.shuffled.beta, p_thresh);
        end
        
    catch ME
        % Silently continue on errors (warnings already printed in main script)
        continue;
    end
end

end
