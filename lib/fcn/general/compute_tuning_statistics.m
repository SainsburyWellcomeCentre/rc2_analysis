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
%           .n_vel_linear       - Number of linear velocity fits
%           .n_vel_quadratic    - Number of quadratic velocity fits
%           .n_acc_linear       - Number of linear acceleration fits
%           .n_acc_quadratic    - Number of quadratic acceleration fits

if nargin < 4
    p_thresh = 0.05;
end

n_clusters = length(clusters);
stats.n_clusters = n_clusters;
stats.n_vel_tuned = 0;
stats.n_acc_tuned = 0;
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
                
                % Classify as linear or quadratic using F-test
                if isfield(tuning_vel, 'tuning') && isfield(tuning_vel, 'bin_centers') && ...
                   isfield(tuning_vel.shuffled, 'beta')
                    x_data = tuning_vel.bin_centers(:);
                    y_data = nanmean(tuning_vel.tuning, 2);
                    
                    is_linear = test_linear_vs_quadratic_fit(x_data, y_data, tuning_vel.shuffled.beta, p_thresh);
                    
                    if is_linear
                        stats.n_vel_linear = stats.n_vel_linear + 1;
                    else
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
                    
                    % Classify as linear or quadratic using F-test
                    if isfield(tuning_acc{1}, 'tuning') && isfield(tuning_acc{1}, 'bin_centers') && ...
                       isfield(tuning_acc{1}.shuffled, 'beta')
                        x_data = tuning_acc{1}.bin_centers(:);
                        y_data = nanmean(tuning_acc{1}.tuning, 2);
                        
                        is_linear = test_linear_vs_quadratic_fit(x_data, y_data, tuning_acc{1}.shuffled.beta, p_thresh);
                        
                        if is_linear
                            stats.n_acc_linear = stats.n_acc_linear + 1;
                        else
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
