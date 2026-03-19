function save_tuning_summary(tuning_curves, accel_tuning_curves, ...
                              tuning_curves_long, tuning_curves_short, ...
                              accel_tuning_curves_long, accel_tuning_curves_short, ...
                              cluster_ids, output_dir, probe_id)
% SAVE_TUNING_SUMMARY Save a text summary of velocity and acceleration tuning results
%
% INPUTS:
%   tuning_curves              - Combined velocity tuning curves
%   accel_tuning_curves        - Combined acceleration tuning curves
%   tuning_curves_long         - Long trial velocity tuning curves
%   tuning_curves_short        - Short trial velocity tuning curves
%   accel_tuning_curves_long   - Long trial acceleration tuning curves
%   accel_tuning_curves_short  - Short trial acceleration tuning curves
%   cluster_ids                - Array of cluster IDs
%   output_dir                 - Directory to save the text file
%   probe_id                   - Probe identifier string

    % Define output file path
    fname = fullfile(output_dir, 'tuning_summary_all_probes.txt');
    
    % Check if file already exists BEFORE opening
    file_exists = exist(fname, 'file');
    
    % Open in append mode
    fid = fopen(fname, 'a');
    
    if fid == -1
        warning('Could not open file for writing: %s', fname);
        return;
    end
    
    % Only write global header if file is new
    if ~file_exists
        fprintf(fid, '========================================================\n');
        fprintf(fid, 'TUNING SUMMARY - ALL PROBES\n');
        fprintf(fid, '========================================================\n\n');
    end
    
    % Always write probe header
    fprintf(fid, '\n========================================================\n');
    fprintf(fid, 'PROBE: %s\n', probe_id);
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, '========================================================\n\n');

    % Write column headers
    fprintf(fid, '%-10s | %-20s | %-20s | %-20s | %-20s | %-20s | %-20s\n', ...
        'Cluster', ...
        'Vel (combined)', ...
        'Vel (long)', ...
        'Vel (short)', ...
        'Accel (combined)', ...
        'Accel (long)', ...
        'Accel (short)');
    fprintf(fid, '%s\n', repmat('-', 1, 130));
    
    % Write each column description
    fprintf(fid, '%-10s | %-20s | %-20s | %-20s | %-20s | %-20s | %-20s\n', ...
        '', ...
        'model/p/r', ...
        'model/p/r', ...
        'model/p/r', ...
        'model/p/r', ...
        'model/p/r', ...
        'model/p/r');
    fprintf(fid, '%s\n\n', repmat('-', 1, 130));
    
    % Write results for each cluster
    for c = 1:length(cluster_ids)
        cluster_id = cluster_ids(c);
        
        % Extract info for each condition
        vel_combined  = extract_tuning_info(tuning_curves{c});
        vel_long      = extract_tuning_info(tuning_curves_long{c});
        vel_short     = extract_tuning_info(tuning_curves_short{c});
        accel_combined = extract_tuning_info(accel_tuning_curves{c});
        accel_long    = extract_tuning_info(accel_tuning_curves_long{c});
        accel_short   = extract_tuning_info(accel_tuning_curves_short{c});
        
        % Write cluster header
        fprintf(fid, '\nCluster %d\n', cluster_id);
        fprintf(fid, '%s\n', repmat('-', 1, 130));
        
        % Write tuning status line
        fprintf(fid, '  Velocity tuning:\n');
        fprintf(fid, '    Combined : %s\n', format_tuning_info(vel_combined));
        fprintf(fid, '    Long     : %s\n', format_tuning_info(vel_long));
        fprintf(fid, '    Short    : %s\n', format_tuning_info(vel_short));
        
        fprintf(fid, '  Acceleration tuning:\n');
        fprintf(fid, '    Combined : %s\n', format_tuning_info(accel_combined));
        fprintf(fid, '    Long     : %s\n', format_tuning_info(accel_long));
        fprintf(fid, '    Short    : %s\n', format_tuning_info(accel_short));
        
        % Write overall tuning summary
        fprintf(fid, '  Summary:\n');
        fprintf(fid, '    Velocity tuned (any condition)      : %s\n', ...
            any_tuned(vel_combined, vel_long, vel_short));
        fprintf(fid, '    Velocity tuned (persists long+short): %s\n', ...
            both_tuned(vel_long, vel_short));
        fprintf(fid, '    Accel tuned (any condition)         : %s\n', ...
            any_tuned(accel_combined, accel_long, accel_short));
        fprintf(fid, '    Accel tuned (persists long+short)   : %s\n', ...
            both_tuned(accel_long, accel_short));
        fprintf(fid, '\n');
    end
    
    % Write overall statistics
    fprintf(fid, '========================================================\n');
    fprintf(fid, 'OVERALL STATISTICS\n');
    fprintf(fid, '========================================================\n');
    
    n_clusters = length(cluster_ids);
    
    % Count tuned clusters
    n_vel_tuned_combined = sum(cellfun(@(x) is_tuned(x), tuning_curves));
    n_vel_tuned_long = sum(cellfun(@(x) is_tuned(x), tuning_curves_long));
    n_vel_tuned_short = sum(cellfun(@(x) is_tuned(x), tuning_curves_short));
    n_accel_tuned_combined = sum(cellfun(@(x) is_tuned(x), accel_tuning_curves));
    n_accel_tuned_long = sum(cellfun(@(x) is_tuned(x), accel_tuning_curves_long));
    n_accel_tuned_short = sum(cellfun(@(x) is_tuned(x), accel_tuning_curves_short));
    
    fprintf(fid, '\nVelocity tuning:\n');
    fprintf(fid, '  Combined : %d/%d (%.1f%%)\n', n_vel_tuned_combined, n_clusters, 100*n_vel_tuned_combined/n_clusters);
    fprintf(fid, '  Long     : %d/%d (%.1f%%)\n', n_vel_tuned_long, n_clusters, 100*n_vel_tuned_long/n_clusters);
    fprintf(fid, '  Short    : %d/%d (%.1f%%)\n', n_vel_tuned_short, n_clusters, 100*n_vel_tuned_short/n_clusters);
    
    fprintf(fid, '\nAcceleration tuning:\n');
    fprintf(fid, '  Combined : %d/%d (%.1f%%)\n', n_accel_tuned_combined, n_clusters, 100*n_accel_tuned_combined/n_clusters);
    fprintf(fid, '  Long     : %d/%d (%.1f%%)\n', n_accel_tuned_long, n_clusters, 100*n_accel_tuned_long/n_clusters);
    fprintf(fid, '  Short    : %d/%d (%.1f%%)\n', n_accel_tuned_short, n_clusters, 100*n_accel_tuned_short/n_clusters);
    
    fclose(fid);
    fprintf('  Saved tuning summary to: %s\n', fname);
end


%% Helper functions

function info = extract_tuning_info(tc)
% Extract tuning info from a tuning curve struct
    info.is_tuned = false;
    info.model = 'N/A';
    info.p = NaN;
    info.r = NaN;
    info.bic = NaN;
    
    if isempty(tc) || ~isstruct(tc)
        return;
    end
    
    if ~isfield(tc, 'shuffled') || isempty(tc.shuffled)
        return;
    end
    
    % New ModelSelectionTuning format
    if isfield(tc.shuffled, 'best_model')
        info.model = tc.shuffled.best_model;
        info.p = tc.shuffled.p;
        info.r = tc.shuffled.best_model_info.fit_metric;
        info.bic = tc.shuffled.best_model_info.bic;
        info.is_tuned = (tc.shuffled.p < 0.05);
        
    % Old ShuffleTuning format
    elseif isfield(tc.shuffled, 'r')
        info.model = 'linear';
        info.p = tc.shuffled.p;
        info.r = tc.shuffled.r;
        info.is_tuned = (tc.shuffled.p < 0.05);
    end
end


function str = format_tuning_info(info)
% Format tuning info as a readable string
    if isnan(info.p)
        str = 'No data';
    elseif info.is_tuned
        str = sprintf('TUNED | model=%s | p=%.3f | r=%.3f | BIC=%.1f', ...
            upper(info.model), info.p, info.r, info.bic);
    else
        str = sprintf('NOT tuned | model=%s | p=%.3f | r=%.3f | BIC=%.1f', ...
            upper(info.model), info.p, info.r, info.bic);
    end
end


function result = any_tuned(combined, long, short)
% Check if tuned in any condition
    if combined.is_tuned || long.is_tuned || short.is_tuned
        result = 'YES';
    else
        result = 'NO';
    end
end


function result = both_tuned(long, short)
% Check if tuning persists in both long and short
    if long.is_tuned && short.is_tuned
        result = 'YES (both)';
    elseif long.is_tuned
        result = 'PARTIAL (long only)';
    elseif short.is_tuned
        result = 'PARTIAL (short only)';
    else
        result = 'NO';
    end
end


function result = is_tuned(tc)
% Return true if tuning curve is significant
    result = false;
    if isempty(tc) || ~isstruct(tc)
        return;
    end
    if isfield(tc, 'shuffled') && isfield(tc.shuffled, 'p')
        result = (tc.shuffled.p < 0.05);
    end
end