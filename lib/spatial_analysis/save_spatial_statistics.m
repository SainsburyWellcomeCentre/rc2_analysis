function save_spatial_statistics(spatial_tuning_stats, probe_id, output_dir, analysis_params, cluster_ids, group_names, group_labels)
% SAVE_SPATIAL_STATISTICS Save statistical results and generate summary reports
%
%   save_spatial_statistics(spatial_tuning_stats, probe_id, output_dir, analysis_params, ...
%                          cluster_ids, group_names, group_labels)
%
%   Saves spatial firing rate statistics to a text summary file for a probe.
%   Creates summary statistics and quick reference table for all clusters.
%
%   Inputs:
%       spatial_tuning_stats  - Struct with statistical test results per cluster and group
%       probe_id          - String identifier for the probe
%       output_dir        - Directory path where files will be saved
%       analysis_params   - Struct containing analysis parameters:
%                             .bin_size_cm, .gauss_sigma_cm, .minSpikes,
%                             .minPeakRate, .fieldFrac, .minFieldBins,
%                             .maxNumFields, .nShuf, .pThresh
%       cluster_ids       - Array of cluster IDs
%       group_names       - Cell array of group names (e.g., {'long', 'short'})
%       group_labels      - Cell array of descriptive labels for groups
%
%   Outputs:
%       Creates one file in output_dir:
%         <probe_id>_statistics_summary.txt - Human-readable summary report

    if isempty(spatial_tuning_stats)
        fprintf('  No statistics to save (spatial_tuning_stats is empty)\n');
        return;
    end

    % Create text summary file
    stats_txt_filename = sprintf('%s_statistics_summary.txt', probe_id);
    stats_txt_filepath = fullfile(output_dir, stats_txt_filename);
    fprintf('  Saving text summary to: %s\n', stats_txt_filename);
    
    fid = fopen(stats_txt_filepath, 'w');
    
    % Write header
    write_header(fid, probe_id);
    
    % Write analysis parameters
    write_parameters(fid, analysis_params);
    
    % Write summary statistics
    [n_pass_long, n_pass_short, n_pass_both] = write_summary_stats(fid, spatial_tuning_stats, group_names);
    
    % Write quick reference table (detailed results removed per user request)
    write_quick_reference_table(fid, spatial_tuning_stats, cluster_ids);
    
    fclose(fid);
end

%% Helper functions

function write_header(fid, probe_id)
    fprintf(fid, '=======================================================\n');
    fprintf(fid, 'SPATIAL FIRING RATE STATISTICS SUMMARY\n');
    fprintf(fid, '=======================================================\n');
    fprintf(fid, 'Probe: %s\n', probe_id);
    fprintf(fid, 'Analysis date: %s\n', datestr(now));
    fprintf(fid, '\n');
end

function write_parameters(fid, params)
    fprintf(fid, '--- Analysis Parameters ---\n');
    fprintf(fid, 'Bin size: %.1f cm\n', params.bin_size_cm);
    fprintf(fid, 'Gaussian sigma: %.1f cm\n', params.gauss_sigma_cm);
    fprintf(fid, 'Min spikes: %d\n', params.minSpikes);
    fprintf(fid, 'Min peak rate: %.1f Hz\n', params.minPeakRate);
    fprintf(fid, 'Field fraction: %.2f\n', params.fieldFrac);
    fprintf(fid, 'Min field bins: %d (%.1f cm)\n', params.minFieldBins, params.minFieldBins * params.bin_size_cm);
    fprintf(fid, 'Max num fields: %d\n', params.maxNumFields);
    fprintf(fid, 'Num shuffles: %d\n', params.nShuf);
    fprintf(fid, 'P-value threshold: %.3f\n', params.pThresh);
    fprintf(fid, '\n');
end

function [n_pass_long, n_pass_short, n_pass_both] = write_summary_stats(fid, spatial_tuning_stats, group_names)
    % Compute summary counts
    cluster_names = fieldnames(spatial_tuning_stats);
    n_total = length(cluster_names);
    n_pass_long = 0;
    n_pass_short = 0;
    n_pass_both = 0;
    
    for c_idx = 1:length(cluster_names)
        cluster_name = cluster_names{c_idx};
        stats = spatial_tuning_stats.(cluster_name);
        pass_long = isfield(stats, 'long') && stats.long.isSpatiallyTuned;
        pass_short = isfield(stats, 'short') && stats.short.isSpatiallyTuned;
        if pass_long; n_pass_long = n_pass_long + 1; end
        if pass_short; n_pass_short = n_pass_short + 1; end
        if pass_long && pass_short; n_pass_both = n_pass_both + 1; end
    end
    
    fprintf(fid, '--- Summary Statistics ---\n');
    fprintf(fid, 'Total clusters analyzed: %d\n', n_total);
    fprintf(fid, 'Clusters with spatial tuning (long trials): %d (%.1f%%)\n', n_pass_long, 100*n_pass_long/n_total);
    fprintf(fid, 'Clusters with spatial tuning (short trials): %d (%.1f%%)\n', n_pass_short, 100*n_pass_short/n_total);
    fprintf(fid, 'Clusters with spatial tuning in both conditions: %d (%.1f%%)\n', n_pass_both, 100*n_pass_both/n_total);
    fprintf(fid, '\n');
end

function write_detailed_results(fid, spatial_tuning_stats, cluster_ids, group_names, group_labels, params)
    fprintf(fid, '=======================================================\n');
    fprintf(fid, 'DETAILED RESULTS BY CLUSTER\n');
    fprintf(fid, '=======================================================\n\n');
    
    cluster_names = fieldnames(spatial_tuning_stats);
    
    for c_idx = 1:length(cluster_names)
        cluster_name = cluster_names{c_idx};
        cluster_id = str2double(strrep(cluster_name, 'cluster_', ''));
        stats = spatial_tuning_stats.(cluster_name);
        
        fprintf(fid, '--- Cluster %d ---\n', cluster_id);
        
        % Process each group
        group_fields = fieldnames(stats);
        for g_idx = 1:length(group_fields)
            group = group_fields{g_idx};
            s = stats.(group);
            
            % Find corresponding group label
            group_label = group;
            for gl_idx = 1:length(group_names)
                if strcmp(group_names{gl_idx}, group)
                    group_label = group_labels{gl_idx};
                    break;
                end
            end
            
            fprintf(fid, '\n  %s:\n', group_label);
            fprintf(fid, '    Number of spikes: %d\n', s.n_spikes);
            
            % Test results
            fprintf(fid, '    \n');
            fprintf(fid, '    Test Results:\n');
            fprintf(fid, '      [%s] Spike count (>= %d): %d spikes\n', ...
                    test_status_str(s.n_spikes >= params.minSpikes), params.minSpikes, s.n_spikes);
            fprintf(fid, '      [%s] Peak rate (>= %.1f Hz): %.2f Hz at %.1f cm\n', ...
                    test_status_str(s.peakRate >= params.minPeakRate), params.minPeakRate, s.peakRate, s.peakPosition);
            fprintf(fid, '      [%s] Information content (p < %.3f): info=%.4f, p=%.4f\n', ...
                    test_status_str(s.pVal < params.pThresh), params.pThresh, s.info, s.pVal);
            fprintf(fid, '      [%s] Field size (>= %d bins = %.1f cm): %d bins = %.1f cm\n', ...
                    test_status_str(s.fieldBins >= params.minFieldBins), params.minFieldBins, ...
                    params.minFieldBins * params.bin_size_cm, s.fieldBins, s.fieldBins * params.bin_size_cm);
            fprintf(fid, '      [%s] Field uniqueness (<= %d fields): %d fields\n', ...
                    test_status_str(s.numFields <= params.maxNumFields), params.maxNumFields, s.numFields);
            fprintf(fid, '    \n');
            
            % Final result
            if s.isSpatiallyTuned
                fprintf(fid, '    >>> FINAL RESULT: PASS (spatially tuned)\n');
            else
                fprintf(fid, '    >>> FINAL RESULT: FAIL (%s)\n', s.reason);
            end
        end
        fprintf(fid, '\n');
    end
end

function write_quick_reference_table(fid, spatial_tuning_stats, cluster_ids)
    fprintf(fid, '=======================================================\n');
    fprintf(fid, 'QUICK REFERENCE TABLE\n');
    fprintf(fid, '=======================================================\n');
    fprintf(fid, '%-12s | %-15s | %-15s\n', 'Cluster ID', 'Long Trials', 'Short Trials');
    fprintf(fid, '-------------------------------------------------------\n');
    
    cluster_names = fieldnames(spatial_tuning_stats);
    
    for c_idx = 1:length(cluster_names)
        cluster_name = cluster_names{c_idx};
        cluster_id = str2double(strrep(cluster_name, 'cluster_', ''));
        stats = spatial_tuning_stats.(cluster_name);
        
        long_result = 'N/A';
        short_result = 'N/A';
        
        if isfield(stats, 'long')
            if stats.long.isSpatiallyTuned
                long_result = 'PASS';
            else
                long_result = sprintf('FAIL (%s)', stats.long.reason);
            end
        end
        
        if isfield(stats, 'short')
            if stats.short.isSpatiallyTuned
                short_result = 'PASS';
            else
                short_result = sprintf('FAIL (%s)', stats.short.reason);
            end
        end
        
        fprintf(fid, '%-12d | %-15s | %-15s\n', cluster_id, long_result, short_result);
    end
end

function str = test_status_str(passed)
    % Return a string indicating pass/fail status
    if passed
        str = 'PASS';
    else
        str = 'FAIL';
    end
end
