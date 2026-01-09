function analyzer = manage_spatial_analysis_cache(analyzer, cache_filepath, re_run_analysis)
% MANAGE_SPATIAL_ANALYSIS_CACHE Load or compute spatial analysis with caching
%
% SYNTAX:
%   analyzer = manage_spatial_analysis_cache(analyzer, cache_filepath, re_run_analysis)
%
% INPUTS:
%   analyzer         - SpatialTuningAnalyzer instance
%   cache_filepath   - Full path to cache file
%   re_run_analysis  - Boolean; true to force recomputation, false to load cache
%
% OUTPUTS:
%   analyzer         - Updated analyzer with loaded or computed data
%
% DESCRIPTION:
%   Manages loading/saving of spatial analysis cache. If cache exists and
%   re_run_analysis is false, loads cached data. Otherwise, runs full analysis
%   and saves to cache.

    % Try to load cached data if not re-running analysis
    if ~re_run_analysis && exist(cache_filepath, 'file')
        [~, cache_filename] = fileparts(cache_filepath);
        fprintf('  Loading cached analysis data from: %s.mat\n', cache_filename);
        analyzer.load_cache(cache_filepath);
        fprintf('  Loaded cached data for %d clusters\n', length(analyzer.cluster_ids));
        
        % Still need trial_groups for plotting
        analyzer.collect_trials();
        
        % Bin config is reconstructed in load_cache, but verify it exists
        if isempty(analyzer.bin_config)
            analyzer.setup_bins();
        end
        
        % Track if we need to update the cache
        cache_needs_update = false;
        
        % Check if TTG analysis was loaded; if not, compute it
        if isempty(analyzer.ttg_analysis) || ~isstruct(analyzer.ttg_analysis) || isempty(fieldnames(analyzer.ttg_analysis))
            fprintf('  TTG analysis missing from cache, computing now...\n');
            analyzer.compute_ttg_analysis();
            cache_needs_update = true;
        end
        
        % Check if distribution comparisons were loaded; if not, compute them
        if isempty(analyzer.distribution_comparisons)
            fprintf('  Distribution comparisons missing from cache, computing now...\n');
            analyzer.compute_distribution_comparisons();
            cache_needs_update = true;
        end
        
        % If we computed missing data, update the cache
        if cache_needs_update
            fprintf('  Updating cache with newly computed data...\n');
            analyzer.save_cache(cache_filepath);
            fprintf('  Cache updated successfully\n');
        end
    else
        if ~re_run_analysis
            fprintf('  No cache found, computing analysis...\n');
        else
            fprintf('  Re-running analysis (re_run_analysis=true)...\n');
        end
        
        % Run full analysis
        analyzer.analyze_all_clusters();
        
        % Compute TTG analysis
        analyzer.compute_ttg_analysis();
        
        % Compute distribution comparisons
        analyzer.compute_distribution_comparisons();
        
        % Print summary
        analyzer.print_summary();
        
        % Save cache
        [~, cache_filename] = fileparts(cache_filepath);
        fprintf('  Saving analysis data to cache: %s.mat\n', cache_filename);
        analyzer.save_cache(cache_filepath);
        fprintf('  Cache saved successfully\n');
    end
end
