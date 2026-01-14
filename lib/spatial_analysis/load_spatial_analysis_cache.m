function analyzer = load_spatial_analysis_cache(clusters, sessions, cache_filepath)
% LOAD_SPATIAL_ANALYSIS_CACHE Load spatial analysis from cache file
%
% SYNTAX:
%   analyzer = load_spatial_analysis_cache(clusters, sessions, cache_filepath)
%
% INPUTS:
%   clusters         - Array of cluster objects
%   sessions         - Cell array of session objects
%   cache_filepath   - Full path to cache file
%
% OUTPUTS:
%   analyzer         - SpatialTuningAnalyzer instance with loaded data
%
% DESCRIPTION:
%   Loads spatial analysis data from cache file and reconstructs the analyzer.
%   All parameters and results are restored from the cached data.

    [~, cache_filename] = fileparts(cache_filepath);
    fprintf('  Loading cached analysis data from: %s.mat\n', cache_filename);
    
    % Create analyzer without parameters (will be loaded from cache)
    analyzer = SpatialTuningAnalyzer(clusters, sessions);
    
    % Load all data from cache (including parameters)
    analyzer.load_cache(cache_filepath);
    fprintf('  Loaded cached data for %d clusters\n', length(analyzer.cluster_ids));
    
    % Reconstruct trial groups for plotting (not stored in cache)
    analyzer.collect_trials();
    
    % Verify bin config exists (should be reconstructed in load_cache)
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
end
