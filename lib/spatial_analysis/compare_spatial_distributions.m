function results = compare_spatial_distributions(rate_per_trial_long, rate_per_trial_short, bin_centers_long, bin_centers_short, ttg_rate_long, ttg_rate_short, rate_per_trial_long_downsampled, bin_centers_long_downsampled)
% COMPARE_SPATIAL_DISTRIBUTIONS Compare firing rate distributions between conditions
%
%   results = compare_spatial_distributions(rate_per_trial_long, rate_per_trial_short, ...
%                                           bin_centers_long, bin_centers_short, ttg_rate_long, ttg_rate_short, ...
%                                           rate_per_trial_long_downsampled, bin_centers_long_downsampled)
%
%   Performs four statistical tests using Kruskal-Wallis (non-parametric ANOVA):
%   1. Absolute position comparison: Second half of long trials vs short trials
%   2. Relative position comparison: Downsampled long trials (normalized) vs short trials
%   3. Absolute shifted comparison: First half of long trials (0-60 cm) vs short trials
%   4. TTG comparison: Time-to-goal firing rates for long vs short trials
%
%   Uses Kruskal-Wallis test to compare firing rate distributions across spatial bins.
%
%   Inputs:
%       rate_per_trial_long  - Matrix (n_trials_long × n_bins_long) of firing rates for long trials
%       rate_per_trial_short - Matrix (n_trials_short × n_bins_short) of firing rates for short trials
%       bin_centers_long     - Vector of bin centers for long trials (cm)
%       bin_centers_short    - Vector of bin centers for short trials (cm)
%       ttg_rate_long        - (Optional) Matrix (n_trials_long × n_ttg_bins) of TTG firing rates for long trials
%       ttg_rate_short       - (Optional) Matrix (n_trials_short × n_ttg_bins) of TTG firing rates for short trials
%       rate_per_trial_long_downsampled - (Optional) Matrix (n_trials_long × n_bins_downsampled) with half spatial resolution
%       bin_centers_long_downsampled    - (Optional) Vector of bin centers for downsampled long trials (cm)
%
%   Outputs:
%       results - Struct with fields:
%           .absolute - Struct with:
%               .same_distribution - Logical, true if distributions are not significantly different
%               .p_value          - P-value from Kruskal-Wallis test
%               .chi2_stat        - Chi-square test statistic
%           .relative - Struct with same fields for relative position comparison
%           .absolute_shifted - Struct with same fields for first half comparison
%           .ttg - Struct with same fields for TTG comparison
%
%   Method:
%       For each comparison, we flatten the trial×bin matrices and perform
%       Kruskal-Wallis test to determine if the firing rate distributions
%       differ between conditions. This is a non-parametric test that doesn't
%       assume normality of the data.
    
    results = struct();
    
    %% Test 1: Absolute position comparison
    % Compare second half of long trials (60-120 cm) with short trials (60-120 cm plotted)
    
    % Find bins corresponding to 60-120 cm in long trials
    second_half_bins = bin_centers_long >= 60;
    
    if sum(second_half_bins) == 0
        % No data in second half
        results.absolute.same_distribution = NaN;
        results.absolute.p_value = NaN;
        results.absolute.chi2_stat = NaN;
        results.absolute.error = 'No bins in second half of long trials';
    else
        % Extract second half of long trials
        rate_long_second_half = rate_per_trial_long(:, second_half_bins);
        n_bins_long_second = size(rate_long_second_half, 2);
        n_bins_short = size(rate_per_trial_short, 2);
        
        % Check if we have matching number of bins
        if n_bins_long_second ~= n_bins_short
            % Interpolate to match bin counts (this should generally match for 2cm bins)
            warning('Bin count mismatch: long second half has %d bins, short has %d bins. Interpolating.', ...
                    n_bins_long_second, n_bins_short);
            
            % Interpolate each trial's profile to match short trial bin count
            rate_long_second_half_interp = zeros(size(rate_long_second_half, 1), n_bins_short);
            x_long = linspace(0, 1, n_bins_long_second);
            x_short = linspace(0, 1, n_bins_short);
            
            for t = 1:size(rate_long_second_half, 1)
                rate_long_second_half_interp(t, :) = interp1(x_long, rate_long_second_half(t, :), x_short, 'linear', 'extrap');
            end
            
            rate_long_second_half = rate_long_second_half_interp;
        end
        
        % Perform Kruskal-Wallis test
        results.absolute = kruskal_wallis_test(rate_long_second_half, rate_per_trial_short);
    end
    
    %% Test 2: Relative position comparison
    % Use pre-downsampled long trials (with half spatial resolution, downsampled BEFORE smoothing)
    % Then normalize position to 0-100% for both conditions
    
    if nargin >= 8 && ~isempty(rate_per_trial_long_downsampled) && ~isempty(bin_centers_long_downsampled)
        % Use the downsampled data that was computed with half spatial resolution
        n_bins_long_ds = size(rate_per_trial_long_downsampled, 2);
        n_bins_short = size(rate_per_trial_short, 2);
        
        % Interpolate both to a common normalized grid (e.g., 50 bins spanning 0-100%)
        n_norm_bins = 50;
        x_norm = linspace(0, 1, n_norm_bins);
        
        % Interpolate downsampled long trials
        rate_long_norm = zeros(size(rate_per_trial_long_downsampled, 1), n_norm_bins);
        x_long_ds = linspace(0, 1, n_bins_long_ds);
        
        for t = 1:size(rate_per_trial_long_downsampled, 1)
            rate_long_norm(t, :) = interp1(x_long_ds, rate_per_trial_long_downsampled(t, :), x_norm, 'linear', 'extrap');
        end
        
        % Interpolate short trials
        rate_short_norm = zeros(size(rate_per_trial_short, 1), n_norm_bins);
        x_short = linspace(0, 1, n_bins_short);
        
        for t = 1:size(rate_per_trial_short, 1)
            rate_short_norm(t, :) = interp1(x_short, rate_per_trial_short(t, :), x_norm, 'linear', 'extrap');
        end
        
        % Perform Kruskal-Wallis test on normalized data
        results.relative = kruskal_wallis_test(rate_long_norm, rate_short_norm);
    else
        % Fallback to old method if downsampled data not provided
        % Extract even bins from long trials (every other bin)
        even_bins = 2:2:length(bin_centers_long);
        
        if length(even_bins) == 0
            results.relative.same_distribution = NaN;
            results.relative.p_value = NaN;
            results.relative.chi2_stat = NaN;
            results.relative.error = 'No even bins in long trials';
        else
            rate_long_even = rate_per_trial_long(:, even_bins);
            n_bins_long_even = size(rate_long_even, 2);
            n_bins_short = size(rate_per_trial_short, 2);
            
            % Interpolate both to a common normalized grid (e.g., 50 bins spanning 0-100%)
            n_norm_bins = 50;
            x_norm = linspace(0, 1, n_norm_bins);
            
            % Interpolate long trials (even bins)
            rate_long_norm = zeros(size(rate_long_even, 1), n_norm_bins);
            x_long_even = linspace(0, 1, n_bins_long_even);
            
            for t = 1:size(rate_long_even, 1)
                rate_long_norm(t, :) = interp1(x_long_even, rate_long_even(t, :), x_norm, 'linear', 'extrap');
            end
            
            % Interpolate short trials
            rate_short_norm = zeros(size(rate_per_trial_short, 1), n_norm_bins);
            x_short = linspace(0, 1, n_bins_short);
            
            for t = 1:size(rate_per_trial_short, 1)
                rate_short_norm(t, :) = interp1(x_short, rate_per_trial_short(t, :), x_norm, 'linear', 'extrap');
            end
            
            % Perform Kruskal-Wallis test on normalized data
            results.relative = kruskal_wallis_test(rate_long_norm, rate_short_norm);
        end
    end
    
    %% Test 3: Absolute shifted comparison
    % Compare first half of long trials (0-60 cm) with short trials (60-120 cm plotted)
    
    % Find bins corresponding to 0-60 cm in long trials
    first_half_bins = bin_centers_long < 60;
    
    if sum(first_half_bins) == 0
        % No data in first half
        results.absolute_shifted.same_distribution = NaN;
        results.absolute_shifted.p_value = NaN;
        results.absolute_shifted.chi2_stat = NaN;
        results.absolute_shifted.error = 'No bins in first half of long trials';
    else
        % Extract first half of long trials
        rate_long_first_half = rate_per_trial_long(:, first_half_bins);
        n_bins_long_first = size(rate_long_first_half, 2);
        n_bins_short = size(rate_per_trial_short, 2);
        
        % Check if we have matching number of bins
        if n_bins_long_first ~= n_bins_short
            % Interpolate to match bin counts
            rate_long_first_half_interp = zeros(size(rate_long_first_half, 1), n_bins_short);
            x_long = linspace(0, 1, n_bins_long_first);
            x_short = linspace(0, 1, n_bins_short);
            
            for t = 1:size(rate_long_first_half, 1)
                rate_long_first_half_interp(t, :) = interp1(x_long, rate_long_first_half(t, :), x_short, 'linear', 'extrap');
            end
            
            rate_long_first_half = rate_long_first_half_interp;
        end
        
        % Perform Kruskal-Wallis test
        results.absolute_shifted = kruskal_wallis_test(rate_long_first_half, rate_per_trial_short);
    end
    
    %% Test 4: TTG comparison
    % Compare time-to-goal firing rates for long vs short trials
    
    if nargin >= 6 && ~isempty(ttg_rate_long) && ~isempty(ttg_rate_short)
        % Both TTG rate matrices should already be on the same normalized time-to-goal axis
        % Perform Kruskal-Wallis test directly
        results.ttg = kruskal_wallis_test(ttg_rate_long, ttg_rate_short);
    else
        % No TTG data provided
        results.ttg.same_distribution = NaN;
        results.ttg.p_value = NaN;
        results.ttg.chi2_stat = NaN;
        results.ttg.error = 'No TTG data provided';
    end
end


function result = kruskal_wallis_test(data1, data2)
    % Perform Kruskal-Wallis test (non-parametric ANOVA)
    %
    % data1: n_trials1 × n_bins
    % data2: n_trials2 × n_bins
    %
    % Tests if the firing rate distributions differ between conditions
    
    result = struct();
    
    % Flatten the matrices to create vectors of all firing rate values
    data1_flat = data1(:);
    data2_flat = data2(:);
    
    % Remove any NaN or Inf values
    valid1 = isfinite(data1_flat);
    valid2 = isfinite(data2_flat);
    data1_flat = data1_flat(valid1);
    data2_flat = data2_flat(valid2);
    
    % Create group labels
    group1 = ones(length(data1_flat), 1);
    group2 = 2 * ones(length(data2_flat), 1);
    
    % Combine data and groups
    combined_data = [data1_flat; data2_flat];
    combined_groups = [group1; group2];
    
    % Perform Kruskal-Wallis test
    try
        [p_value, tbl, stats] = kruskalwallis(combined_data, combined_groups, 'off');
        % The chi-square statistic is in the table output, not the stats struct
        if istable(tbl)
            chi2_stat = tbl{1, 5};  % Chi-square value is in column 5, row 1 (Groups row)
        else
            chi2_stat = tbl{2, 5};  % For cell array, Groups row is typically row 2
        end
    catch ME
        warning('Kruskal-Wallis test failed: %s', ME.message);
        p_value = NaN;
        chi2_stat = NaN;
    end
    
    % Determine if distributions are the same (fail to reject null)
    alpha = 0.05;
    same_distribution = p_value > alpha;
    
    % Store results
    result.same_distribution = same_distribution;
    result.p_value = p_value;
    result.chi2_stat = chi2_stat;
end
