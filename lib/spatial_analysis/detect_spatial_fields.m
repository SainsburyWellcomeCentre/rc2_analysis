function [longestField, numFields, fieldStart, fieldEnd] = detect_spatial_fields(rate_smooth, bin_centers, bin_size_cm, fieldFrac, minFieldBins)
% DETECT_SPATIAL_FIELDS Identify firing fields in spatial rate map
%
%   [longestField, numFields, fieldStart, fieldEnd] = detect_spatial_fields(...
%       rate_smooth, bin_centers, bin_size_cm, fieldFrac, minFieldBins)
%
%   Detects contiguous spatial regions where firing rate exceeds a threshold.
%   The threshold is defined as a fraction of the rate range (min to max).
%
%   FIELD DETECTION ALGORITHM:
%   1. Set threshold = min_rate + fieldFrac * (max_rate - min_rate)
%      This adapts to the modulation depth rather than using absolute peak
%      Example: rates 10-20 Hz with fieldFrac=0.7 → threshold = 10 + 0.7*10 = 17 Hz
%               rates 0-20 Hz with fieldFrac=0.7 → threshold = 0 + 0.7*20 = 14 Hz
%   2. Find all bins where firing rate >= threshold
%   3. Group consecutive bins into "fields" (contiguous regions above threshold)
%   4. Count the number of bins in each field
%   5. Record the longest field and total number of separate fields
%
%   EXAMPLE: If bins [5,6,7,8] and [15,16] are above threshold:
%     - numFields = 2 (two separate regions)
%     - fieldBins for field 1 = 4 bins (covering 8 cm with 2cm bins)
%     - fieldBins for field 2 = 2 bins (covering 4 cm with 2cm bins)
%     - longestField = 4 bins (the larger of the two)
%
%   For spatial tuning classification, typical requirements:
%     - longestField >= minFieldBins (e.g., >= 5 bins = 10 cm)
%     - numFields <= maxNumFields (typically 1, for single-field criterion)
%
%   Inputs:
%       rate_smooth  - Smoothed firing rate per spatial bin (Hz), 1D array
%       bin_centers  - Center position of each bin (cm), 1D array
%       bin_size_cm  - Size of each spatial bin (cm)
%       fieldFrac    - Fraction of rate range for threshold (0-1), typically 0.7
%       minFieldBins - Minimum bins required for a valid field (optional, not used in detection)
%
%   Outputs:
%       longestField - Number of bins in the longest field (0 if none found)
%       numFields    - Total number of separate fields detected
%       fieldStart   - Start position (cm) of the longest field (NaN if none)
%       fieldEnd     - End position (cm) of the longest field (NaN if none)
%
%   Example:
%       [longest, nFields, start, stop] = detect_spatial_fields(...
%           rate, bin_centers, 2, 0.7, 5);
%       if longest >= 5 && nFields == 1
%           fprintf('Single field from %.1f to %.1f cm\n', start, stop);
%       end
%
%   See also: skaggs_info, compute_shuffled_spatial_info

    % Compute threshold as fraction of rate range
    minRate = min(rate_smooth);
    maxRate = max(rate_smooth);
    thresh = minRate + fieldFrac * (maxRate - minRate);
    
    % Find bins above threshold
    above = find(rate_smooth >= thresh);
    
    if isempty(above)
        longestField = 0;
        numFields = 0;
        fieldStart = NaN;
        fieldEnd = NaN;
        return;
    end
    
    % Find contiguous segments by looking for breaks in the sequence
    % diff(above) > 1 indicates a gap between consecutive bins
    d = diff(above);
    breaks = [0, find(d > 1), length(above)];
    numFields = length(breaks) - 1;
    
    % Measure each field
    segLens = zeros(numFields, 1);
    segStarts = zeros(numFields, 1);
    segEnds = zeros(numFields, 1);
    
    for k = 1:numFields
        idx = (breaks(k)+1) : breaks(k+1);
        seg = above(idx);
        segLens(k) = length(seg);
        segStarts(k) = seg(1);   % First bin index in this field
        segEnds(k) = seg(end);   % Last bin index in this field
    end
    
    % Find longest field
    [longestField, longestIdx] = max(segLens);
    
    % Store the boundaries of the longest field (in cm, bin edges not centers)
    fieldStartIdx = segStarts(longestIdx);
    fieldEndIdx = segEnds(longestIdx);
    fieldStart = bin_centers(fieldStartIdx) - bin_size_cm/2;
    fieldEnd = bin_centers(fieldEndIdx) + bin_size_cm/2;
end
