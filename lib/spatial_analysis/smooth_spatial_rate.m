function rate_smooth = smooth_spatial_rate(rate, bin_size_cm, gauss_sigma_cm)
% SMOOTH_SPATIAL_RATE Apply Gaussian smoothing with reflection padding
%
%   rate_smooth = smooth_spatial_rate(rate, bin_size_cm, gauss_sigma_cm)
%
%   Applies Gaussian smoothing to spatial firing rate data using reflection
%   padding at the edges to minimize boundary artifacts.
%
%   Inputs:
%       rate           - 1D array of firing rates per spatial bin (Hz)
%       bin_size_cm    - Size of each spatial bin (cm)
%       gauss_sigma_cm - Standard deviation of Gaussian kernel (cm)
%
%   Outputs:
%       rate_smooth    - Smoothed firing rate array (same size as input)
%
%   Example:
%       rate_smooth = smooth_spatial_rate(rate, 2, 8);  % 2-cm bins, 8-cm sigma

    % Handle NaN values by replacing with 0 for smoothing
    rate_for_smooth = rate;
    rate_for_smooth(isnan(rate_for_smooth)) = 0;
    
    % Create Gaussian kernel
    kernel_size_cm = 8 * gauss_sigma_cm;  % Total kernel size (8 sigma covers >99.99%)
    kernel_samples = round(kernel_size_cm / bin_size_cm);
    
    % Ensure kernel has odd number of samples for symmetry
    if mod(kernel_samples, 2) == 0
        kernel_samples = kernel_samples + 1;
    end
    
    x = linspace(-kernel_samples/2, kernel_samples/2, kernel_samples);
    sigma_in_bins = gauss_sigma_cm / bin_size_cm;
    kernel = exp(-x.^2 / (2 * sigma_in_bins^2));
    kernel = kernel / sum(kernel);  % Normalize to sum to 1
    
    % Pad the rate data using reflection to avoid edge effects
    pad_size = floor(length(kernel) / 2);
    
    % Handle case where rate array is shorter than pad size
    if length(rate_for_smooth) <= pad_size
        % For very short arrays, just use available data for padding
        pad_left = fliplr(rate_for_smooth);
        pad_right = fliplr(rate_for_smooth);
        rate_padded = [pad_left, rate_for_smooth, pad_right];
    else
        rate_padded = [fliplr(rate_for_smooth(1:pad_size)), ...
                       rate_for_smooth, ...
                       fliplr(rate_for_smooth(end-pad_size+1:end))];
    end
    
    % Apply convolution
    rate_smooth_padded = conv(rate_padded, kernel, 'same');
    
    % Extract the central portion (remove padding)
    if length(rate_for_smooth) <= pad_size
        % Handle short array case
        rate_smooth = rate_smooth_padded(length(pad_left)+1 : length(pad_left)+length(rate_for_smooth));
    else
        rate_smooth = rate_smooth_padded(pad_size + 1 : end - pad_size);
    end
end
