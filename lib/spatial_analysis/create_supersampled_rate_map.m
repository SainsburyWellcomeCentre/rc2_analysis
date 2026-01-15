function [supersampled_map, supersampled_alpha, x_coords, y_coords] = create_supersampled_rate_map(rate_data, alpha_matrix, x_edges, y_edges, resolution_factor)
% CREATE_SUPERSAMPLED_RATE_MAP Create high-resolution map for variable-width bins
%
% SYNTAX:
%   [supersampled_map, supersampled_alpha, x_coords, y_coords] = ...
%       create_supersampled_rate_map(rate_data, alpha_matrix, x_edges, y_edges, resolution_factor)
%
% DESCRIPTION:
%   Creates a super-sampled matrix where each bin is expanded to occupy
%   a number of pixels proportional to its actual size. This allows
%   imagesc() to display variable-width bins correctly (e.g., for
%   equal-occupancy velocity/acceleration bins).
%
% INPUTS:
%   rate_data         - n_x_bins × n_y_bins matrix of firing rates
%   alpha_matrix      - n_x_bins × n_y_bins matrix of alpha values (0-1)
%   x_edges           - (n_x_bins+1) × 1 vector of x-axis bin edges
%   y_edges           - (n_y_bins+1) × 1 vector of y-axis bin edges  
%   resolution_factor - (Optional) Multiplier for pixel resolution (default: 10)
%
% OUTPUTS:
%   supersampled_map   - High-res matrix with replicated rate values
%   supersampled_alpha - High-res matrix with replicated alpha values
%   x_coords           - X-coordinates for imagesc (bin centers)
%   y_coords           - Y-coordinates for imagesc (bin centers)
%
% EXAMPLE:
%   % Create supersampled map
%   [map, alpha, x, y] = create_supersampled_rate_map(rates, alphas, pos_edges, vel_edges);
%   % Display with imagesc
%   h = imagesc(x, y, map');
%   set(h, 'AlphaData', alpha');
%   set(gca, 'YDir', 'normal');

    if nargin < 5
        resolution_factor = 10;  % Default: 10x oversampling
    end
    
    % Calculate bin widths
    x_widths = diff(x_edges);
    y_widths = diff(y_edges);
    
    % Find minimum bin width as base resolution
    min_x_width = min(x_widths);
    min_y_width = min(y_widths);
    
    % Calculate pixel resolution (finer than smallest bin)
    x_pixel_size = min_x_width / resolution_factor;
    y_pixel_size = min_y_width / resolution_factor;
    
    % Calculate how many pixels each bin should occupy
    x_pixels_per_bin = round(x_widths / x_pixel_size);
    y_pixels_per_bin = round(y_widths / y_pixel_size);
    
    % Total size of supersampled matrix
    total_x_pixels = sum(x_pixels_per_bin);
    total_y_pixels = sum(y_pixels_per_bin);
    
    % Pre-allocate supersampled matrices
    supersampled_map = nan(total_x_pixels, total_y_pixels);
    supersampled_alpha = zeros(total_x_pixels, total_y_pixels);
    
    % Fill supersampled matrix
    x_pixel_start = 1;
    for i = 1:length(x_pixels_per_bin)
        x_pixel_end = x_pixel_start + x_pixels_per_bin(i) - 1;
        
        y_pixel_start = 1;
        for j = 1:length(y_pixels_per_bin)
            y_pixel_end = y_pixel_start + y_pixels_per_bin(j) - 1;
            
            % Replicate the rate value and alpha across this bin's pixels
            supersampled_map(x_pixel_start:x_pixel_end, y_pixel_start:y_pixel_end) = rate_data(i, j);
            supersampled_alpha(x_pixel_start:x_pixel_end, y_pixel_start:y_pixel_end) = alpha_matrix(i, j);
            
            y_pixel_start = y_pixel_end + 1;
        end
        
        x_pixel_start = x_pixel_end + 1;
    end
    
    % Create coordinate vectors for imagesc (pixel centers)
    % These represent the actual physical coordinates
    x_coords = zeros(1, total_x_pixels);
    x_pixel_idx = 1;
    for i = 1:length(x_pixels_per_bin)
        bin_start = x_edges(i);
        bin_width = x_widths(i);
        n_pixels = x_pixels_per_bin(i);
        
        % Calculate center position for each pixel within this bin
        pixel_centers = linspace(bin_start + bin_width/(2*n_pixels), ...
                                 bin_start + bin_width - bin_width/(2*n_pixels), ...
                                 n_pixels);
        x_coords(x_pixel_idx:x_pixel_idx+n_pixels-1) = pixel_centers;
        x_pixel_idx = x_pixel_idx + n_pixels;
    end
    
    y_coords = zeros(1, total_y_pixels);
    y_pixel_idx = 1;
    for j = 1:length(y_pixels_per_bin)
        bin_start = y_edges(j);
        bin_width = y_widths(j);
        n_pixels = y_pixels_per_bin(j);
        
        % Calculate center position for each pixel within this bin
        pixel_centers = linspace(bin_start + bin_width/(2*n_pixels), ...
                                 bin_start + bin_width - bin_width/(2*n_pixels), ...
                                 n_pixels);
        y_coords(y_pixel_idx:y_pixel_idx+n_pixels-1) = pixel_centers;
        y_pixel_idx = y_pixel_idx + n_pixels;
    end
end
