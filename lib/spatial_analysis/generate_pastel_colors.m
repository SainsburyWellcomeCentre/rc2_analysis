function colors = generate_pastel_colors(n)
% GENERATE_PASTEL_COLORS Generate n distinct pastel colors
%
%   colors = generate_pastel_colors(n)
%
%   Generates n visually distinct pastel colors using HSV color space.
%   Pastel colors have high value (brightness) and low-to-medium saturation.
%
%   Inputs:
%       n - Number of colors to generate
%
%   Outputs:
%       colors - n x 3 matrix of RGB colors (values in [0, 1])
%
%   Example:
%       colors = generate_pastel_colors(5);
%       for i = 1:5
%           plot(1:10, i*ones(1,10), 'Color', colors(i,:), 'LineWidth', 3);
%       end

    % Generate evenly spaced hues
    hues = linspace(0, 1, n+1);
    hues = hues(1:end-1); % Remove duplicate at end
    
    colors = zeros(n, 3);
    for i = 1:n
        % Use moderate saturation (0.4-0.6) and high value (0.85-0.95) for pastel effect
        hsv_color = [hues(i), 0.5, 0.9];
        colors(i, :) = hsv2rgb(hsv_color);
    end
end
