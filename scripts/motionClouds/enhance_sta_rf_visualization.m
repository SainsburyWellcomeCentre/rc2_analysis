% Enhance STA/RF images for better visualization (Niell & Stryker-inspired)
% - Robust z-scoring (median/MAD)
% - Difference-of-Gaussians (center-surround enhancement)
% - Mild smoothing
% - Contrast clipping and scaling to [0,1] or diverging BWR colormap
% - Optional contour overlays at z-levels
%
% This script scans the STA PNGs produced by sta_rf_analysis.m and writes
% enhanced versions into a sibling folder.

close all;

% Parameters
experiment_groups    = {'passive_same_luminance_mc'}; % kept for default path layout
output_subdirs_img   = {'motionClouds','passive_same_luminance_mc','sta_rf','images'};

% If you want to override the input directory, set input_dir_override to a path
input_dir_override = '';

% Visualization parameters
gauss_sigma_small = 1.0;   % pixels (center)
gauss_sigma_large = 5.0;   % pixels (surround)
gauss_sigma_post  = 0.8;   % additional smoothing
z_clip            = 3.0;   % clip range [-z_clip, z_clip]
color_mode        = 'bwr'; % 'gray' or 'bwr' (blue-negative, white-zero, red-positive)
contour_levels    = [2, 3, 4]; % z-levels for contour overlays
write_contours    = true;

% Resolve paths
script_dir = fileparts(mfilename('fullpath'));
if isempty(input_dir_override)
	input_dir = fullfile(script_dir, output_subdirs_img{:});
else
	input_dir = input_dir_override;
end
if ~isfolder(input_dir)
	error('Input directory not found: %s', input_dir);
end

out_dir_enhanced = fullfile(input_dir, '..', 'images_enhanced');
out_dir_enhanced = char(java.io.File(out_dir_enhanced).getCanonicalPath()); %#ok<J2F>
if ~isfolder(out_dir_enhanced), mkdir(out_dir_enhanced); end

if write_contours
	out_dir_contours = fullfile(input_dir, '..', 'images_contours');
	out_dir_contours = char(java.io.File(out_dir_contours).getCanonicalPath()); %#ok<J2F>
	if ~isfolder(out_dir_contours), mkdir(out_dir_contours); end
end

% Gather PNGs
pngs = dir(fullfile(input_dir, '*.png'));
fprintf('Found %d PNG STA images in %s\n', numel(pngs), input_dir);
if isempty(pngs)
	return;
end

% Processing helpers
use_imgauss = exist('imgaussfilt','file') == 2;

for i = 1:numel(pngs)
	try
		in_path = fullfile(input_dir, pngs(i).name);
		I = imread(in_path);
		if ndims(I) == 3
			I = rgb2gray(I);
		end
		I = im2single(I);

		% Robust center and scale (z-score)
		v = I(:);
		med = median(v);
		mad = median(abs(v - med));
		scale = 1.4826 * mad; % ~= std for Gaussian
		if scale <= eps
			scale = std(v);
			if scale <= eps, scale = 1; end
		end
		Z = (I - med) / scale;

		% Difference-of-Gaussians to enhance center-surround structure
		if use_imgauss
			Zc = imgaussfilt(Z, gauss_sigma_small);
			Zs = imgaussfilt(Z, gauss_sigma_large);
		else
			Zc = local_gaussfilt(Z, gauss_sigma_small);
			Zs = local_gaussfilt(Z, gauss_sigma_large);
		end
		Zdog = Zc - Zs;

		% Optional light extra smoothing
		if gauss_sigma_post > 0
			if use_imgauss
				Zdog = imgaussfilt(Zdog, gauss_sigma_post);
			else
				Zdog = local_gaussfilt(Zdog, gauss_sigma_post);
			end
		end

		% Clip
		Zclip = max(min(Zdog, z_clip), -z_clip);

		% Write enhanced image (grayscale or BWR)
		[~, base, ~] = fileparts(pngs(i).name);
		out_path = fullfile(out_dir_enhanced, sprintf('%s_enhanced.png', base));
		if strcmpi(color_mode,'bwr')
			M = local_bwr_colormap(256);
			A = (Zclip + z_clip) / (2 * z_clip); % map [-z_clip,z_clip] -> [0,1]
			A = min(max(A,0),1);
			RGB = ind2rgb(1 + floor(A*255), M);
			imwrite(RGB, out_path);
		else
			E = (Zclip + z_clip) / (2 * z_clip);
			imwrite(E, out_path);
		end

		% Optional contour overlays at specified z-levels
		if write_contours
			f = figure('Visible','off','Color','w');
			subplot(1,1,1);
			imagesc(Zdog); axis image off; 
			if strcmpi(color_mode,'bwr')
				colormap(gca, local_bwr_colormap(256));
				caxis([-z_clip z_clip]);
			else
				colormap(gca, 'gray');
			end
			hold on;
			for lv = contour_levels
				contour(Zdog, [lv lv], 'r', 'LineWidth', 1.0);
				contour(Zdog, [-lv -lv], 'b', 'LineWidth', 1.0);
			end
			title(sprintf('%s  (z-DoG, contours: %s)', base, mat2str(contour_levels)),'Interpreter','none');
			outc = fullfile(out_dir_contours, sprintf('%s_contours.png', base));
			exportgraphics(gca, outc, 'Resolution', 200);
			close(f);
		end

		fprintf('Processed %d/%d: %s\n', i, numel(pngs), pngs(i).name);
	catch ME
		warning('Failed on %s: %s', pngs(i).name, ME.message);
	end
end

fprintf('Done. Enhanced images in: %s\n', out_dir_enhanced);
if write_contours
	fprintf('Contour overlays in: %s\n', out_dir_contours);
end


function O = local_gaussfilt(I, sigma)
% Fallback Gaussian filter using fspecial if imgaussfilt is absent
	if sigma <= 0
		O = I; return;
	end
	sz = max(3, 2*ceil(3*sigma) + 1);
	h = fspecial('gaussian', [sz sz], sigma);
	O = conv2(I, h, 'same');
end

function M = local_bwr_colormap(n)
% Blue-White-Red diverging colormap with white at zero
	if nargin < 1 || isempty(n), n = 256; end
	n = max(3, round(n));
	half = floor(n/2);
	% Blue to white
	blue = [0 0 1]; white = [1 1 1]; red = [1 0 0];
	B = [linspace(blue(1), white(1), half)', linspace(blue(2), white(2), half)', linspace(blue(3), white(3), half)'];
	% White to red (ensure total n)
	rem = n - half;
	R = [linspace(white(1), red(1), rem)', linspace(white(2), red(2), rem)', linspace(white(3), red(3), rem)'];
	M = [B; R];
end



