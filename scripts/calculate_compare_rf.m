% Calculate and compare receptive field (RF) sizes under two conditions:
%   1) Goggles (session rec1, no warp applied)
%   2) Screens (session rec4, warped to mp_300)
%
% Both sessions belong to the same probe recording (experiment group:
% sparse_noise_goggles). RFs are computed as spike-triggered averages (STA)
% for white (on) and black (off) sparse noise rectangles.
%
% How are RFs calculated?
% 1. The starting time of every new matrix presentation is found from the
%    photodiode signal.
% 2. We calculate the total spike count R per cluster at every stimulus
%    presentation.
% 3. The STA is calculated for white and black rectangles.
% 4. A mean and std shuffled STA is generated (100 shuffles).
% 5. RF_raw = STA > (STA_shuffled_mean + threshold * STA_shuffled_std)
%    Connected squares are retained with bwareaopen.
% 6. Width and height of the RF are measured and converted to degrees.
%
% Output PDFs:
%   Goggles (2x2 layout, no warped panel):
%       C:\Users\lee\Documents\mvelez\figures\new_RFs\goggles\
%   Screens (2x3 layout, includes warped-to-screen panel):
%       C:\Users\lee\Documents\mvelez\figures\new_RFs\screens\

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
degs_per_square = 5;   % degrees per square in the stimulus (constant)
threshold       = 2.5;   % multiplicative factor applied to shuffled STA std

% ON/OFF colormaps: white→red (responses to white/ON squares)
%                  white→blue (responses to black/OFF squares)
cmap_on  = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)'];  % white→red
cmap_off = [linspace(1,0,256)', linspace(1,0,256)', linspace(1,1,256)'];  % white→blue

% -------------------------------------------------------------------------
% Stimulus and warp file paths
% -------------------------------------------------------------------------
goggles_stim_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_wisecoco.mat';
screens_stim_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_mp_300_20260323.mat';
screens_warp_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\warp\warp_mp_300_eye_to_screen_dist_15cm.mat';

% -------------------------------------------------------------------------
% Output directories
% -------------------------------------------------------------------------
goggles_save_dir = 'C:\Users\lee\Documents\mvelez\figures\new_RFs\goggles';
screens_save_dir = 'C:\Users\lee\Documents\mvelez\figures\new_RFs\screens';
if ~exist(goggles_save_dir, 'dir'), mkdir(goggles_save_dir); end
if ~exist(screens_save_dir, 'dir'), mkdir(screens_save_dir); end

% -------------------------------------------------------------------------
% Select recordings
% -------------------------------------------------------------------------
ctl       = RC2Analysis();
probe_ids = ctl.get_probe_ids('sparse_noise_goggles');

% -------------------------------------------------------------------------
% Storage for RF dimensions across all probes
% -------------------------------------------------------------------------
RF_w_goggles = [];
RF_h_goggles = [];
RF_w_screens = [];
RF_h_screens = [];

% -------------------------------------------------------------------------
% Pre-compute stimulus matrices (same for all probes)
% -------------------------------------------------------------------------

% --- Goggles stimulus ---
% The wisecoco format does not store h_pix/w_pix explicitly; derive from
% n_y/n_x (number of squares) * pix_square (pixels per square).
stim_g      = load(goggles_stim_fname);
n_stim_g    = length(stim_g.cols);
stim_mtx_g  = 0.5 * ones([stim_g.grid_size, n_stim_g]);
for stim_i_g = 1 : n_stim_g
    t = 0.5 * ones(stim_g.grid_size);
    t(stim_g.locations{stim_i_g}) = stim_g.cols{stim_i_g};
    stim_mtx_g(:, :, stim_i_g) = t;
end
% x_border/y_border hold pixel positions of grid lines; the last value
% is the far edge of the final square (~399.5 for a 400-pixel display).
h_pix_g      = round(stim_g.y_border(end));
w_pix_g      = round(stim_g.x_border(end));
n_y_g        = size(stim_mtx_g, 1);
n_x_g        = size(stim_mtx_g, 2);
stim_mtx_g   = reshape(stim_mtx_g, [], size(stim_mtx_g, 3));  % (n_y*n_x) x t
stim_white_g = stim_mtx_g == 1;
stim_black_g = stim_mtx_g == 0;

% --- Screens stimulus ---
% Same format as wisecoco: no h_pix/w_pix field, derive from x/y_border.
stim_s      = load(screens_stim_fname);
n_stim_s    = length(stim_s.cols);
stim_mtx_s  = 0.5 * ones([stim_s.grid_size, n_stim_s]);
for stim_i_s = 1 : n_stim_s
    t = 0.5 * ones(stim_s.grid_size);
    t(stim_s.locations{stim_i_s}) = stim_s.cols{stim_i_s};
    stim_mtx_s(:, :, stim_i_s) = t;
end
h_pix_s      = round(stim_s.y_border(end));
w_pix_s      = round(stim_s.x_border(end));
n_y_s        = size(stim_mtx_s, 1);
n_x_s        = size(stim_mtx_s, 2);
stim_mtx_s   = reshape(stim_mtx_s, [], size(stim_mtx_s, 3));  % (n_y*n_x) x t
stim_white_s = stim_mtx_s == 1;
stim_black_s = stim_mtx_s == 0;

% -------------------------------------------------------------------------
% Pre-compute screens warping indices (same geometry for all probes)
% -------------------------------------------------------------------------
load(screens_warp_fname, 'scal');
[Xq_s, Yq_s] = meshgrid(1 : w_pix_s + 1, 1 : h_pix_s + 1);
orixi_s = round(interp2(scal.vcoords(:, :, 1), scal.vcoords(:, :, 2), ...
                        scal.tcoords(:, :, 1), Xq_s, Yq_s));
oriyi_s = round(interp2(scal.vcoords(:, :, 1), scal.vcoords(:, :, 2), ...
                        scal.tcoords(:, :, 2), Xq_s, Yq_s));
I_s = sub2ind(size(orixi_s), oriyi_s, orixi_s);

% =========================================================================
% Main loop over probe recordings
% =========================================================================
for probe_i = 1 : length(probe_ids)

    this_probe = probe_ids(probe_i)

    % ---------------------------------------------------------------------
    % Load data and select sessions
    % ---------------------------------------------------------------------
    data            = ctl.load_formatted_data(this_probe{1});
    goggles_session = data.sessions{1, 1};   % rec1 - goggles
    screens_session = data.sessions{1, 4};   % rec4 - screens

    % Get animal ID for photodiode parameter lookup
    animal_id = get_animal_id(this_probe);

    % Find stimulus onset times from photodiode signal
    starts_g = photodiode_times(animal_id, goggles_session, "sparse_noise_goggles");
    starts_s = photodiode_times(animal_id, screens_session, "sparse_noise_screens");

    % ---------------------------------------------------------------------
    % Initialise STA storage for this probe
    % ---------------------------------------------------------------------
    sta_white_g_store     = [];   sta_black_g_store     = [];
    sta_white_g_shuf_mean = [];   sta_black_g_shuf_mean = [];
    sta_white_g_shuf_std  = [];   sta_black_g_shuf_std  = [];

    sta_white_s_store     = [];   sta_black_s_store     = [];
    sta_white_s_shuf_mean = [];   sta_black_s_shuf_mean = [];
    sta_white_s_shuf_std  = [];   sta_black_s_shuf_std  = [];

    % ---------------------------------------------------------------------
    % Compute STAs for every VISp cluster
    % NOTE: restrict_to_visp=false includes clusters without anatomy data
    % (labelled 'unknownLocation'). Update to true once anatomy is available.
    % ---------------------------------------------------------------------
    clusters = data.VISp_clusters(true, 'any', false);

    for clust_i = 1 : length(clusters)

        all_spike_times = sort(clusters(clust_i).spike_times);

        % ----- Goggles: spike count per stimulus frame -----
        R_g = zeros(1, length(starts_g));
        for i = 1 : length(starts_g) - 1
            R_g(i) = sum(all_spike_times > goggles_session.probe_t(starts_g(i)) & ...
                         all_spike_times < goggles_session.probe_t(starts_g(i + 1)));
        end
        R_g(end) = sum(all_spike_times > goggles_session.probe_t(starts_g(end)) & ...
                       all_spike_times < goggles_session.probe_t(starts_g(end)) + 0.25);

        % Goggles STA - white
        sta_white_g = bsxfun(@times, stim_white_g, R_g(:)');
        sta_white_g = reshape(mean(sta_white_g, 2), n_y_g, n_x_g);

        sta_white_g_shuffled = nan(n_y_g, n_x_g, 100);
        for sh_i = 1 : 100
            I = randperm(length(R_g));
            st = mean(bsxfun(@times, stim_white_g, R_g(I)), 2);
            sta_white_g_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
        end

        % Goggles STA - black
        sta_black_g = bsxfun(@times, stim_black_g, R_g(:)');
        sta_black_g = reshape(mean(sta_black_g, 2), n_y_g, n_x_g);

        sta_black_g_shuffled = nan(n_y_g, n_x_g, 100);
        for sh_i = 1 : 100
            I = randperm(length(R_g));
            st = mean(bsxfun(@times, stim_black_g, R_g(I)), 2);
            sta_black_g_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
        end

        sta_white_g_store     = cat(3, sta_white_g_store,     sta_white_g);
        sta_black_g_store     = cat(3, sta_black_g_store,     sta_black_g);
        sta_white_g_shuf_mean = cat(3, sta_white_g_shuf_mean, mean(sta_white_g_shuffled, 3));
        sta_black_g_shuf_mean = cat(3, sta_black_g_shuf_mean, mean(sta_black_g_shuffled, 3));
        sta_white_g_shuf_std  = cat(3, sta_white_g_shuf_std,  std(sta_white_g_shuffled, [], 3));
        sta_black_g_shuf_std  = cat(3, sta_black_g_shuf_std,  std(sta_black_g_shuffled, [], 3));

        % ----- Screens: spike count per stimulus frame -----
        R_s = zeros(1, length(starts_s));
        for i = 1 : length(starts_s) - 1
            R_s(i) = sum(all_spike_times > screens_session.probe_t(starts_s(i)) & ...
                         all_spike_times < screens_session.probe_t(starts_s(i + 1)));
        end
        R_s(end) = sum(all_spike_times > screens_session.probe_t(starts_s(end)) & ...
                       all_spike_times < screens_session.probe_t(starts_s(end)) + 0.25);

        % Screens STA - white
        sta_white_s = bsxfun(@times, stim_white_s, R_s(:)');
        sta_white_s = reshape(mean(sta_white_s, 2), n_y_s, n_x_s);

        sta_white_s_shuffled = nan(n_y_s, n_x_s, 100);
        for sh_i = 1 : 100
            I = randperm(length(R_s));
            st = mean(bsxfun(@times, stim_white_s, R_s(I)), 2);
            sta_white_s_shuffled(:, :, sh_i) = reshape(st, n_y_s, n_x_s);
        end

        % Screens STA - black
        sta_black_s = bsxfun(@times, stim_black_s, R_s(:)');
        sta_black_s = reshape(mean(sta_black_s, 2), n_y_s, n_x_s);

        sta_black_s_shuffled = nan(n_y_s, n_x_s, 100);
        for sh_i = 1 : 100
            I = randperm(length(R_s));
            st = mean(bsxfun(@times, stim_black_s, R_s(I)), 2);
            sta_black_s_shuffled(:, :, sh_i) = reshape(st, n_y_s, n_x_s);
        end

        sta_white_s_store     = cat(3, sta_white_s_store,     sta_white_s);
        sta_black_s_store     = cat(3, sta_black_s_store,     sta_black_s);
        sta_white_s_shuf_mean = cat(3, sta_white_s_shuf_mean, mean(sta_white_s_shuffled, 3));
        sta_black_s_shuf_mean = cat(3, sta_black_s_shuf_mean, mean(sta_black_s_shuffled, 3));
        sta_white_s_shuf_std  = cat(3, sta_white_s_shuf_std,  std(sta_white_s_shuffled, [], 3));
        sta_black_s_shuf_std  = cat(3, sta_black_s_shuf_std,  std(sta_black_s_shuffled, [], 3));

    end  % cluster loop

    % =====================================================================
    % Figures: Goggles (2x2 layout - no warp panel)
    % =====================================================================
    goggles_fnames = {};

    for i = 1 : size(sta_white_g_store, 3)

        c = clusters(i);

        rf_white_g         = sta_white_g_store(:, :, i);
        subtracted_white_g = rf_white_g - sta_white_g_shuf_mean(:, :, i);
        rf_black_g         = sta_black_g_store(:, :, i);
        subtracted_black_g = rf_black_g - sta_black_g_shuf_mean(:, :, i);

        max_response = max(max([subtracted_white_g subtracted_black_g]));
        min_response = min(min([subtracted_white_g subtracted_black_g]));
        if max_response == min_response
            max_response = min_response + 0.1;
        end

        figure;

        % White, raw
        subplot(2, 2, 1);
        imagesc(subtracted_white_g);
        colormap(gca, cmap_on);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('White, raw')

        % White RF detection and size measurement
        RF_raw_g       = rf_white_g > sta_white_g_shuf_mean(:, :, i) + threshold * sta_white_g_shuf_std(:, :, i);
        RF_g           = bwareaopen(RF_raw_g, 4);
        [len_w, len_h] = count_squares_for_RF(RF_g);
        size_RF_w      = len_w * degs_per_square
        size_RF_h      = len_h * degs_per_square
        if len_w > 0, RF_w_goggles(end + 1) = size_RF_w; end
        if len_h > 0, RF_h_goggles(end + 1) = size_RF_h; end

        % White, processed
        RF_g = imgaussfilt(double(imresize(RF_g, [h_pix_g + 1, w_pix_g + 1])), 10);
        subplot(2, 2, 2);
        imagesc(RF_g);
        colormap(gca, cmap_on);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
        title('White, proc')

        % Black, raw
        subplot(2, 2, 3);
        imagesc(subtracted_black_g);
        colormap(gca, cmap_off);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('Black, raw')

        % Black RF detection and size measurement
        RF_raw_g       = rf_black_g > sta_black_g_shuf_mean(:, :, i) + threshold * sta_black_g_shuf_std(:, :, i);
        RF_g           = bwareaopen(RF_raw_g, 4);
        [len_w, len_h] = count_squares_for_RF(RF_g);
        size_RF_w      = len_w * degs_per_square
        size_RF_h      = len_h * degs_per_square
        if len_w > 0, RF_w_goggles(end + 1) = size_RF_w; end
        if len_h > 0, RF_h_goggles(end + 1) = size_RF_h; end

        % Black, processed
        RF_g = imgaussfilt(double(imresize(RF_g, [h_pix_g + 1, w_pix_g + 1])), 10);
        subplot(2, 2, 4);
        imagesc(RF_g);
        colormap(gca, cmap_off);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
        title('Black, proc')

        FigureTitle(gcf, sprintf('Cluster %i (%s) - Goggles', c.id, c.region_str));
        goggles_fnames{i} = fullfile(goggles_save_dir, ...
            sprintf('%s_cluster_%i_rf_goggles.pdf', probe_ids{probe_i}, c.id));
        print(goggles_fnames{i}, '-dpdf');
        close(gcf);

    end  % goggles cluster figure loop

    new_fname_g = fullfile(goggles_save_dir, ...
        sprintf('%s_rf_by_cluster_goggles.pdf', probe_ids{probe_i}));
    join_pdfs(goggles_fnames, new_fname_g, true);

    % =====================================================================
    % Figures: Screens (2x3 layout - includes warped-to-screen panel)
    % =====================================================================
    screens_fnames = {};

    for i = 1 : size(sta_white_s_store, 3)

        c = clusters(i);

        rf_white_s         = sta_white_s_store(:, :, i);
        subtracted_white_s = rf_white_s - sta_white_s_shuf_mean(:, :, i);
        rf_black_s         = sta_black_s_store(:, :, i);
        subtracted_black_s = rf_black_s - sta_black_s_shuf_mean(:, :, i);

        max_response = max(max([subtracted_white_s subtracted_black_s]));
        min_response = min(min([subtracted_white_s subtracted_black_s]));
        if max_response == min_response
            max_response = min_response + 0.1;
        end

        figure;

        % White, raw
        subplot(2, 3, 1);
        imagesc(subtracted_white_s);
        colormap(gca, cmap_on);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('White, raw')

        % White RF detection and size measurement
        RF_raw_s       = rf_white_s > sta_white_s_shuf_mean(:, :, i) + threshold * sta_white_s_shuf_std(:, :, i);
        RF_s           = bwareaopen(RF_raw_s, 4);
        [len_w, len_h] = count_squares_for_RF(RF_s);
        size_RF_w      = len_w * degs_per_square
        size_RF_h      = len_h * degs_per_square
        if len_w > 0, RF_w_screens(end + 1) = size_RF_w; end
        if len_h > 0, RF_h_screens(end + 1) = size_RF_h; end

        % White, processed
        RF_s = imgaussfilt(double(imresize(RF_s, [h_pix_s + 1, w_pix_s + 1])), 10);
        subplot(2, 3, 2);
        imagesc(RF_s);
        colormap(gca, cmap_on);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_s, n_y_s, 1]);
        title('White, proc')

        % White, warped to screen
        RF_s_warped = imresize(RF_s, [h_pix_s + 1, w_pix_s + 1], 'nearest');
        warped_rf_s = RF_s_warped(I_s);
        subplot(2, 3, 3);
        imagesc(warped_rf_s);
        colormap(gca, cmap_on);
        set(gca, 'plotboxaspectratio', [30, 18, 1], 'xtick', [], 'ytick', [], 'xdir', 'reverse');
        title('White, on screen')

        % Black, raw
        subplot(2, 3, 4);
        imagesc(subtracted_black_s);
        colormap(gca, cmap_off);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('Black, raw')

        % Black RF detection and size measurement
        RF_raw_s       = rf_black_s > sta_black_s_shuf_mean(:, :, i) + threshold * sta_black_s_shuf_std(:, :, i);
        RF_s           = bwareaopen(RF_raw_s, 4);
        [len_w, len_h] = count_squares_for_RF(RF_s);
        size_RF_w      = len_w * degs_per_square
        size_RF_h      = len_h * degs_per_square
        if len_w > 0, RF_w_screens(end + 1) = size_RF_w; end
        if len_h > 0, RF_h_screens(end + 1) = size_RF_h; end

        % Black, processed
        RF_s = imgaussfilt(double(imresize(RF_s, [h_pix_s + 1, w_pix_s + 1])), 10);
        subplot(2, 3, 5);
        imagesc(RF_s);
        colormap(gca, cmap_off);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_s, n_y_s, 1]);
        title('Black, proc')

        % Black, warped to screen
        RF_s_warped = imresize(RF_s, [h_pix_s + 1, w_pix_s + 1], 'nearest');
        warped_rf_s = RF_s_warped(I_s);
        subplot(2, 3, 6);
        imagesc(warped_rf_s);
        colormap(gca, cmap_off);
        set(gca, 'plotboxaspectratio', [30, 18, 1], 'xtick', [], 'ytick', [], 'xdir', 'reverse');
        title('Black, on screen')

        FigureTitle(gcf, sprintf('Cluster %i (%s) - Screens', c.id, c.region_str));
        screens_fnames{i} = fullfile(screens_save_dir, ...
            sprintf('%s_cluster_%i_rf_screens.pdf', probe_ids{probe_i}, c.id));
        print(screens_fnames{i}, '-dpdf');
        close(gcf);

    end  % screens cluster figure loop

    new_fname_s = fullfile(screens_save_dir, ...
        sprintf('%s_rf_by_cluster_screens.pdf', probe_ids{probe_i}));
    join_pdfs(screens_fnames, new_fname_s, true);

end  % probe loop

% =========================================================================
% Display comparison statistics
% =========================================================================
fprintf('\n--- Goggles RF sizes (degrees) ---\n');
max_RF_w_goggles  = max(RF_w_goggles)
mean_RF_w_goggles = mean(RF_w_goggles)
min_RF_w_goggles  = min(RF_w_goggles)
max_RF_h_goggles  = max(RF_h_goggles)
mean_RF_h_goggles = mean(RF_h_goggles)
min_RF_h_goggles  = min(RF_h_goggles)

fprintf('\n--- Screens RF sizes (degrees) ---\n');
max_RF_w_screens  = max(RF_w_screens)
mean_RF_w_screens = mean(RF_w_screens)
min_RF_w_screens  = min(RF_w_screens)
max_RF_h_screens  = max(RF_h_screens)
mean_RF_h_screens = mean(RF_h_screens)
min_RF_h_screens  = min(RF_h_screens)
