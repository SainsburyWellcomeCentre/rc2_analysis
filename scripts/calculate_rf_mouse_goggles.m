% Calculate receptive field (RF) sizes from combined goggles sessions.
%   Sessions rec1 (2500 stimuli) and rec3 (2500 stimuli) are concatenated
%   to compute a single RF from 5000 total presentations.
%
% Both sessions belong to the same probe recording (experiment group:
% sparse_noise_goggles). RFs are computed as spike-triggered averages (STA)
% for white (on) and black (off) sparse noise rectangles.
%
% How are RFs calculated?
% 1. The starting time of every new matrix presentation is found from the
%    photodiode signal for both rec1 and rec3.
% 2. We calculate the total spike count R per cluster at every stimulus
%    presentation in both sessions.
% 3. Spike counts from rec1 and rec3 are concatenated (5000 total).
% 4. The STA is calculated for white and black rectangles using combined data.
% 5. A mean and std shuffled STA is generated (100 shuffles).
% 6. RF_raw = STA > (STA_shuffled_mean + threshold * STA_shuffled_std)
%    Connected squares are retained with bwareaopen.
% 7. Width and height of the RF are measured and converted to degrees.
%
% Output PDFs:
%   Combined (2x2 layout):
%       C:\Users\lee\Documents\mvelez\figures\new_RFs\goggles_combined\
%   Comparison (2x4 layout, if enabled):
%       C:\Users\lee\Documents\mvelez\figures\new_RFs\goggles_comparison\

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
degs_per_square = 5;     % degrees per square in the stimulus (constant)
threshold       = 2.5;   % multiplicative factor applied to shuffled STA std

% Optional: Generate comparison plots showing rec1 vs rec3 separately
show_individual_sessions = true;  % Set to false to only generate combined RFs

% ON/OFF colormaps: white→red (responses to white/ON squares)
%                  white→blue (responses to black/OFF squares)
cmap_on  = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)'];  % white→red
cmap_off = [linspace(1,0,256)', linspace(1,0,256)', linspace(1,1,256)'];  % white→blue

% -------------------------------------------------------------------------
% Stimulus file path
% -------------------------------------------------------------------------
goggles_stim_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_wisecoco_2500stim.mat';

% -------------------------------------------------------------------------
% Output directories
% -------------------------------------------------------------------------
goggles_combined_dir = 'C:\Users\lee\Documents\mvelez\figures\new_RFs\goggles_combined';
goggles_comparison_dir = 'C:\Users\lee\Documents\mvelez\figures\new_RFs\goggles_comparison';
if ~exist(goggles_combined_dir, 'dir'), mkdir(goggles_combined_dir); end
if show_individual_sessions && ~exist(goggles_comparison_dir, 'dir')
    mkdir(goggles_comparison_dir);
end

% -------------------------------------------------------------------------
% Select recordings
% -------------------------------------------------------------------------
ctl       = RC2Analysis();
probe_ids = ctl.get_probe_ids('sparse_noise_goggles');

% Subset selection: Leave empty to process all probes, or specify probe IDs
probe_ids_to_analyze = {'CAA-1124370_rec1_rec2_rec3', 'CAA-1124371_rec1_rec2_rec3'};  % Empty = all probes

% Apply probe filtering if specified
if ~isempty(probe_ids_to_analyze)
    probe_ids = intersect(probe_ids, probe_ids_to_analyze, 'stable');
    fprintf('Restricting analysis to %d specified probes\n', length(probe_ids));
end

% -------------------------------------------------------------------------
% Storage for RF dimensions across all probes
% -------------------------------------------------------------------------
RF_w_combined = [];
RF_h_combined = [];

% -------------------------------------------------------------------------
% Pre-compute stimulus matrices (same for all probes)
% -------------------------------------------------------------------------

% --- Goggles stimulus (2500 stimuli for rec1 and rec3) ---
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

% For combined analysis: duplicate stimulus matrix (5000 presentations)
stim_white_combined = [stim_white_g, stim_white_g];  % (n_y*n_x) x 5000
stim_black_combined = [stim_black_g, stim_black_g];  % (n_y*n_x) x 5000

fprintf('Loaded stimulus with %d presentations (%d x 2 sessions = %d combined)\n', ...
    n_stim_g, n_stim_g, 2*n_stim_g);

% =========================================================================
% Main loop over probe recordings
% =========================================================================
for probe_i = 1 : length(probe_ids)

    this_probe = probe_ids(probe_i)

    % ---------------------------------------------------------------------
    % Load data and select sessions
    % ---------------------------------------------------------------------
    data         = ctl.load_formatted_data(this_probe{1});
    session_rec1 = data.sessions{1, 1};   % rec1 - goggles (2500 stimuli)
    session_rec3 = data.sessions{1, 3};   % rec3 - goggles (2500 stimuli)

    % Get animal ID for photodiode parameter lookup
    animal_id = get_animal_id(this_probe);

    % Find stimulus onset times from photodiode signal
    starts_rec1 = photodiode_times(animal_id, session_rec1, "sparse_noise_goggles");
    starts_rec3 = photodiode_times(animal_id, session_rec3, "sparse_noise_goggles");

    fprintf('  rec1: %d stimulus onsets detected\n', length(starts_rec1));
    fprintf('  rec3: %d stimulus onsets detected\n', length(starts_rec3));

    % ---------------------------------------------------------------------
    % Initialise STA storage for this probe
    % ---------------------------------------------------------------------
    sta_white_combined_store     = [];   sta_black_combined_store     = [];
    sta_white_combined_shuf_mean = [];   sta_black_combined_shuf_mean = [];
    sta_white_combined_shuf_std  = [];   sta_black_combined_shuf_std  = [];

    % Optional: individual session storage for comparison
    if show_individual_sessions
        sta_white_rec1_store     = [];   sta_black_rec1_store     = [];
        sta_white_rec1_shuf_mean = [];   sta_black_rec1_shuf_mean = [];
        sta_white_rec1_shuf_std  = [];   sta_black_rec1_shuf_std  = [];
        
        sta_white_rec3_store     = [];   sta_black_rec3_store     = [];
        sta_white_rec3_shuf_mean = [];   sta_black_rec3_shuf_mean = [];
        sta_white_rec3_shuf_std  = [];   sta_black_rec3_shuf_std  = [];
    end

    % ---------------------------------------------------------------------
    % Compute STAs for every VISp cluster
    % NOTE: restrict_to_visp=false includes clusters without anatomy data
    % (labelled 'unknownLocation'). Update to true once anatomy is available.
    % ---------------------------------------------------------------------
    clusters = data.VISp_clusters(true, 'any', false);

    for clust_i = 1 : length(clusters)

        all_spike_times = sort(clusters(clust_i).spike_times);

        % ----- rec1: spike count per stimulus frame -----
        R_rec1 = zeros(1, length(starts_rec1));
        for i = 1 : length(starts_rec1) - 1
            R_rec1(i) = sum(all_spike_times > session_rec1.probe_t(starts_rec1(i)) & ...
                            all_spike_times < session_rec1.probe_t(starts_rec1(i + 1)));
        end
        R_rec1(end) = sum(all_spike_times > session_rec1.probe_t(starts_rec1(end)) & ...
                          all_spike_times < session_rec1.probe_t(starts_rec1(end)) + 0.25);

        % ----- rec3: spike count per stimulus frame -----
        R_rec3 = zeros(1, length(starts_rec3));
        for i = 1 : length(starts_rec3) - 1
            R_rec3(i) = sum(all_spike_times > session_rec3.probe_t(starts_rec3(i)) & ...
                            all_spike_times < session_rec3.probe_t(starts_rec3(i + 1)));
        end
        R_rec3(end) = sum(all_spike_times > session_rec3.probe_t(starts_rec3(end)) & ...
                          all_spike_times < session_rec3.probe_t(starts_rec3(end)) + 0.25);

        % ----- Combined: concatenate spike counts from rec1 and rec3 -----
        R_combined = [R_rec1, R_rec3];  % 1 x 5000

        % Combined STA - white
        sta_white_comb = bsxfun(@times, stim_white_combined, R_combined(:)');
        sta_white_comb = reshape(mean(sta_white_comb, 2), n_y_g, n_x_g);

        sta_white_comb_shuffled = nan(n_y_g, n_x_g, 100);
        for sh_i = 1 : 100
            I = randperm(length(R_combined));
            st = mean(bsxfun(@times, stim_white_combined, R_combined(I)), 2);
            sta_white_comb_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
        end

        % Combined STA - black
        sta_black_comb = bsxfun(@times, stim_black_combined, R_combined(:)');
        sta_black_comb = reshape(mean(sta_black_comb, 2), n_y_g, n_x_g);

        sta_black_comb_shuffled = nan(n_y_g, n_x_g, 100);
        for sh_i = 1 : 100
            I = randperm(length(R_combined));
            st = mean(bsxfun(@times, stim_black_combined, R_combined(I)), 2);
            sta_black_comb_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
        end

        sta_white_combined_store     = cat(3, sta_white_combined_store,     sta_white_comb);
        sta_black_combined_store     = cat(3, sta_black_combined_store,     sta_black_comb);
        sta_white_combined_shuf_mean = cat(3, sta_white_combined_shuf_mean, mean(sta_white_comb_shuffled, 3));
        sta_black_combined_shuf_mean = cat(3, sta_black_combined_shuf_mean, mean(sta_black_comb_shuffled, 3));
        sta_white_combined_shuf_std  = cat(3, sta_white_combined_shuf_std,  std(sta_white_comb_shuffled, [], 3));
        sta_black_combined_shuf_std  = cat(3, sta_black_combined_shuf_std,  std(sta_black_comb_shuffled, [], 3));

        % ----- Optional: Individual session STAs for comparison -----
        if show_individual_sessions
            % rec1 STA - white
            sta_white_r1 = bsxfun(@times, stim_white_g, R_rec1(:)');
            sta_white_r1 = reshape(mean(sta_white_r1, 2), n_y_g, n_x_g);

            sta_white_r1_shuffled = nan(n_y_g, n_x_g, 100);
            for sh_i = 1 : 100
                I = randperm(length(R_rec1));
                st = mean(bsxfun(@times, stim_white_g, R_rec1(I)), 2);
                sta_white_r1_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
            end

            % rec1 STA - black
            sta_black_r1 = bsxfun(@times, stim_black_g, R_rec1(:)');
            sta_black_r1 = reshape(mean(sta_black_r1, 2), n_y_g, n_x_g);

            sta_black_r1_shuffled = nan(n_y_g, n_x_g, 100);
            for sh_i = 1 : 100
                I = randperm(length(R_rec1));
                st = mean(bsxfun(@times, stim_black_g, R_rec1(I)), 2);
                sta_black_r1_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
            end

            sta_white_rec1_store     = cat(3, sta_white_rec1_store,     sta_white_r1);
            sta_black_rec1_store     = cat(3, sta_black_rec1_store,     sta_black_r1);
            sta_white_rec1_shuf_mean = cat(3, sta_white_rec1_shuf_mean, mean(sta_white_r1_shuffled, 3));
            sta_black_rec1_shuf_mean = cat(3, sta_black_rec1_shuf_mean, mean(sta_black_r1_shuffled, 3));
            sta_white_rec1_shuf_std  = cat(3, sta_white_rec1_shuf_std,  std(sta_white_r1_shuffled, [], 3));
            sta_black_rec1_shuf_std  = cat(3, sta_black_rec1_shuf_std,  std(sta_black_r1_shuffled, [], 3));

            % rec3 STA - white
            sta_white_r3 = bsxfun(@times, stim_white_g, R_rec3(:)');
            sta_white_r3 = reshape(mean(sta_white_r3, 2), n_y_g, n_x_g);

            sta_white_r3_shuffled = nan(n_y_g, n_x_g, 100);
            for sh_i = 1 : 100
                I = randperm(length(R_rec3));
                st = mean(bsxfun(@times, stim_white_g, R_rec3(I)), 2);
                sta_white_r3_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
            end

            % rec3 STA - black
            sta_black_r3 = bsxfun(@times, stim_black_g, R_rec3(:)');
            sta_black_r3 = reshape(mean(sta_black_r3, 2), n_y_g, n_x_g);

            sta_black_r3_shuffled = nan(n_y_g, n_x_g, 100);
            for sh_i = 1 : 100
                I = randperm(length(R_rec3));
                st = mean(bsxfun(@times, stim_black_g, R_rec3(I)), 2);
                sta_black_r3_shuffled(:, :, sh_i) = reshape(st, n_y_g, n_x_g);
            end

            sta_white_rec3_store     = cat(3, sta_white_rec3_store,     sta_white_r3);
            sta_black_rec3_store     = cat(3, sta_black_rec3_store,     sta_black_r3);
            sta_white_rec3_shuf_mean = cat(3, sta_white_rec3_shuf_mean, mean(sta_white_r3_shuffled, 3));
            sta_black_rec3_shuf_mean = cat(3, sta_black_rec3_shuf_mean, mean(sta_black_r3_shuffled, 3));
            sta_white_rec3_shuf_std  = cat(3, sta_white_rec3_shuf_std,  std(sta_white_r3_shuffled, [], 3));
            sta_black_rec3_shuf_std  = cat(3, sta_black_rec3_shuf_std,  std(sta_black_r3_shuffled, [], 3));
        end

    end  % cluster loop

    % =====================================================================
    % Figures: Combined rec1+rec3 (2x2 layout)
    % =====================================================================
    combined_fnames = {};

    for i = 1 : size(sta_white_combined_store, 3)

        c = clusters(i);

        rf_white_comb         = sta_white_combined_store(:, :, i);
        subtracted_white_comb = rf_white_comb - sta_white_combined_shuf_mean(:, :, i);
        rf_black_comb         = sta_black_combined_store(:, :, i);
        subtracted_black_comb = rf_black_comb - sta_black_combined_shuf_mean(:, :, i);

        max_response = max(max([subtracted_white_comb subtracted_black_comb]));
        min_response = min(min([subtracted_white_comb subtracted_black_comb]));
        if max_response == min_response
            max_response = min_response + 0.1;
        end

        figure;

        % White, raw
        subplot(2, 2, 1);
        imagesc(subtracted_white_comb);
        colormap(gca, cmap_on);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('White, raw')

        % White RF detection and size measurement
        RF_raw_comb    = rf_white_comb > sta_white_combined_shuf_mean(:, :, i) + threshold * sta_white_combined_shuf_std(:, :, i);
        RF_comb        = bwareaopen(RF_raw_comb, 4);
        [len_w, len_h] = count_squares_for_RF(RF_comb);
        size_RF_w      = len_w * degs_per_square;
        size_RF_h      = len_h * degs_per_square;
        if len_w > 0, RF_w_combined(end + 1) = size_RF_w; end
        if len_h > 0, RF_h_combined(end + 1) = size_RF_h; end

        % White, processed
        RF_comb = imgaussfilt(double(imresize(RF_comb, [h_pix_g + 1, w_pix_g + 1])), 10);
        subplot(2, 2, 2);
        imagesc(RF_comb);
        colormap(gca, cmap_on);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
        title('White, proc')

        % Black, raw
        subplot(2, 2, 3);
        imagesc(subtracted_black_comb);
        colormap(gca, cmap_off);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('Black, raw')

        % Black RF detection and size measurement
        RF_raw_comb    = rf_black_comb > sta_black_combined_shuf_mean(:, :, i) + threshold * sta_black_combined_shuf_std(:, :, i);
        RF_comb        = bwareaopen(RF_raw_comb, 4);
        [len_w, len_h] = count_squares_for_RF(RF_comb);
        size_RF_w      = len_w * degs_per_square;
        size_RF_h      = len_h * degs_per_square;
        if len_w > 0, RF_w_combined(end + 1) = size_RF_w; end
        if len_h > 0, RF_h_combined(end + 1) = size_RF_h; end

        % Black, processed
        RF_comb = imgaussfilt(double(imresize(RF_comb, [h_pix_g + 1, w_pix_g + 1])), 10);
        subplot(2, 2, 4);
        imagesc(RF_comb);
        colormap(gca, cmap_off);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
        title('Black, proc')

        FigureTitle(gcf, sprintf('Cluster %i (%s) - Combined rec1+rec3', c.id, c.region_str));
        combined_fnames{i} = fullfile(goggles_combined_dir, ...
            sprintf('%s_cluster_%i_rf_goggles_combined.pdf', probe_ids{probe_i}, c.id));
        print(combined_fnames{i}, '-dpdf');
        close(gcf);

    end  % combined cluster figure loop

    new_fname_combined = fullfile(goggles_combined_dir, ...
        sprintf('%s_rf_by_cluster_goggles_combined.pdf', probe_ids{probe_i}));
    join_pdfs(combined_fnames, new_fname_combined, true);

    % =====================================================================
    % Figures: Comparison rec1 vs rec3 (2x4 layout, optional)
    % =====================================================================
    if show_individual_sessions
        comparison_fnames = {};

        for i = 1 : size(sta_white_rec1_store, 3)

            c = clusters(i);

            % Prepare data for rec1
            rf_white_r1         = sta_white_rec1_store(:, :, i);
            subtracted_white_r1 = rf_white_r1 - sta_white_rec1_shuf_mean(:, :, i);
            rf_black_r1         = sta_black_rec1_store(:, :, i);
            subtracted_black_r1 = rf_black_r1 - sta_black_rec1_shuf_mean(:, :, i);

            % Prepare data for rec3
            rf_white_r3         = sta_white_rec3_store(:, :, i);
            subtracted_white_r3 = rf_white_r3 - sta_white_rec3_shuf_mean(:, :, i);
            rf_black_r3         = sta_black_rec3_store(:, :, i);
            subtracted_black_r3 = rf_black_r3 - sta_black_rec3_shuf_mean(:, :, i);

            % Common color limits across all panels
            max_response = max(max([subtracted_white_r1 subtracted_black_r1 subtracted_white_r3 subtracted_black_r3]));
            min_response = min(min([subtracted_white_r1 subtracted_black_r1 subtracted_white_r3 subtracted_black_r3]));
            if max_response == min_response
                max_response = min_response + 0.1;
            end

            figure;

            % --- Row 1: White RFs ---
            % rec1 white, raw
            subplot(2, 4, 1);
            imagesc(subtracted_white_r1);
            colormap(gca, cmap_on);
            axis image;
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
            title('rec1 White, raw')

            % rec1 white, processed
            RF_raw_r1 = rf_white_r1 > sta_white_rec1_shuf_mean(:, :, i) + threshold * sta_white_rec1_shuf_std(:, :, i);
            RF_r1     = bwareaopen(RF_raw_r1, 4);
            RF_r1     = imgaussfilt(double(imresize(RF_r1, [h_pix_g + 1, w_pix_g + 1])), 10);
            subplot(2, 4, 2);
            imagesc(RF_r1);
            colormap(gca, cmap_on);
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
            title('rec1 White, proc')

            % rec3 white, raw
            subplot(2, 4, 3);
            imagesc(subtracted_white_r3);
            colormap(gca, cmap_on);
            axis image;
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
            title('rec3 White, raw')

            % rec3 white, processed
            RF_raw_r3 = rf_white_r3 > sta_white_rec3_shuf_mean(:, :, i) + threshold * sta_white_rec3_shuf_std(:, :, i);
            RF_r3     = bwareaopen(RF_raw_r3, 4);
            RF_r3     = imgaussfilt(double(imresize(RF_r3, [h_pix_g + 1, w_pix_g + 1])), 10);
            subplot(2, 4, 4);
            imagesc(RF_r3);
            colormap(gca, cmap_on);
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
            title('rec3 White, proc')

            % --- Row 2: Black RFs ---
            % rec1 black, raw
            subplot(2, 4, 5);
            imagesc(subtracted_black_r1);
            colormap(gca, cmap_off);
            axis image;
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
            title('rec1 Black, raw')

            % rec1 black, processed
            RF_raw_r1 = rf_black_r1 > sta_black_rec1_shuf_mean(:, :, i) + threshold * sta_black_rec1_shuf_std(:, :, i);
            RF_r1     = bwareaopen(RF_raw_r1, 4);
            RF_r1     = imgaussfilt(double(imresize(RF_r1, [h_pix_g + 1, w_pix_g + 1])), 10);
            subplot(2, 4, 6);
            imagesc(RF_r1);
            colormap(gca, cmap_off);
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
            title('rec1 Black, proc')

            % rec3 black, raw
            subplot(2, 4, 7);
            imagesc(subtracted_black_r3);
            colormap(gca, cmap_off);
            axis image;
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
            title('rec3 Black, raw')

            % rec3 black, processed
            RF_raw_r3 = rf_black_r3 > sta_black_rec3_shuf_mean(:, :, i) + threshold * sta_black_rec3_shuf_std(:, :, i);
            RF_r3     = bwareaopen(RF_raw_r3, 4);
            RF_r3     = imgaussfilt(double(imresize(RF_r3, [h_pix_g + 1, w_pix_g + 1])), 10);
            subplot(2, 4, 8);
            imagesc(RF_r3);
            colormap(gca, cmap_off);
            set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x_g, n_y_g, 1]);
            title('rec3 Black, proc')

            FigureTitle(gcf, sprintf('Cluster %i (%s) - rec1 vs rec3 Comparison', c.id, c.region_str));
            comparison_fnames{i} = fullfile(goggles_comparison_dir, ...
                sprintf('%s_cluster_%i_rf_goggles_comparison.pdf', probe_ids{probe_i}, c.id));
            print(comparison_fnames{i}, '-dpdf');
            close(gcf);

        end  % comparison cluster figure loop

        new_fname_comparison = fullfile(goggles_comparison_dir, ...
            sprintf('%s_rf_by_cluster_goggles_comparison.pdf', probe_ids{probe_i}));
        join_pdfs(comparison_fnames, new_fname_comparison, true);
    end

end  % probe loop

% =========================================================================
% Display statistics
% =========================================================================
fprintf('\n--- Combined (rec1+rec3) RF sizes (degrees) ---\n');
if ~isempty(RF_w_combined)
    max_RF_w_combined  = max(RF_w_combined)
    mean_RF_w_combined = mean(RF_w_combined)
    min_RF_w_combined  = min(RF_w_combined)
else
    fprintf('No RF width measurements available\n');
end

if ~isempty(RF_h_combined)
    max_RF_h_combined  = max(RF_h_combined)
    mean_RF_h_combined = mean(RF_h_combined)
    min_RF_h_combined  = min(RF_h_combined)
else
    fprintf('No RF height measurements available\n');
end
