% Calculate receptive field (RF) size of each cluster
% Receptive fields are calculated by analysing the responses to a sparse noise stimulus,
% i.e. a matrix of flickering squares, black and white for off and on RFs.
%
% How are RFs calculated?
% 1. The starting time of every new matrix presentation (stimulus) is found by reading
% the photodiode signal
% 2. We calculate the total amount of spikes R per cluster at every stimulus presentation
% 3. The spike triggered average (STA) is calculated for white and black rectangles
% 4. A mean shuffled STA (STA_shuffled_mean) matrix is generated for the black and white conditions
% 5. Receptive fields are found where this condition is satisfied:
%    RF_raw = (STA - STA_shuffled_mean) > t_matrix
%    where t_matrix = STA_shuffled_mean * threshold; threshold is defined at the beginning of this script.
%    We then retain only the connected squares with bwareaopen, obtaining RF.
% 6. We finally measure the x and y dimensions of the RF.



% How many degrees in each square in the stimulus
degs_per_square = 5; % constant, from the setup
threshold       = 4; % multiplicative factor to be applied to the shuffled STA

% SelecRF_rawall datasets for which there is a Sparse Noise session
ctl       = RC2Analysis();
probe_ids = ctl.get_probe_ids('visual_flow', ...
                              'mismatch_nov20', ...
                              'mismatch_jul21', ...
                              'mismatch_darkness_oct21', ...
                              'darkness');

% initialize variable to hold receptive field width and height
RF_w = [];
RF_h = [];
                           
for probe_i = 1 : length(probe_ids)
    this_probe = probe_ids(probe_i)
    if strcmp(this_probe, 'CAA-1113219_rec1_rec2_rec3')
        % skip this probe as the recording started in the middle of the
        % sparse noise session
        continue
    end
    
    % find the correct warp files for a given experiment
    if strcmp(probe_ids{probe_i}, 'CA_176_1_rec1_rec2_rec3') || strcmp(probe_ids{probe_i}, 'CA_176_3_rec1_rec2_rec3')

        sparse_noise_info_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_warped_rc2_20200302.mat';
        warp_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\warp\warp_sony_projector_300x180.mat';
    elseif strcmp(this_probe, 'CAA-1114977_rec1_rec2_rec3')  || strcmp(this_probe, 'CAA-1114978_rec1_rec2') || ...
           strcmp(this_probe, 'CAA-1114979_rec1_rec2')       || strcmp(this_probe, 'CAA-1114980_rec1_rec2') || ...
           strcmp(this_probe, 'CAA-1115689_rec1_rec2')       || strcmp(this_probe, 'CAA-1115691_rec1_rec2') 

        sparse_noise_info_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_warped_mp_300_20210827.mat';
        warp_fname              = 'C:\Users\lee\Documents\mvelez\visual_stimuli\warp\warp_mp_300.mat';
    else 
        sparse_noise_info_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_warped_rc2_20200708.mat';
        warp_fname              = 'C:\Users\lee\Documents\mvelez\visual_stimuli\warp\warp_sony_projector_300x180.mat';
    end
    
    % get the 3D stimulus matrix (n_y * n_x * t), in which t is stimulus repetitions
    % every stimulus is a grid of rectangles, which shape (n_y, n_x) depends on the monitor size
    % We also get the pixel shape of the monitor, useful for later warping calculations.
    [stim_mtx, h_pix, w_pix] = stimulus_matrix(sparse_noise_info_fname);
    n_y                      = size(stim_mtx, 1);
    n_x                      = size(stim_mtx, 2);
    % reshape to 2D matrx
    % new size: (n_y * n_x) x t
    stim_mtx       = reshape(stim_mtx, [], size(stim_mtx, 3));
    % make new matrices on where and when the squares are white or black.
    stim_mtx_white = stim_mtx == 1;
    stim_mtx_black = stim_mtx == 0;

    % load data containing clusters firing rate
    data = ctl.load_formatted_data(this_probe{1});
    
    % select the Sparse Noise session
    if strcmp(this_probe, 'CAA-1112872_rec1_rec1b_rec2_rec3') || strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3') ||...
       strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3')

        % in some datasets it is in a different position
        sparse_noise_session = data.sessions{1, 3};
    else
        sparse_noise_session = data.sessions{1, 2};
    end
    
    % Get animal ID - useful to obtain custom parameters to find starting times from the photodiode signal
    animal_id = get_animal_id(this_probe)

    % find starting times of the stimuli from photodiode
    starts = photodiode_times(animal_id, sparse_noise_session, "sparse_noise");
    
    % initialize arrays
    sta_white_store = [];
    sta_black_store = [];
    
    sta_white_shuffle_mean_store = [];
    sta_black_shuffle_mean_store = [];
    sta_white_shuffle_std_store  = [];
    sta_black_shuffle_std_store  = [];
    
    % Calculate spike triggered average (STA) for each cluster 
    clusters = data.VISp_clusters();
    for clust_i = 1 : length(clusters)
        % for each cluster get and sort spike_times
        all_spike_times = clusters(clust_i).spike_times;
        all_spike_times = sort(all_spike_times);

        % Identify total count of spikes in time windows between starting times
        R = [];
        for i = 1 : length(starts)-1
            R(i) = sum(all_spike_times > sparse_noise_session.probe_t(starts(i)) & ...
                       all_spike_times < sparse_noise_session.probe_t(starts(i + 1)));
        end
        R(end+1) = sum(all_spike_times > sparse_noise_session.probe_t(starts(end)) & ...
                       all_spike_times < sparse_noise_session.probe_t(starts(end)) + 0.25);

        % Calculate spike triggered average (STA) for every white rectangle
        % and reshape the matrix
        sta_white = bsxfun(@times, stim_mtx_white, R(:)');
        sta_white = mean(sta_white, 2);
        sta_white = reshape(sta_white, n_y, n_x);
        
        % shuffle the data
        sta_white_shuffled = nan(n_y, n_x, 100);
        tic;
        for sh_i = 1 : 100
            I = randperm(length(R));
            st = bsxfun(@times, stim_mtx_white, R(I));
            st = mean(st, 2);
            sta_white_shuffled(:, :, sh_i) = reshape(st, n_y, n_x);
        end
        toc;

        % calculate again STA but for black rectangles
        sta_black = bsxfun(@times, stim_mtx_black, R(:)');
        sta_black = mean(sta_black, 2);
        sta_black = reshape(sta_black, n_y, n_x);

        sta_black_shuffled = nan(n_y, n_x, 100);
        for sh_i = 1 : 100
            I = randperm(length(R));
            st = bsxfun(@times, stim_mtx_black, R(I));
            st = mean(st, 2);
            sta_black_shuffled(:, :, sh_i) = reshape(st, n_y, n_x);
        end

        % store the STA results
        sta_white_store = cat(3, sta_white_store, sta_white);
        sta_black_store = cat(3, sta_black_store, sta_black);

        sta_white_shuffle_mean_store = cat(3, sta_white_shuffle_mean_store, mean(sta_white_shuffled, 3));
        sta_black_shuffle_mean_store = cat(3, sta_black_shuffle_mean_store, mean(sta_black_shuffled, 3));
        sta_white_shuffle_std_store  = cat(3, sta_white_shuffle_std_store,  std(sta_white_shuffled, [], 3));
        sta_black_shuffle_std_store  = cat(3, sta_black_shuffle_std_store,  std(sta_black_shuffled, [], 3));
    end

    
    % Load warping data and get the linear indices based on warping information
    load(warp_fname, 'scal');
    
    [X, Y]   = meshgrid(1:size(scal.vcoords, 2), 1:size(scal.vcoords, 1));
    [Xq, Yq] = meshgrid(1: w_pix + 1, 1: h_pix + 1);
    orixi    = round(interp2(scal.vcoords(:, :, 1), scal.vcoords(:, :, 2), scal.tcoords(:, :, 1), Xq, Yq));
    oriyi    = round(interp2(scal.vcoords(:, :, 1), scal.vcoords(:, :, 2), scal.tcoords(:, :, 2), Xq, Yq));
    I        = sub2ind(size(orixi), oriyi, orixi);
    
    screen_white_bw = [];
    screen_black_bw = [];
    
    save_fnames = {};
    
    for i = 1 : size(sta_white_store, 3)
        % for every cluster
        c = clusters(i);

        % calculate the delta STA from black and white stimulus from the shuffled dataset.
        rf_white         = sta_white_store(:, :, i);
        subtracted_white = rf_white - sta_white_shuffle_mean_store(:, :, i);
        rf_black         = sta_black_store(:, :, i);
        subtracted_black = rf_black - sta_black_shuffle_mean_store(:, :, i);
        
        % get the maximum response to normalize plots
        max_response = max(max([subtracted_white subtracted_black]));
        min_response = min(min([subtracted_white subtracted_black]));
        if max_response == min_response
            max_response = min_response + 0.1;
        end
        

        % Make a figure to visualize the shpe and location of on and off
        % receptive fields and calculate their size.

        % setup figure
        figure;
        subplot(2, 3, 1);
        imagesc(subtracted_white);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('White, raw')
        
        % get RFs from subtracted STA (on)
        % the threshold to be applied to the shuffled STA can be controlled at the beginning of the script
        RF_raw         = rf_white > sta_white_shuffle_mean_store(:, :, i) + threshold * sta_white_shuffle_std_store(:, :, i);
        % retain connected squares
        RF             = bwareaopen(RF_raw, 4);
        % measure shape of the receptive field and convert it to degrees
        [len_w, len_h] = count_squares_for_RF(RF);
        size_RF_w      = len_w * degs_per_square
        size_RF_h      = len_h * degs_per_square
        % store calculated receptive fields shapes.
        if len_w > 0
            RF_w(end+1) = size_RF_w;
        end
        if len_h > 0 
            RF_h(end+1) = size_RF_h;
        end
        
        % plot found RF
        RF = imresize(RF, [h_pix + 1, w_pix + 1]);
        RF = imgaussfilt(double(RF), 10);
        subplot(2, 3, 2);
        imagesc(RF);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x, n_y, 1]);
        title('White, proc')
        
        % plot warped RF
        RF = imresize(RF, [h_pix + 1, w_pix + 1], 'nearest');
        warped_rf = RF(I);
        subplot(2, 3, 3);
        imagesc(warped_rf);
        set(gca, 'plotboxaspectratio', [30, 18, 1], 'xtick', [], 'ytick', [], 'xdir', 'reverse');
        title('White, on screen')
        

        % Again, similar plots and calculations for off RF
        subplot(2, 3, 4);
        imagesc(subtracted_black);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('Black, raw')
        
        % Calculate RF
        RF_raw            = rf_black > sta_black_shuffle_mean_store(:, :, i) + threshold * sta_black_shuffle_std_store(:, :, i);
        RF             = bwareaopen(RF_raw, 4);
        [len_w, len_h] = count_squares_for_RF(RF);
        size_RF_w      = len_w * degs_per_square
        size_RF_h      = len_h * degs_per_square
        
        if len_w > 0
            RF_w(end+1) = size_RF_w;
        end
        if len_h > 0 
            RF_h(end+1) = size_RF_h;
        end
        
        RF = imresize(RF, [h_pix + 1, w_pix + 1]);
        RF = imgaussfilt(double(RF), 10);
        
        subplot(2, 3, 5)
        imagesc(RF);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x, n_y, 1]);
        title('Black, proc')
        
        RF = imresize(RF, [h_pix + 1, w_pix + 1], 'nearest');
        warped_rf = RF(I);
        
        subplot(2, 3, 6)
        imagesc(warped_rf)
        set(gca, 'plotboxaspectratio', [30, 18, 1], 'xtick', [], 'ytick', [], 'xdir', 'reverse');
        title('Black, on screen')
        
        % Save figure
        FigureTitle(gcf, sprintf('Cluster %i (%s)', c.id, c.region_str));
        save_fnames{i} = fullfile('C:\Users\lee\Documents\mvelez\figures\new_RFs', sprintf('%s_cluster_%i_rf.pdf', probe_ids{probe_i}, c.id));
        print(save_fnames{i}, '-dpdf');
    end
    
    new_fname = fullfile('C:\Users\lee\Documents\mvelez\figures\new_RFs', sprintf('%s_rf_by_cluster.pdf', probe_ids{probe_i}));
    join_pdfs(save_fnames, new_fname, true);
    
end  

% Display statistics on the shape of receptive fields
max_RF_w  = max(RF_w)
mean_RF_w = mean(RF_w)
min_RF_w  = min(RF_w)

max_RF_h  = max(RF_h)
mean_RF_h = mean(RF_h)
min_RF_h  = min(RF_h)


