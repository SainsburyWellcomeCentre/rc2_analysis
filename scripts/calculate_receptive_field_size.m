close all

ctl = RC2Analysis();
probe_ids  = ctl.get_probe_ids('visual_flow', ...
                               'mismatch_nov20', ...
                               'mismatch_jul21', ...
                               'mismatch_darkness_oct21', ...
                               'darkness');
RF_w = [];
RF_h = [];
                           
for probe_i = 1 : length(probe_ids)
    this_probe = probe_ids(probe_i)
    if strcmp(this_probe, 'CAA-1113219_rec1_rec2_rec3')
        % skip this probe as the recording started in the middle of the
        % sparse noise session
        continue
    end
    
    if strcmp(probe_ids{probe_i}, 'CA_176_1_rec1_rec2_rec3') || strcmp(probe_ids{probe_i}, 'CA_176_3_rec1_rec2_rec3')
        sparse_noise_info_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_warped_rc2_20200302.mat';
        warp_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\warp\warp_sony_projector_300x180.mat';
    elseif strcmp(this_probe, 'CAA-1114977_rec1_rec2_rec3') || strcmp(this_probe, 'CAA-1114978_rec1_rec2') || ...
        strcmp(this_probe, 'CAA-1114979_rec1_rec2') || strcmp(this_probe, 'CAA-1114980_rec1_rec2') || ...
        strcmp(this_probe, 'CAA-1115689_rec1_rec2') || strcmp(this_probe, 'CAA-1115691_rec1_rec2') 
        sparse_noise_info_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_warped_mp_300_20210827.mat';
        warp_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\warp\warp_mp_300.mat';
    else
        sparse_noise_info_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\user\mvelez\sparse_noise\sparse_noise_warped_rc2_20200708.mat';
        warp_fname = 'C:\Users\lee\Documents\mvelez\visual_stimuli\warp\warp_sony_projector_300x180.mat';
    end
    
    [stim_mtx, h_pix, w_pix] = stimulus_matrix(sparse_noise_info_fname);
    n_y = size(stim_mtx, 1);
    n_x = size(stim_mtx, 2);
    stim_mtx = reshape(stim_mtx, [], size(stim_mtx, 3));
    stim_mtx_white = stim_mtx == 1;
    stim_mtx_black = stim_mtx == 0;

    
    data = ctl.load_formatted_data(this_probe{1});
    
    if strcmp(this_probe, 'CAA-1112872_rec1_rec1b_rec2_rec3') || strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3') || strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3')
        sparse_noise_session = data.sessions{1, 3};
    else
        sparse_noise_session = data.sessions{1, 2};
    end
    
    animal_id = split(this_probe, "_");
    if strcmp(this_probe, 'CA_176_1_rec1_rec2_rec3') || strcmp(this_probe, 'CA_176_3_rec1_rec2_rec3')
        animal_id = join([animal_id(1), animal_id(2), animal_id(3)], '_');
    else
        animal_id = animal_id(1);
    end
    starts = pd_times(animal_id, sparse_noise_session);
       
    sta_white_store = [];
    sta_black_store = [];
    
    sta_white_shuffle_mean_store = [];
    sta_black_shuffle_mean_store = [];
    sta_white_shuffle_std_store = [];
    sta_black_shuffle_std_store = [];
    
    % Collect spikes for all clusters and sort them by time - not sure why
    % we pool all clusters together
    clusters = data.VISp_clusters();
    for clust_i = 1 : length(clusters)
        all_spike_times = [];
        all_spike_times = [all_spike_times; clusters(clust_i).spike_times];

        all_spike_times = sort(all_spike_times);

        % Identify total spikes in the time interval of interest
        R = [];
        for i = 1 : length(starts)-1
            R(i) = sum(all_spike_times > sparse_noise_session.probe_t(starts(i)) & ...
                all_spike_times < sparse_noise_session.probe_t(starts(i+1)));
        end
        R(end+1) = sum(all_spike_times > sparse_noise_session.probe_t(starts(end)) & ...
            all_spike_times < sparse_noise_session.probe_t(starts(end))+0.25);

        % Calculate spike triggered average (STA)
        sta_white = bsxfun(@times, stim_mtx_white, R(:)');
        sta_white = mean(sta_white, 2);
        sta_white = reshape(sta_white, n_y, n_x);

        sta_white_shuffled = nan(n_y, n_x, 100);
        tic;
        for sh_i = 1 : 100
            I = randperm(length(R));
            st = bsxfun(@times, stim_mtx_white, R(I));
            st = mean(st, 2);
            sta_white_shuffled(:, :, sh_i) = reshape(st, n_y, n_x);
        end
        toc;

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

        sta_white_store = cat(3, sta_white_store, sta_white);
        sta_black_store = cat(3, sta_black_store, sta_black);

        sta_white_shuffle_mean_store = cat(3, sta_white_shuffle_mean_store, mean(sta_white_shuffled, 3));
        sta_black_shuffle_mean_store = cat(3, sta_black_shuffle_mean_store, mean(sta_black_shuffled, 3));
        sta_white_shuffle_std_store = cat(3, sta_white_shuffle_std_store, std(sta_white_shuffled, [], 3));
        sta_black_shuffle_std_store = cat(3, sta_black_shuffle_std_store, std(sta_black_shuffled, [], 3));
    end

    
    %% account for warp
    load(warp_fname, 'scal');
    
    [X, Y] = meshgrid(1:size(scal.vcoords, 2), 1:size(scal.vcoords, 1));
    [Xq, Yq] = meshgrid(1: w_pix + 1, 1: h_pix + 1);
    orixi = round(interp2(scal.vcoords(:, :, 1), scal.vcoords(:, :, 2), scal.tcoords(:, :, 1), Xq, Yq));
    oriyi = round(interp2(scal.vcoords(:, :, 1), scal.vcoords(:, :, 2), scal.tcoords(:, :, 2), Xq, Yq));
    I = sub2ind(size(orixi), oriyi, orixi);
    
    screen_white_bw = [];
    screen_black_bw = [];
    
    save_fnames = {};
    
    for i = 1 : size(sta_white_store, 3)
        
        c = clusters(i);

        degs_per_square = 5;
        
        rf_white = sta_white_store(:, :, i);
        subtracted_white = rf_white - sta_white_shuffle_mean_store(:, :, i);
        rf_black = sta_black_store(:, :, i);
        subtracted_black = rf_black - sta_black_shuffle_mean_store(:, :, i);
        
        max_response = max(max([subtracted_white subtracted_black]));
        min_response = min(min([subtracted_white subtracted_black]));
        
        if max_response == min_response
            max_response = min_response + 0.1;
        end
        
        figure;
        subplot(2, 3, 1);
        imagesc(subtracted_white);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('White, raw')
        
        threshold = 4;
        
        t = rf_white > sta_white_shuffle_mean_store(:, :, i) + threshold * sta_white_shuffle_std_store(:, :, i);
        t2 = bwareaopen(t, 4);
        [len_w, len_h] = count_squares_for_RF(t2);
        size_RF_w = len_w * degs_per_square
        size_RF_h = len_h * degs_per_square
        
        if len_w > 0
            RF_w(end+1) = size_RF_w;
        end
        if len_h > 0 
            RF_h(end+1) = size_RF_h;
        end
        
        t2 = imresize(t2, [h_pix + 1, w_pix + 1]);
        t2 = imgaussfilt(double(t2), 10);
        
        
        subplot(2, 3, 2);
        imagesc(t2);

        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x, n_y, 1]);
        title('White, proc')
        
        t2 = imresize(t2, [h_pix + 1, w_pix + 1], 'nearest');
        warped_rf = t2(I);
        
        subplot(2, 3, 3);
        imagesc(warped_rf);
        set(gca, 'plotboxaspectratio', [30, 18, 1], 'xtick', [], 'ytick', [], 'xdir', 'reverse');
        title('White, on screen')
        
        if sum(warped_rf(:)) > 1
            screen_white_bw = cat(3, screen_white_bw, warped_rf);
        end
        
        subplot(2, 3, 4);
        imagesc(subtracted_black);
        axis image;
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'CLim', [min_response max_response]);
        title('Black, raw')
        
        
        t = rf_black > sta_black_shuffle_mean_store(:, :, i) + threshold * sta_black_shuffle_std_store(:, :, i);
        t2 = bwareaopen(t, 4);
        [len_w, len_h] = count_squares_for_RF(t2);
        size_RF_w = len_w * degs_per_square
        size_RF_h = len_h * degs_per_square
        
        if len_w > 0
            RF_w(end+1) = size_RF_w;
        end
        if len_h > 0 
            RF_h(end+1) = size_RF_h;
        end
        
        t2 = imresize(t2, [h_pix + 1, w_pix + 1]);
        t2 = imgaussfilt(double(t2), 10);
        
        subplot(2, 3, 5)
        imagesc(t2);
        set(gca, 'xtick', [], 'ytick', [], 'xdir', 'reverse', 'plotboxaspectratio', [n_x, n_y, 1]);
        title('Black, proc')
        
        t2 = imresize(t2, [h_pix + 1, w_pix + 1], 'nearest');
        warped_rf = t2(I);
        
        subplot(2, 3, 6)
        imagesc(warped_rf)
        set(gca, 'plotboxaspectratio', [30, 18, 1], 'xtick', [], 'ytick', [], 'xdir', 'reverse');
        title('Black, on screen')
        
        if sum(warped_rf(:)) > 1
            screen_black_bw = cat(3, screen_black_bw, warped_rf);
        end
        
        FigureTitle(gcf, sprintf('Cluster %i (%s)', c.id, c.region_str));
        
        save_fnames{i} = fullfile('C:\Users\lee\Documents\mvelez\figures\new_RFs', sprintf('%s_cluster_%i_rf.pdf', probe_ids{probe_i}, c.id));
        print(save_fnames{i}, '-dpdf');
    end
    
    new_fname = fullfile('C:\Users\lee\Documents\mvelez\figures\new_RFs', sprintf('%s_rf_by_cluster.pdf', probe_ids{probe_i}));
    join_pdfs(save_fnames, new_fname, true);
    
end  

max_RF_w = max(RF_w)
mean_RF_w = mean(RF_w)
min_RF_w = min(RF_w)

max_RF_h = max(RF_h)
mean_RF_h = mean(RF_h)
min_RF_h = min(RF_h)


