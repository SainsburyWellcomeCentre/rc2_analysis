function [mtx, h_pix, w_pix] = stimulus_matrix(stim_fname)

stim = load(stim_fname);

n_stim = length(stim.cols);

mtx = 0.5*ones([stim.grid_size, n_stim]);

for stim_i = 1 : n_stim
    
    t = 0.5 * ones(stim.grid_size);
    t(stim.locations{stim_i}) = stim.cols{stim_i};
    
    mtx(:, :, stim_i) = t;
end

h_pix = stim.h_pix;
w_pix = stim.w_pix;

