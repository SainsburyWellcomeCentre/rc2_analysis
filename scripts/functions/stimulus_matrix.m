function [mtx, h_pix, w_pix]  = stimulus_matrix(stim_fname)
%%STIMULUS_MATRIX Load and compute sparse noise stimulus matrix
%
%   [MTX, H_PIX, W_PIX] = stimulus_matrix(stim_fname) 
%   loads the stimulus file, computes and returns the stimulus matrix
%   for the sparse_noise experimental protocol. It returns also the 
%   shape in number of pixels as H_PIX and W_PIX.

% load sparse_noise stimulus data
stim = load(stim_fname);

n_stim = length(stim.cols);

% initialize matrix array
% stim.grid_size correponds to monitor shape
mtx = 0.5*ones([stim.grid_size, n_stim]);

for stim_i = 1 : n_stim
    % at every time point t
    t = 0.5 * ones(stim.grid_size);
    
    % fill the t matrix with the correct sparse_noise stimulus
    t(stim.locations{stim_i}) = stim.cols{stim_i};

    % include t in the main matrix mtx
    mtx(:, :, stim_i) = t;
end

% Get the stimulus matrix dimension in pixels
h_pix = stim.h_pix;
w_pix = stim.w_pix;
