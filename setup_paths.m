% sets up paths for use of the analysis scripts

path_config     = path_config();

addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(path_config.npy_matlab_dir));
addpath(genpath(path_config.pdf_dir));
addpath(genpath(path_config.spikes_dir));
