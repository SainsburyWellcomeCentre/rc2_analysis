% Plot velocity/camera/population FR information for each trial
%
%   Specify options:
%
%       experiment_groups:      Will generate plots for all trials
%                               and all probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21'
%                               Should be a cell array of strings with each
%                               entry an experiment group
%
%
%       save_figs:              true or false, whether to save the figures to pdf
%
%       overwrite:              true or false. If figure pdf's already exist,
%                               whether to overwrite 
%       
%       figure_dir:             cell array of strings specifying which
%                               directory to save pdf's. The directory will
%                               be relative to the directory specified by
%                               path_config.figure_dir (in
%                               `path_config.m`), so that {'one', 'two',
%                               'three'} will save .pdfs to:
%                               <path_config.figure_dir>\one\two\three\      
%
% If `save_figs` is true, one pdf will be created for each probe recording,
% and contain a A4 page for each trial, containing the information about
% that trial including:
%
%   R, position, T velocity traces
%   Population FR traces
%   Stationary and motion masks
%   Analysis window mask
%   Solenoid state
%   Mismatch window

experiment_groups       = {'training_running'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'trial_structure', 'training_running'};


%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.selected_clusters();
    trials      = data.motion_trials();
    
    for jj = 1 : length(trials)
        
        ts = TrialStructurePosition();
        ts.plot(trials{jj}.to_aligned, clusters);
        
        ctl.figs.save_fig_to_join();
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end
