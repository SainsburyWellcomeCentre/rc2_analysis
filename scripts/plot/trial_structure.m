% plot information for each trial
experiment_groups       = {'mismatch_jul21'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'trial_structure', 'mismatch_jul21'};


ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.VISp_clusters();
    trials      = data.motion_trials();
    
    for jj = 1 : length(trials)
        
        ts = TrialStructure();
        ts.plot(trials{jj}.to_aligned, clusters);
        
        ctl.figs.save_fig_to_join();
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end
