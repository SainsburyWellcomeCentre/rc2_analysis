close all

experiment_groups       = {'darkness', 'mismatch_darkness_oct21'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'trial_structure', 'darkness_vel_acc'};


%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.selected_clusters();
    trials      = data.motion_trials();
    
    for jj = 1 : length(trials)
        
        ts = TrialStructureVelAcc();
        ts.plot(trials{jj}.to_aligned, clusters);
        
        ctl.figs.save_fig_to_join();
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end
