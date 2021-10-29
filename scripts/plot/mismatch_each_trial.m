% plot information for each trial
experiment_groups       = {'mismatch_jul21'};
trial_group_labels      = {'RVT_gain_up', 'RVT_gain_down', 'RV_gain_up', 'RV_gain_down'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'mismatch_trials'};

limits                  = [-1, 1];
n_traces_per_fig        = 5;

ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);


for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    
    for jj = 1 : length(trial_group_labels)
        
        trials      = data.get_trials_with_trial_group_label(trial_group_labels{ii});
        
        for kk = 1 : length(trials)
            
            axis_n = mod(kk-1, n_traces_per_fig) + 1;
            
            if axis_n == 1
                % start a new figure
                h_fig = ctl.figs.a4figure();
            end
            
            h_ax                    = subplot(n_traces_per_fig, 1, axis_n); hold on;
            
            mm_onset_t              = trials{kk}.mismatch_onset_t();
            mm_offset_t             = trials{kk}.mismatch_offset_t();
            idx                     = data.time_idx_around_trigger_time(trials{kk}, mm_onset_t, limits);
            treadmill_speed         = trials{kk}.treadmill_speed(idx);
            multiplexer_speed       = trials{kk}.multiplexer_speed(idx);
            t                       = trials{kk}.probe_t(idx) - mm_onset_t;
            
            plot(h_ax, t, multiplexer_speed, 'r');
            plot(h_ax, t, treadmill_speed, 'k');
            line(h_ax, [0, 0], get(h_ax, 'ylim'), 'color', 'k', 'linestyle', '--');
            line(h_ax, (mm_offset_t - mm_onset_t)*[1, 1], get(h_ax, 'ylim'), 'color', 'k', 'linestyle', '--');
            
            title(sprintf('Trial %i', trials{kk}.trial_id));
            
            set(h_ax, 'plotboxaspectratio', [5, 1, 1], 'box', 'off');
            
            if axis_n == n_traces_per_fig
                ylabel(h_ax, '(cm/s)');
                xlabel(h_ax, 'Time (s)');
                FigureTitle(h_fig, trial_group_labels{jj});
                ctl.figs.save_fig_to_join();
            end
        end
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end
