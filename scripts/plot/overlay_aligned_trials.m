% Plot velocity traces for replay and original trials, where the replay has
% been aligned to the original.
%
%   Specify options:
%
%       experiment_groups:      Will generate plots for all replay trials for all probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21'
%                               Should be a cell array of strings with each
%                               entry an experiment group
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
%       n_traces_per_fig:       number of traces to plot on one A4 page
%
% If `save_figs` is true, then one .pdf will be saved for each probe
% recording, and contain traces for all replay trials.


experiment_groups       = {'mismatch_darkness_oct21'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'alignment'};
n_traces_per_fig        = 5;


%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);


for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    trials      = data.motion_trials();
    trial_count = 0;
    
    for jj = 1 : length(trials)
        
        if ~trials{jj}.is_replay
            continue
        end
        
        trial_count = trial_count + 1;
        axis_n = mod(trial_count-1, n_traces_per_fig) + 1;
        
        if axis_n == 1
            % start a new figure
            h_fig = ctl.figs.a4figure();
        end
        
        h_ax                    = subplot(n_traces_per_fig, 1, axis_n); hold on;
        
        replay_trial_speed      = trials{jj}.multiplexer_speed;
        
        aligned_trial           = trials{jj}.to_aligned;
        original_trial_speed    = aligned_trial.original_trial.treadmill_speed;
        overlay                 = nan(size(replay_trial_speed));
        n_points                = min(length(overlay) - aligned_trial.offset, length(original_trial_speed));
        overlay(aligned_trial.offset + (1:n_points)) = original_trial_speed(1:n_points);
        
        plot(h_ax, trials{jj}.probe_t, replay_trial_speed);
        plot(h_ax, trials{jj}.probe_t, overlay);
        
        title(sprintf('Trial %i, replay of trial %i', trials{jj}.trial_id, trials{jj}.original_trial_id));
        
        set(h_ax, 'plotboxaspectratio', [5, 1, 1], 'box', 'off');
        
        if axis_n == n_traces_per_fig
            ctl.figs.save_fig_to_join();
        end
    end
    
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end
