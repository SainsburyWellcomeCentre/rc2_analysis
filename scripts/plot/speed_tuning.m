% Plot speed tuning curve for an experimental group
%  Tuning cuves must have been generated and saved using
%  RC2Analysis.create_tuning_curves
%
%   Specify options:
%
%       experiment_groups:      Will generate plots for all clusters
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
%       trial_group_labels:     Will generate a plots with speed tuning
%                               curves for the trials specified. 
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%                               e.g. {'R', {'T_bank', 'T_RT', 'T_R'}, 'RT'}
%                               will plot three sets of tuning curves, the first
%                               for all 'R' (running) trials, the
%                               second for all trials of any of 
%                               'T_bank', 'T_RT' or 'T_R' type, and the
%                               third for all 'RT'
%                               (running+translation) trials.
%                       Importantly, the provided entries must match the
%                       entries used to create the tuning curve information
%                       (see RC2Analysis.create_tuning_curves, and
%                       FormattedData.create_tuning_curves).
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
% and contain a A4 page for each cluster, containing the raster data for
% several conditions.

%%
% experiment_groups       = {'mismatch_darkness_oct21'};
% trial_group_labels      = {'R', 'T', 'RT_gain_up'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'tuning', 'mismatch_darkness_oct21'};


%%
experiment_groups       = {'darkness', 'mismatch_darkness_oct21'};
trial_group_labels      = {{'T_bank', 'T_RT', 'T_R'}, 'T'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'tuning_curves', 'darkness'};



%%
ctl                     = RC2Analysis();
ctl.setup_figures(figure_dir, save_figs);


for ll = 1 : length(trial_group_labels)
    tuning          = {};
    p_svm           = [];
    direction       = [];
    probe_id        = {};
    cluster_id      = [];
    cluster_region  = {};
    
    probe_ids           = ctl.get_probe_ids(experiment_groups{ll});

    for ii = 1 : length(probe_ids)

        data        = ctl.load_formatted_data(probe_ids{ii});
        clusters    = data.VISp_clusters();

        for jj = 1 : length(clusters)


          [~, p_svm(ii, jj), direction(ii, jj)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{ll});
          tuning{ii}{jj} = data.load_tuning_curves(clusters(jj).id, trial_group_labels{ll});

            % extra info to plot on the figure
            probe_id{ii, jj} = probe_ids{ii};
            cluster_id(ii, jj) = clusters(jj).id;
            cluster_region{ii, jj} = clusters(jj).region_str;
        end
    end

        

    %% plot

    for ii = 1 : length(probe_ids)  

        for jj = 1 : length(tuning{ii})

            h_fig                   = ctl.figs.a4figure();
            plot_array              = PlotArray(3, 2);

            tuning_curve_plot       = {};
            shuff_histogram         = {};


            pos         = plot_array.get_position(2-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);

            tuning_curve_plot = TuningCurvePlot(h_ax);

            multicol = lines(2);

            if p_svm(ii, jj) < 0.05 && direction(ii, jj) == 1
                % tonic increase
                main_col = [0.85, 0.32, 0.1];
            elseif p_svm(ii, jj) < 0.05 && direction(ii, jj) == -1
                % tonic decrease
                main_col = [0, 0.45, 0.75];
            else
                main_col = 'k';
            end

            tuning_curve_plot.main_col = main_col;
            tuning_curve_plot.plot(tuning{ii}{jj});
            title(gca, trial_group_labels{ll}, 'interpreter', 'none');

            tuning_curve_plot.xlabel('Speed (cm/s)');
            tuning_curve_plot.ylabel('Firing rate (Hz)');

            pos         = plot_array.get_position(2);
            h_ax2        = axes('units', 'centimeters', 'position', pos);

            shuff_histogram = TuningCurveHistogram(h_ax2);
            shuff_histogram.plot(tuning{ii}{jj});

            mx              = tuning_curve_plot.xmin;
            Mx              = tuning_curve_plot.xmax;
            My              = tuning_curve_plot.ymax;

            tuning_curve_plot.xlim([mx, Mx]);
            tuning_curve_plot.ylim([0, My]);

            FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                    probe_id{ii, jj}, ...
                    cluster_id(ii, jj), ...
                    cluster_region{ii, jj}));

            ctl.figs.save_fig_to_join();
        end

        fname = sprintf('%s.pdf', probe_ids{ii});
        ctl.figs.join_figs(fname, overwrite);
        ctl.figs.clear_figs();
    end

end

