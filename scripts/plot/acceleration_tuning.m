%TODO

close all
%%
experiment_groups       = {'darkness', 'mismatch_darkness_oct21'};
trial_group_labels      = {{'T_bank', 'T_RT', 'T_R'}, 'T'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'tuning_curves', 'darkness_acc'};

modalities = ["all", "acc", "dec"];

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
          tuning{ii}{jj} = data.load_tuning_curves_acceleration(clusters(jj).id, trial_group_labels{ll});

            % extra info to plot on the figure
            probe_id{ii, jj} = probe_ids{ii};
            cluster_id(ii, jj) = clusters(jj).id;
            cluster_region{ii, jj} = clusters(jj).region_str;
        end
    end

        

    %% plot

    for ii = 1 : length(probe_ids)  

        for jj = 1 : length(tuning{ii})

            h_fig  = ctl.figs.a4figure();

            for kk = 1 : length(tuning{ii}{jj})

                plot_array              = PlotArray(3, 2);

                tuning_curve_plot       = {};
                shuff_histogram         = {};


                pos         = plot_array.get_position(kk * 2 - 1);
                pos(2) = pos(2) - 5;
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
                tuning_curve_plot.plot(tuning{ii}{jj}{kk});
                title(gca, trial_group_labels{ll}, 'interpreter', 'none');

                tuning_curve_plot.xlabel('Acceleration (cm/s^2)');
                tuning_curve_plot.ylabel('Firing rate (Hz)');

                pos         = plot_array.get_position(kk * 2);
                pos(2) = pos(2) - 5;
                h_ax2        = axes('units', 'centimeters', 'position', pos);

                shuff_histogram = TuningCurveHistogram(h_ax2);
                shuff_histogram.plot(tuning{ii}{jj}{kk});

                mx              = tuning_curve_plot.xmin;
                Mx              = tuning_curve_plot.xmax;
                My              = tuning_curve_plot.ymax;

                tuning_curve_plot.xlim([mx, Mx]);
                tuning_curve_plot.ylim([0, My]);



                FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                        probe_id{ii, jj}, ...
                        cluster_id(ii, jj), ...
                        cluster_region{ii, jj}));
            end
            ctl.figs.save_fig_to_join();
        end

        fname = sprintf('%s.pdf', probe_ids{ii});
        ctl.figs.join_figs(fname, overwrite);
        ctl.figs.clear_figs();
    end

end

