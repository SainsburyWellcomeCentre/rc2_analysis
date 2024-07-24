% Plot acceleration tuning curve for an experimental group
% Tuning cuves must have been generated and saved using
% RC2Analysis.create_tuning_curves_acceleration

% Acceleration data consist of negative (deceleration) and positive ("real" acceleration) values, 
% and they can be binned in different ways (modalities).
% 1. All: maintaining both acceletation and deceleration in a matrix and then binning
% 2. Acc: retaining only positive acceleration values, and then binning
% 3. Dec: retaining only negative acceleration values, and then binning

% Initialize parameters
experiment_groups       = {'passive_same_luminance'};
trial_group_labels      = {'T_Vstatic'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'tuning_curves', 'passive_same_luminance_acc'};

modalities = ["all", "acc", "dec"];

% Creating controller object
ctl                     = RC2Analysis();
ctl.setup_figures(figure_dir, save_figs);

% Iterate across coditions
for ll = 1 : length(trial_group_labels)
    % Collect all the data befor plotting
    % Initialize empty matrices and cell arrays
    tuning          = {};
    p_svm           = [];
    direction       = [];
    probe_id        = {};
    cluster_id      = [];
    cluster_region  = {};

    % Get probe ids for a given experimental group
    probe_ids           = ctl.get_probe_ids(experiment_groups{ll});

    % Iterate thorugh probe ids
    for ii = 1 : length(probe_ids)
        % Get the data for VISp clusters
        data        = ctl.load_formatted_data(probe_ids{ii});
        clusters    = data.VISp_clusters();

        % Loop across clusters
        for jj = 1 : length(clusters)
            % Get statistics on stationary vs motion (read from svm table) 
            [~, p_svm(ii, jj), direction(ii, jj)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{ll});
            
            % Store tuning table (contains the three modalities)
            tuning{ii}{jj} = data.load_tuning_curves_acceleration(clusters(jj).id, trial_group_labels{ll});

            % Collect further info for plotting
            probe_id{ii, jj} = probe_ids{ii};
            cluster_id(ii, jj) = clusters(jj).id;
            cluster_region{ii, jj} = clusters(jj).region_str;
        end
    end

        
    % Make the plots
    % Loop across clusters
    for ii = 1 : length(probe_ids)  
        % Loop across speed tuning bins
        for jj = 1 : length(tuning{ii})
            % Setup figure
            h_fig  = ctl.figs.a4figure();

            % loop across modalities
            for kk = 1 : length(tuning{ii}{jj})
                % create a PlotArray  (object to derive positions in A4 space)
                plot_array = PlotArray(3, 2);

                tuning_curve_plot = {};
                shuff_histogram   = {};

                % Set up first subplot: the tuning curve plot
                % Set position according to which modality we are plottinf
                pos    = plot_array.get_position(kk * 2 - 1);
                pos(2) = pos(2) - 5;
                h_ax   = axes('units', 'centimeters', 'position', pos);
                
                % Instantiate TuningCurvePlot object
                tuning_curve_plot = TuningCurvePlot(h_ax);

                if p_svm(ii, jj) < 0.05 && direction(ii, jj) == 1
                    % Color the plot orange if the cluster was positively modulated (tonic increase)
                    main_col = [0.85, 0.32, 0.1];
                elseif p_svm(ii, jj) < 0.05 && direction(ii, jj) == -1
                    % Color the plot cyan if the cluster was negatively modulated (tonic decrease)
                    main_col = [0, 0.45, 0.75];
                else
                    % Color the plot black if the cluster was not modulated 
                    main_col = 'k';
                end

                % Evaluate TuningCurvePlot attributes (color, tuning, title, labels)
                tuning_curve_plot.main_col = main_col;
                tuning_curve_plot.plot(tuning{ii}{jj}{kk});
                title(gca, trial_group_labels{ll}, 'interpreter', 'none');

                tuning_curve_plot.xlabel('Acceleration (cm/s^2)');
                tuning_curve_plot.ylabel('Firing rate (Hz)');

                % Set up second subplot: the tuning curve histogram
                pos    = plot_array.get_position(kk * 2);
                pos(2) = pos(2) - 5;
                h_ax2  = axes('units', 'centimeters', 'position', pos);

                % Instantiate TuningCurveHistogram object
                shuff_histogram = TuningCurveHistogram(h_ax2);

                % Evaluate TuningCurveHistogram attributes (tuning, max min limits)
                shuff_histogram.plot(tuning{ii}{jj}{kk});

                mx              = tuning_curve_plot.xmin;
                Mx              = tuning_curve_plot.xmax;
                My              = tuning_curve_plot.ymax;

                tuning_curve_plot.xlim([mx, Mx]);
                tuning_curve_plot.ylim([0, My]);

                % Set figure title
                FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                        probe_id{ii, jj}, ...
                        cluster_id(ii, jj), ...
                        cluster_region{ii, jj}));
            end
            % Save plot as a separated pdf file
            ctl.figs.save_fig_to_join();
        end

        % Join pdfs from a unique probe
        fname = sprintf('%s.pdf', probe_ids{ii});
        ctl.figs.join_figs(fname, overwrite);
        ctl.figs.clear_figs();
    end

end

