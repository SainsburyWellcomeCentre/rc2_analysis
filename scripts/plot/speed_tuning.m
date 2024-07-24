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


% Initialize parameters
experiment_groups       = {'darkness', 'mismatch_darkness_oct21'};
trial_group_labels      = {{'T_bank', 'T_RT', 'T_R'}, 'T'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'tuning_curves', 'darkness'};


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

            % Store tuning table
            tuning{ii}{jj} = data.load_tuning_curves(clusters(jj).id, trial_group_labels{ll});

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
            h_fig                   = ctl.figs.a4figure();
            plot_array              = PlotArray(3, 2);

            tuning_curve_plot       = {};
            shuff_histogram         = {};

            % Set up first subplot: the tuning curve plot
            pos         = plot_array.get_position(2-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);

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
            tuning_curve_plot.plot(tuning{ii}{jj});
            title(gca, trial_group_labels{ll}, 'interpreter', 'none');

            tuning_curve_plot.xlabel('Speed (cm/s)');
            tuning_curve_plot.ylabel('Firing rate (Hz)');

            % Set up second subplot: the tuning curve histogram
            pos         = plot_array.get_position(2);
            h_ax2        = axes('units', 'centimeters', 'position', pos);

            % Instantiate TuningCurveHistogram object
            shuff_histogram = TuningCurveHistogram(h_ax2);

            % Evaluate TuningCurveHistogram attributes (tuning, max min limits)
            shuff_histogram.plot(tuning{ii}{jj});

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

            % Save plot as a separated pdf file
            ctl.figs.save_fig_to_join();
        end

        % Join pdfs from a unique probe
        fname = sprintf('%s.pdf', probe_ids{ii});
        ctl.figs.join_figs(fname, overwrite);
        ctl.figs.clear_figs();
    end
end

