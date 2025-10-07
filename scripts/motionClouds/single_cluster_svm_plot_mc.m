% Single-cluster unity plots for Motion Clouds (passive_same_luminance_mc)
% Exact behavior as scripts/plot/unity_plots/single_cluster_svm_unity_plot.m
% but with experiment_groups and trial_group_labels set for Motion Clouds.
%
% Plot stationary vs. motion unity plot for an experiment group and trial
% group for each cluster individually. i.e. show the individual trials on
% the unity plot (c.f. population_svm_unity_and_mi_plot).

close all;

%% Parameters (Motion Clouds)
experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};
marker_style            = {'o', 'o', 'o'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'motionClouds', 'passive_same_luminance_mc', 'stationary_vs_motion', 'single_cluster'};

%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

x_all                   = cell(1, length(probe_ids));
y_all                   = cell(1, length(probe_ids));
p_val                   = cell(1, length(probe_ids));
direction               = cell(1, length(probe_ids));
cluster_ids             = cell(1, length(probe_ids));
region_str              = cell(1, length(probe_ids));

for ii = 1 : length(probe_ids)
    data                = ctl.load_formatted_data(probe_ids{ii});
    clusters            = data.VISp_clusters;

    % get list of cluster IDs and the regions
    cluster_ids{ii}     = data.VISp_cluster_ids;
    region_str{ii}      = {clusters(:).region_str};

    for jj = 1 : length(trial_group_labels)
        % skip if the trial group label is not in the experiment
        if ~data.check_trial_group(trial_group_labels{jj})
            continue
        end

        for kk = 1 : length(cluster_ids{ii})
            x_all{ii}{kk}{jj} = data.stationary_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            y_all{ii}{kk}{jj} = data.motion_fr_for_trial_group(cluster_ids{ii}(kk), trial_group_labels{jj});
            [~, p_val{ii}{kk}(jj), direction{ii}{kk}(jj)] = ...
                data.is_stationary_vs_motion_significant(cluster_ids{ii}(kk), trial_group_labels{jj});
        end
    end
end

%%
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    for jj = 1 : length(cluster_ids{ii})
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        axs                     = gobjects(0);
        panel_mins              = [];
        panel_maxs              = [];

        for kk = 1 : length(trial_group_labels)
            pos         = plot_array.get_position(kk);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            axs(end+1)  = h_ax;
            hold(h_ax, 'on');

            xb = x_all{ii}{jj}{kk};
            yr = y_all{ii}{jj}{kk};
            if isempty(xb), xb = []; end
            if isempty(yr), yr = []; end

            % Paired faint lines
            n_pairs = min(numel(xb), numel(yr));
            for t = 1:n_pairs
                if isfinite(xb(t)) && isfinite(yr(t))
                    line(h_ax, [1 2], [xb(t) yr(t)], 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
                end
            end

            % Baseline points
            if ~isempty(xb)
                plot(h_ax, ones(size(xb))*1, xb, 'ko', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 4);
            end

            % Response points with significance coloring
            if ~isempty(yr)
                resp_color = [0.7 0.7 0.7];
                if numel(p_val{ii})>=jj && numel(p_val{ii}{jj})>=kk && isfinite(p_val{ii}{jj}(kk))
                    pv = p_val{ii}{jj}(kk);
                    dirn = direction{ii}{jj}(kk);
                    if pv < 0.05 && dirn > 0
                        resp_color = [1.0 0.0 0.0];
                    elseif pv < 0.05 && dirn < 0
                        resp_color = [0.0 0.2 1.0];
                    end
                end
                plot(h_ax, ones(size(yr))*2, yr, 'o', 'Color', resp_color, 'MarkerFaceColor', resp_color, 'MarkerSize', 4);
            end

            xlim(h_ax, [0.5 2.5]);
            set(h_ax, 'XTick', [1 2], 'XTickLabel', {'Baseline','Response'});
            ylabel(h_ax, 'Firing rate (Hz)');
            title(h_ax, trial_group_labels{kk});
            grid(h_ax, 'on');

            vals = [xb(:); yr(:)];
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                panel_mins(end+1) = min(vals);
                panel_maxs(end+1) = max(vals);
            else
                panel_mins(end+1) = 0;
                panel_maxs(end+1) = 1;
            end
        end

        % Sync y-limits across panels
        if ~isempty(panel_mins)
            m = min(panel_mins);
            M = max(panel_maxs);
            if m == M, M = m + 1; end
            for kk = 1 : numel(axs)
                ylim(axs(kk), [m, M]);
            end
        end

        % page title
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', probe_ids{ii}, cluster_ids{ii}(jj), region_str{ii}{jj}));

        ctl.figs.save_fig_to_join();
    end

    ctl.figs.join_figs(sprintf('%s.pdf', probe_ids{ii}), overwrite);
    ctl.figs.clear_figs();
end


