% Plot beta values from velocity and acceleration tunings across probes
% Creates a figure with rows for each probe and columns for velocity/acceleration:
%   - Rows: one per probe
%   - Columns: Velocity (left) and Acceleration-all (right)
%   - X-axis: categories of beta (beta0, beta1, beta2, beta3) from 3rd order polynomial fit
%   - Y-axis: values of beta coefficients
%   - Dots connected by lines to track beta values for every cluster
%   - Black lines if significant, gray if not significant
%
% Requires:
%   - Pre-computed tuning curves from create_tables.m (velocity)
%   - Pre-computed tuning curves from create_tables_acceleration.m (acceleration)

%%
% Configuration
experiment_groups    = {'ambient_light'};
trial_group_label    = 'RT';  % Must match the label used in create_tables.m
save_figs            = true;
overwrite            = true;
figure_dir           = {'beta_comparison', 'ambient_light'};

% Significance threshold
p_thresh = 0.05;

% Initialize controller
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));

% Subplot layout: each row is a probe, 2 columns (velocity, acceleration)
n_probes = length(probe_ids);
n_rows = n_probes;
n_cols = 2;  % Velocity and Acceleration

% Create figure
fig = ctl.figs.a4figure('portrait');

for pid = 1:length(probe_ids)
    fprintf('\nProcessing probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    
    % Load data
    data = ctl.load_formatted_data(probe_ids{pid});
    clusters = data.selected_clusters();
    
    fprintf('  Found %d clusters\n', length(clusters));
    
    % Load and classify velocity tuning
    [beta_velocity, sig_velocity, is_linear_velocity] = load_and_classify_tuning(...
        data, clusters, trial_group_label, 'velocity', p_thresh);
    
    % Load and classify acceleration tuning
    [beta_acceleration, sig_acceleration, is_linear_acceleration] = load_and_classify_tuning(...
        data, clusters, trial_group_label, 'acceleration', p_thresh);
    
    % Plot velocity tuning (left column)
    subplot(n_rows, n_cols, (pid-1)*n_cols + 1);
    plot_beta_coefficients_subplot(beta_velocity, sig_velocity, is_linear_velocity, ...
        sprintf('%s - Velocity (n=%d)', probe_ids{pid}, length(clusters)));
    
    % Plot acceleration tuning (right column)
    subplot(n_rows, n_cols, (pid-1)*n_cols + 2);
    plot_beta_coefficients_subplot(beta_acceleration, sig_acceleration, is_linear_acceleration, ...
        sprintf('%s - Acceleration (n=%d)', probe_ids{pid}, length(clusters)));
    
    fprintf('  Completed probe %d/%d\n', pid, length(probe_ids));
end

% Add legend to the figure
% Create invisible dummy plots for legend
fig_handle = gcf;
ax_legend = axes(fig_handle, 'Visible', 'off');
hold(ax_legend, 'on');
h_linear = plot(ax_legend, nan, nan, '-o', 'Color', [0, 0.6, 0.3], ...
    'MarkerFaceColor', [0, 0.6, 0.3], 'MarkerSize', 6, 'LineWidth', 1.5);
h_quadratic = plot(ax_legend, nan, nan, '-o', 'Color', [0.5, 0, 0.5], ...
    'MarkerFaceColor', [0.5, 0, 0.5], 'MarkerSize', 6, 'LineWidth', 1.5);
h_nonsig = plot(ax_legend, nan, nan, '-o', 'Color', [0.7, 0.7, 0.7], ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerSize', 6, 'LineWidth', 1.5);
legend(ax_legend, [h_linear, h_quadratic, h_nonsig], ...
    {'Linear (significant)', 'Quadratic (significant)', 'Not significant'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
hold(ax_legend, 'off');

% Collect statistics for summary
probe_stats = struct();
for pid = 1:length(probe_ids)
    data = ctl.load_formatted_data(probe_ids{pid});
    clusters = data.selected_clusters();
    
    % Compute statistics using helper function
    stats = compute_tuning_statistics(data, clusters, trial_group_label, p_thresh);
    
    probe_stats(pid).probe_id = probe_ids{pid};
    probe_stats(pid).n_clusters = stats.n_clusters;
    probe_stats(pid).n_vel_tuned = stats.n_vel_tuned;
    probe_stats(pid).n_acc_tuned = stats.n_acc_tuned;
    probe_stats(pid).pct_vel_tuned = (stats.n_vel_tuned / stats.n_clusters) * 100;
    probe_stats(pid).pct_acc_tuned = (stats.n_acc_tuned / stats.n_clusters) * 100;
    probe_stats(pid).n_vel_linear = stats.n_vel_linear;
    probe_stats(pid).n_vel_quadratic = stats.n_vel_quadratic;
    probe_stats(pid).n_acc_linear = stats.n_acc_linear;
    probe_stats(pid).n_acc_quadratic = stats.n_acc_quadratic;
end

% Write summary to text file
summary_file = fullfile(ctl.figs.curr_dir, 'tuning_summary.txt');
fid = fopen(summary_file, 'w');

fprintf(fid, '=== Velocity and Acceleration Tuning Summary ===\n');
fprintf(fid, 'Experiment group: %s\n', strjoin(experiment_groups, ', '));
fprintf(fid, 'Trial group label: %s\n', trial_group_label);
fprintf(fid, 'Significance threshold: p < %.2f\n\n', p_thresh);

fprintf(fid, '--- Per Probe Statistics ---\n');
fprintf(fid, '%-20s %10s %15s %15s %12s %12s %15s %15s %12s %12s\n', ...
    'Probe ID', 'N Clusters', 'Vel Tuned', 'Vel Tuned %%', 'Vel Linear%%', 'Vel Quad%%', ...
    'Acc Tuned', 'Acc Tuned %%', 'Acc Linear%%', 'Acc Quad%%');
fprintf(fid, '%s\n', repmat('-', 1, 140));

% Separate probes by region (A=RSP, B=CA1)
rsp_stats = [];
ca1_stats = [];

for pid = 1:length(probe_stats)
    % Calculate percentages for linear/quadratic
    pct_vel_linear = 0;
    pct_vel_quad = 0;
    pct_acc_linear = 0;
    pct_acc_quad = 0;
    
    if probe_stats(pid).n_vel_tuned > 0
        pct_vel_linear = (probe_stats(pid).n_vel_linear / probe_stats(pid).n_vel_tuned) * 100;
        pct_vel_quad = (probe_stats(pid).n_vel_quadratic / probe_stats(pid).n_vel_tuned) * 100;
    end
    if probe_stats(pid).n_acc_tuned > 0
        pct_acc_linear = (probe_stats(pid).n_acc_linear / probe_stats(pid).n_acc_tuned) * 100;
        pct_acc_quad = (probe_stats(pid).n_acc_quadratic / probe_stats(pid).n_acc_tuned) * 100;
    end
    
    fprintf(fid, '%-20s %10d %15d %14.1f%% %11.1f%% %11.1f%% %15d %14.1f%% %11.1f%% %11.1f%%\n', ...
        probe_stats(pid).probe_id, probe_stats(pid).n_clusters, ...
        probe_stats(pid).n_vel_tuned, probe_stats(pid).pct_vel_tuned, pct_vel_linear, pct_vel_quad, ...
        probe_stats(pid).n_acc_tuned, probe_stats(pid).pct_acc_tuned, pct_acc_linear, pct_acc_quad);
    
    % Determine region based on probe ID naming convention (a_rec = RSP, b_rec = CA1)
    probe_name = probe_stats(pid).probe_id;
    if contains(probe_name, 'a_rec', 'IgnoreCase', true)
        rsp_stats = [rsp_stats, pid];
    elseif contains(probe_name, 'b_rec', 'IgnoreCase', true)
        ca1_stats = [ca1_stats, pid];
    end
end

fprintf(fid, '\n--- Grouped by Region ---\n');

% RSP (A) statistics
if ~isempty(rsp_stats)
    total_clusters_rsp = sum([probe_stats(rsp_stats).n_clusters]);
    total_vel_rsp = sum([probe_stats(rsp_stats).n_vel_tuned]);
    total_acc_rsp = sum([probe_stats(rsp_stats).n_acc_tuned]);
    pct_vel_rsp = (total_vel_rsp / total_clusters_rsp) * 100;
    pct_acc_rsp = (total_acc_rsp / total_clusters_rsp) * 100;
    
    total_vel_linear_rsp = sum([probe_stats(rsp_stats).n_vel_linear]);
    total_vel_quadratic_rsp = sum([probe_stats(rsp_stats).n_vel_quadratic]);
    total_acc_linear_rsp = sum([probe_stats(rsp_stats).n_acc_linear]);
    total_acc_quadratic_rsp = sum([probe_stats(rsp_stats).n_acc_quadratic]);
    
    fprintf(fid, '\nRSP (A) - %d probes:\n', length(rsp_stats));
    fprintf(fid, '  Total clusters: %d\n', total_clusters_rsp);
    fprintf(fid, '  Velocity tuned: %d (%.1f%%)\n', total_vel_rsp, pct_vel_rsp);
    if total_vel_rsp > 0
        fprintf(fid, '    Linear: %d (%.1f%% of tuned)\n', total_vel_linear_rsp, (total_vel_linear_rsp/total_vel_rsp)*100);
        fprintf(fid, '    Quadratic: %d (%.1f%% of tuned)\n', total_vel_quadratic_rsp, (total_vel_quadratic_rsp/total_vel_rsp)*100);
    end
    fprintf(fid, '  Acceleration tuned: %d (%.1f%%)\n', total_acc_rsp, pct_acc_rsp);
    if total_acc_rsp > 0
        fprintf(fid, '    Linear: %d (%.1f%% of tuned)\n', total_acc_linear_rsp, (total_acc_linear_rsp/total_acc_rsp)*100);
        fprintf(fid, '    Quadratic: %d (%.1f%% of tuned)\n', total_acc_quadratic_rsp, (total_acc_quadratic_rsp/total_acc_rsp)*100);
    end
end

% CA1 (B) statistics
if ~isempty(ca1_stats)
    total_clusters_ca1 = sum([probe_stats(ca1_stats).n_clusters]);
    total_vel_ca1 = sum([probe_stats(ca1_stats).n_vel_tuned]);
    total_acc_ca1 = sum([probe_stats(ca1_stats).n_acc_tuned]);
    pct_vel_ca1 = (total_vel_ca1 / total_clusters_ca1) * 100;
    pct_acc_ca1 = (total_acc_ca1 / total_clusters_ca1) * 100;
    
    total_vel_linear_ca1 = sum([probe_stats(ca1_stats).n_vel_linear]);
    total_vel_quadratic_ca1 = sum([probe_stats(ca1_stats).n_vel_quadratic]);
    total_acc_linear_ca1 = sum([probe_stats(ca1_stats).n_acc_linear]);
    total_acc_quadratic_ca1 = sum([probe_stats(ca1_stats).n_acc_quadratic]);
    
    fprintf(fid, '\nCA1 (B) - %d probes:\n', length(ca1_stats));
    fprintf(fid, '  Total clusters: %d\n', total_clusters_ca1);
    fprintf(fid, '  Velocity tuned: %d (%.1f%%)\n', total_vel_ca1, pct_vel_ca1);
    if total_vel_ca1 > 0
        fprintf(fid, '    Linear: %d (%.1f%% of tuned)\n', total_vel_linear_ca1, (total_vel_linear_ca1/total_vel_ca1)*100);
        fprintf(fid, '    Quadratic: %d (%.1f%% of tuned)\n', total_vel_quadratic_ca1, (total_vel_quadratic_ca1/total_vel_ca1)*100);
    end
    fprintf(fid, '  Acceleration tuned: %d (%.1f%%)\n', total_acc_ca1, pct_acc_ca1);
    if total_acc_ca1 > 0
        fprintf(fid, '    Linear: %d (%.1f%% of tuned)\n', total_acc_linear_ca1, (total_acc_linear_ca1/total_acc_ca1)*100);
        fprintf(fid, '    Quadratic: %d (%.1f%% of tuned)\n', total_acc_quadratic_ca1, (total_acc_quadratic_ca1/total_acc_ca1)*100);
    end
end

% Overall statistics
fprintf(fid, '\n--- Overall Statistics ---\n');
total_clusters_all = sum([probe_stats.n_clusters]);
total_vel_all = sum([probe_stats.n_vel_tuned]);
total_acc_all = sum([probe_stats.n_acc_tuned]);
pct_vel_all = (total_vel_all / total_clusters_all) * 100;
pct_acc_all = (total_acc_all / total_clusters_all) * 100;

total_vel_linear_all = sum([probe_stats.n_vel_linear]);
total_vel_quadratic_all = sum([probe_stats.n_vel_quadratic]);
total_acc_linear_all = sum([probe_stats.n_acc_linear]);
total_acc_quadratic_all = sum([probe_stats.n_acc_quadratic]);

fprintf(fid, 'Total probes: %d\n', length(probe_stats));
fprintf(fid, 'Total clusters: %d\n', total_clusters_all);
fprintf(fid, 'Velocity tuned: %d (%.1f%%)\n', total_vel_all, pct_vel_all);
if total_vel_all > 0
    fprintf(fid, '  Linear: %d (%.1f%% of tuned)\n', total_vel_linear_all, (total_vel_linear_all/total_vel_all)*100);
    fprintf(fid, '  Quadratic: %d (%.1f%% of tuned)\n', total_vel_quadratic_all, (total_vel_quadratic_all/total_vel_all)*100);
end
fprintf(fid, 'Acceleration tuned: %d (%.1f%%)\n', total_acc_all, pct_acc_all);
if total_acc_all > 0
    fprintf(fid, '  Linear: %d (%.1f%% of tuned)\n', total_acc_linear_all, (total_acc_linear_all/total_acc_all)*100);
    fprintf(fid, '  Quadratic: %d (%.1f%% of tuned)\n', total_acc_quadratic_all, (total_acc_quadratic_all/total_acc_all)*100);
end

fclose(fid);
fprintf('\nTuning summary saved to: %s\n', summary_file);

% Add figure title
FigureTitle(fig, sprintf('Beta Values Comparison - %s', strjoin(experiment_groups, ', ')));

% Save figure
if save_figs
    fprintf('\nSaving figure...\n');
    ctl.figs.save_fig_to_join();
    fname = sprintf('beta_values_comparison_%s.pdf', strjoin(experiment_groups, '_'));
    ctl.figs.join_figs(fname, overwrite);
    ctl.figs.clear_figs();
end

fprintf('\n=== Analysis complete ===\n');
