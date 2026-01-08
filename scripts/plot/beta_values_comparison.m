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
    
    % Initialize storage for this probe (3 betas per cluster)
    n_clusters = length(clusters);
    beta_velocity = nan(n_clusters, 3);  % 3 coefficients from 2nd order polynomial
    beta_acceleration = nan(n_clusters, 3);
    
    sig_velocity = false(n_clusters, 1);
    sig_acceleration = false(n_clusters, 1);
    is_linear_velocity = false(n_clusters, 1);
    is_linear_acceleration = false(n_clusters, 1);
    
    % Load tuning data for each cluster
    for c = 1:n_clusters
        cluster_id = clusters(c).id;
        
        try
            % Load velocity tuning
            tuning_vel = data.load_tuning_curves(cluster_id, trial_group_label);
            if isstruct(tuning_vel) && isfield(tuning_vel, 'shuffled') && ...
               isstruct(tuning_vel.shuffled) && isfield(tuning_vel.shuffled, 'beta') && ...
               length(tuning_vel.shuffled.beta) >= 3
                beta_velocity(c, :) = tuning_vel.shuffled.beta;  % All 3 coefficients
                sig_velocity(c) = tuning_vel.shuffled.p < p_thresh;
                
                % F-test to determine if quadratic is better than linear
                % Use raw data from tuning structure
                if sig_velocity(c) && isfield(tuning_vel, 'tuning') && isfield(tuning_vel, 'bin_centers')
                    x_data = tuning_vel.bin_centers(:);
                    y_data = nanmean(tuning_vel.tuning, 2);  % Average firing rate per bin
                    
                    % Remove NaN values
                    valid_idx = ~isnan(x_data) & ~isnan(y_data);
                    x_data = x_data(valid_idx);
                    y_data = y_data(valid_idx);
                    n = length(y_data);
                    
                    if n > 3
                        % Linear fit
                        beta_lin = polyfit(x_data, y_data, 1);
                        yfit_lin = polyval(beta_lin, x_data);
                        ss_res_lin = sum((y_data - yfit_lin).^2);
                        
                        % Quadratic fit (already computed)
                        yfit_quad = polyval(tuning_vel.shuffled.beta, x_data);
                        ss_res_quad = sum((y_data - yfit_quad).^2);
                        
                        % F-test for nested models
                        df1 = 1;  % difference in parameters (3 - 2)
                        df2 = n - 3;  % residual degrees of freedom for quadratic
                        
                        if ss_res_quad > 0 && df2 > 0
                            F = ((ss_res_lin - ss_res_quad) / df1) / (ss_res_quad / df2);
                            p_quadratic = 1 - fcdf(F, df1, df2);
                            
                            % If p > 0.05, quadratic doesn't significantly improve fit
                            is_linear_velocity(c) = (p_quadratic > 0.05);
                        end
                    end
                end
            end
        catch ME
            fprintf('  Warning: Could not load velocity tuning for cluster %d: %s\n', cluster_id, ME.message);
        end
        
        try
            % Load acceleration tuning (returns cell array {all, acc, dec})
            % Use index 1 for "all" acceleration data
            tuning_acc = data.load_tuning_curves_acceleration(cluster_id, trial_group_label);
            
            if iscell(tuning_acc) && length(tuning_acc) >= 1
                % Acceleration "all" (index 1)
                if isstruct(tuning_acc{1}) && isfield(tuning_acc{1}, 'shuffled') && ...
                   isstruct(tuning_acc{1}.shuffled) && isfield(tuning_acc{1}.shuffled, 'beta') && ...
                   length(tuning_acc{1}.shuffled.beta) >= 3
                    beta_acceleration(c, :) = tuning_acc{1}.shuffled.beta;  % All 3 coefficients
                    sig_acceleration(c) = tuning_acc{1}.shuffled.p < p_thresh;
                    
                    % F-test to determine if quadratic is better than linear
                    if sig_acceleration(c) && isfield(tuning_acc{1}, 'tuning') && isfield(tuning_acc{1}, 'bin_centers')
                        x_data = tuning_acc{1}.bin_centers(:);
                        y_data = nanmean(tuning_acc{1}.tuning, 2);  % Average firing rate per bin
                        
                        % Remove NaN values
                        valid_idx = ~isnan(x_data) & ~isnan(y_data);
                        x_data = x_data(valid_idx);
                        y_data = y_data(valid_idx);
                        n = length(y_data);
                        
                        if n > 3
                            % Linear fit
                            beta_lin = polyfit(x_data, y_data, 1);
                            yfit_lin = polyval(beta_lin, x_data);
                            ss_res_lin = sum((y_data - yfit_lin).^2);
                            
                            % Quadratic fit (already computed)
                            yfit_quad = polyval(tuning_acc{1}.shuffled.beta, x_data);
                            ss_res_quad = sum((y_data - yfit_quad).^2);
                            
                            % F-test for nested models
                            df1 = 1;  % difference in parameters (3 - 2)
                            df2 = n - 3;  % residual degrees of freedom for quadratic
                            
                            if ss_res_quad > 0 && df2 > 0
                                F = ((ss_res_lin - ss_res_quad) / df1) / (ss_res_quad / df2);
                                p_quadratic = 1 - fcdf(F, df1, df2);
                                
                                % If p > 0.05, quadratic doesn't significantly improve fit
                                is_linear_acceleration(c) = (p_quadratic > 0.05);
                            end
                        end
                    end
                end
            end
        catch ME
            fprintf('  Warning: Could not load acceleration tuning for cluster %d: %s\n', cluster_id, ME.message);
        end
    end
    
    % X positions for the 3 beta coefficients
    x_positions = [1, 2, 3];  % beta0, beta1, beta2
    
    % Plot velocity tuning (left column)
    subplot(n_rows, n_cols, (pid-1)*n_cols + 1);
    hold on;
    for c = 1:n_clusters
        % Reverse order: polyfit gives [x², x¹, x⁰], we want [x⁰, x¹, x²]
        beta_values = beta_velocity(c, [3, 2, 1]);
        
        % Skip if all values are NaN
        if all(isnan(beta_values))
            continue;
        end
        
        % Set color based on significance and linearity
        if sig_velocity(c)
            if is_linear_velocity(c)
                line_color = [0, 0.6, 0.3];  % Green for significant linear
                marker_face = [0, 0.6, 0.3];
            else
                line_color = [0.5, 0, 0.5];  % Purple for significant quadratic
                marker_face = [0.5, 0, 0.5];
            end
        else
            line_color = [0.7, 0.7, 0.7];  % Gray for not significant
            marker_face = [0.7, 0.7, 0.7];
        end
        
        % Plot line connecting the points
        plot(x_positions, beta_values, '-o', ...
            'Color', line_color, ...
            'MarkerFaceColor', marker_face, ...
            'MarkerEdgeColor', line_color, ...
            'MarkerSize', 4, ...
            'LineWidth', 1);
    end
    hold off;
    xlim([0.5, 3.5]);
    set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'\beta_0', '\beta_1', '\beta_2'});
    set(gca, 'YScale', 'log');
    ylabel('Beta value (log scale)');
    title(sprintf('%s - Velocity (n=%d)', probe_ids{pid}, n_clusters), 'Interpreter', 'none');
    grid on;
    
    % Plot acceleration tuning (right column)
    subplot(n_rows, n_cols, (pid-1)*n_cols + 2);
    hold on;
    for c = 1:n_clusters
        % Reverse order: polyfit gives [x², x¹, x⁰], we want [x⁰, x¹, x²]
        beta_values = beta_acceleration(c, [3, 2, 1]);
        
        % Skip if all values are NaN
        if all(isnan(beta_values))
            continue;
        end
        
        % Set color based on significance and linearity
        if sig_acceleration(c)
            if is_linear_acceleration(c)
                line_color = [0, 0.6, 0.3];  % Green for significant linear
                marker_face = [0, 0.6, 0.3];
            else
                line_color = [0.5, 0, 0.5];  % Purple for significant quadratic
                marker_face = [0.5, 0, 0.5];
            end
        else
            line_color = [0.7, 0.7, 0.7];  % Gray for not significant
            marker_face = [0.7, 0.7, 0.7];
        end
        
        % Plot line connecting the points
        plot(x_positions, beta_values, '-o', ...
            'Color', line_color, ...
            'MarkerFaceColor', marker_face, ...
            'MarkerEdgeColor', line_color, ...
            'MarkerSize', 4, ...
            'LineWidth', 1);
    end
    hold off;
    xlim([0.5, 3.5]);
    set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'\beta_0', '\beta_1', '\beta_2'});
    set(gca, 'YScale', 'log');
    ylabel('Beta value (log scale)');
    title(sprintf('%s - Acceleration (n=%d)', probe_ids{pid}, n_clusters), 'Interpreter', 'none');
    grid on;
    
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
    n_clusters = length(clusters);
    
    n_vel_tuned = 0;
    n_acc_tuned = 0;
    n_vel_linear = 0;
    n_vel_quadratic = 0;
    n_acc_linear = 0;
    n_acc_quadratic = 0;
    
    for c = 1:n_clusters
        cluster_id = clusters(c).id;
        
        try
            tuning_vel = data.load_tuning_curves(cluster_id, trial_group_label);
            if isstruct(tuning_vel) && isfield(tuning_vel, 'shuffled') && ...
               isstruct(tuning_vel.shuffled) && isfield(tuning_vel.shuffled, 'p')
                if tuning_vel.shuffled.p < p_thresh
                    n_vel_tuned = n_vel_tuned + 1;
                    
                    % Classify as linear or quadratic using F-test
                    if isfield(tuning_vel, 'tuning') && isfield(tuning_vel, 'bin_centers')
                        x_data = tuning_vel.bin_centers(:);
                        y_data = nanmean(tuning_vel.tuning, 2);
                        
                        % Remove NaN values
                        valid_idx = ~isnan(x_data) & ~isnan(y_data);
                        x_data = x_data(valid_idx);
                        y_data = y_data(valid_idx);
                        n = length(y_data);
                        
                        if n > 3
                            % Linear fit
                            beta_lin = polyfit(x_data, y_data, 1);
                            yfit_lin = polyval(beta_lin, x_data);
                            ss_res_lin = sum((y_data - yfit_lin).^2);
                            
                            % Quadratic fit
                            yfit_quad = polyval(tuning_vel.shuffled.beta, x_data);
                            ss_res_quad = sum((y_data - yfit_quad).^2);
                            
                            % F-test
                            df1 = 1;
                            df2 = n - 3;
                            
                            if ss_res_quad > 0 && df2 > 0
                                F = ((ss_res_lin - ss_res_quad) / df1) / (ss_res_quad / df2);
                                p_quadratic = 1 - fcdf(F, df1, df2);
                                
                                if p_quadratic > 0.05
                                    n_vel_linear = n_vel_linear + 1;
                                else
                                    n_vel_quadratic = n_vel_quadratic + 1;
                                end
                            end
                        end
                    end
                end
            end
        catch
        end
        
        try
            tuning_acc = data.load_tuning_curves_acceleration(cluster_id, trial_group_label);
            if iscell(tuning_acc) && length(tuning_acc) >= 1
                if isstruct(tuning_acc{1}) && isfield(tuning_acc{1}, 'shuffled') && ...
                   isstruct(tuning_acc{1}.shuffled) && isfield(tuning_acc{1}.shuffled, 'p')
                    if tuning_acc{1}.shuffled.p < p_thresh
                        n_acc_tuned = n_acc_tuned + 1;
                        
                        % Classify as linear or quadratic using F-test
                        if isfield(tuning_acc{1}, 'tuning') && isfield(tuning_acc{1}, 'bin_centers')
                            x_data = tuning_acc{1}.bin_centers(:);
                            y_data = nanmean(tuning_acc{1}.tuning, 2);
                            
                            % Remove NaN values
                            valid_idx = ~isnan(x_data) & ~isnan(y_data);
                            x_data = x_data(valid_idx);
                            y_data = y_data(valid_idx);
                            n = length(y_data);
                            
                            if n > 3
                                % Linear fit
                                beta_lin = polyfit(x_data, y_data, 1);
                                yfit_lin = polyval(beta_lin, x_data);
                                ss_res_lin = sum((y_data - yfit_lin).^2);
                                
                                % Quadratic fit
                                yfit_quad = polyval(tuning_acc{1}.shuffled.beta, x_data);
                                ss_res_quad = sum((y_data - yfit_quad).^2);
                                
                                % F-test
                                df1 = 1;
                                df2 = n - 3;
                                
                                if ss_res_quad > 0 && df2 > 0
                                    F = ((ss_res_lin - ss_res_quad) / df1) / (ss_res_quad / df2);
                                    p_quadratic = 1 - fcdf(F, df1, df2);
                                    
                                    if p_quadratic > 0.05
                                        n_acc_linear = n_acc_linear + 1;
                                    else
                                        n_acc_quadratic = n_acc_quadratic + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        catch
        end
    end
    
    probe_stats(pid).probe_id = probe_ids{pid};
    probe_stats(pid).n_clusters = n_clusters;
    probe_stats(pid).n_vel_tuned = n_vel_tuned;
    probe_stats(pid).n_acc_tuned = n_acc_tuned;
    probe_stats(pid).pct_vel_tuned = (n_vel_tuned / n_clusters) * 100;
    probe_stats(pid).pct_acc_tuned = (n_acc_tuned / n_clusters) * 100;
    probe_stats(pid).n_vel_linear = n_vel_linear;
    probe_stats(pid).n_vel_quadratic = n_vel_quadratic;
    probe_stats(pid).n_acc_linear = n_acc_linear;
    probe_stats(pid).n_acc_quadratic = n_acc_quadratic;
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
