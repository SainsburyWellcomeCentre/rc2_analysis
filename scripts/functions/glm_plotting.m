classdef glm_plotting
%GLM_PLOTTING Plotting functions for GLM single cluster analysis
%   Static methods for generating population-level summary figures,
%   per-cluster visualizations, and tuning curve comparison plots.
%
%   Usage:
%       fig = glm_plotting.plot_forward_selection_summary(results, n_clusters, threshold);
%       fig = glm_plotting.plot_correlation_summary(tc_corr_results, model_categories);
%
%   Categories:
%       Population Summaries:   plot_forward_selection_summary
%       Per-Cluster Plots:      plot_cluster_model_overview, plot_cluster_tuning_curves
%       Comparison Plots:       plot_correlation_summary, plot_selected_vs_additive
%       Helpers:                get_condition_colors, get_main_effect_colors,
%                               draw_boxplot, get_model_category_markers

    methods (Static)
        
        %% ================================================================
        %  Color and Style Helpers
        %  ================================================================
        
        function colors = get_condition_colors()
        %GET_CONDITION_COLORS Return standard colors for experimental conditions
        %   colors = glm_plotting.get_condition_colors()
        %
        %   Returns a containers.Map with keys: 'T_Vstatic', 'V', 'VT'
            colors = containers.Map();
            colors('T_Vstatic') = [0.2 0.7 0.2];  % Green
            colors('V') = [0.9 0.6 0.1];          % Orange-yellow
            colors('VT') = [0 0.4 0.8];           % Blue
        end
        
        
        function colors = get_main_effect_colors()
        %GET_MAIN_EFFECT_COLORS Return standard colors for main effects
        %   colors = glm_plotting.get_main_effect_colors()
        %
        %   Returns a struct with fields: speed, tf, sf, or_color, null, interaction
            colors = struct();
            colors.speed = [0.17 0.63 0.17];       % Green
            colors.tf = [1.0 0.50 0.05];           % Orange
            colors.sf = [0.95 0.85 0.10];          % Yellow
            colors.or_color = [0.84 0.15 0.16];    % Red
            colors.null = [0.75 0.75 0.75];        % Light gray
            colors.fourplus = [0.35 0.35 0.35];    % Dark gray
            colors.interaction = [0.20 0.40 0.80]; % Blue
        end
        
        
        function markers = get_model_category_markers()
        %GET_MODEL_CATEGORY_MARKERS Return marker styles for model categories
        %   markers = glm_plotting.get_model_category_markers()
        %
        %   Returns a containers.Map with model category names as keys
            markers = containers.Map();
            markers('None') = 'o';           % Circle - no variables selected
            markers('Speed') = '^';          % Triangle up - Speed only
            markers('TF') = 'v';             % Triangle down - TF only
            markers('Speed+TF') = 's';       % Square - Speed+TF without interaction
            markers('Speed+TF+Int') = 'd';   % Diamond - Speed+TF with interaction
            markers('Visual(SF/OR)') = 'p';  % Pentagon - SF/OR only
            markers('Speed+Visual') = 'h';   % Hexagram - Speed + SF/OR
            markers('TF+Visual') = '+';      % Plus - TF + SF/OR
            markers('Other') = 'x';          % Cross - other combinations
        end
        
        
        function sf_colors = get_sf_colormap()
        %GET_SF_COLORMAP Return colors for SF levels
        %   sf_colors = glm_plotting.get_sf_colormap()
        %
        %   Returns 3x3 matrix for SF levels [0.003, 0.006, 0.012]
            sf_colors = [0.2 0.6 0.2;    % 0.003 cpd = green
                         0.2 0.4 0.9;    % 0.006 cpd = blue
                         0.9 0.2 0.2];   % 0.012 cpd = red
        end
        
        
        function or_colors = get_or_colormap()
        %GET_OR_COLORMAP Return colors for orientation levels
        %   or_colors = glm_plotting.get_or_colormap()
        %
        %   Returns 4x3 matrix for OR levels [-pi/4, 0, pi/4, pi/2]
            or_colors = [0.9 0.1 0.1;    % -45° = bright red
                         0.1 0.8 0.9;    % 0° = cyan
                         0.95 0.55 0.0;  % 45° = orange
                         0.6 0.1 0.85];  % 90° = purple
        end
        
        
        function c_dark = darken_color(c, factor)
        %DARKEN_COLOR Darken a color by a factor
        %   c_dark = glm_plotting.darken_color(c, factor)
        %
        %   Inputs:
        %       c      - RGB color [1x3]
        %       factor - Darkening factor (default: 0.65)
            if nargin < 2, factor = 0.65; end
            c_dark = max(c * factor, 0);
        end
        
        
        %% ================================================================
        %  Box Plot Helper
        %  ================================================================
        
        function draw_boxplot(ax, x_pos, values, clr, box_width)
        %DRAW_BOXPLOT Draw a simple boxplot at the specified position
        %   glm_plotting.draw_boxplot(ax, x_pos, values, clr, box_width)
        %
        %   Inputs:
        %       ax        - Axes handle
        %       x_pos     - X position for the box
        %       values    - Data values for the boxplot
        %       clr       - RGB color [1x3]
        %       box_width - Width of the box (default: 0.6)
            if nargin < 5, box_width = 0.6; end
            
            valid = values(~isnan(values));
            if length(valid) < 2, return; end
            
            axes(ax); %#ok<LAXES>
            hold on;
            
            q = quantile(valid, [0.25, 0.5, 0.75]);
            iqr_val = q(3) - q(1);
            whisker_lo = max(min(valid), q(1) - 1.5*iqr_val);
            whisker_hi = min(max(valid), q(3) + 1.5*iqr_val);
            
            % Box
            rectangle('Position', [x_pos-box_width/2, q(1), box_width, max(q(3)-q(1), 0.01)], ...
                'EdgeColor', clr, 'LineWidth', 1.5, 'FaceColor', [clr, 0.3]);
            % Median line
            plot([x_pos-box_width/2, x_pos+box_width/2], [q(2), q(2)], ...
                'Color', clr, 'LineWidth', 2);
            % Whiskers
            plot([x_pos, x_pos], [whisker_lo, q(1)], 'Color', clr, 'LineWidth', 1);
            plot([x_pos, x_pos], [q(3), whisker_hi], 'Color', clr, 'LineWidth', 1);
            plot([x_pos-0.2, x_pos+0.2], [whisker_lo, whisker_lo], 'Color', clr, 'LineWidth', 1);
            plot([x_pos-0.2, x_pos+0.2], [whisker_hi, whisker_hi], 'Color', clr, 'LineWidth', 1);
        end
        
        
        %% ================================================================
        %  Population Summary Plots
        %  ================================================================
        
        function fig = plot_forward_selection_summary(results, n_clusters, threshold)
        %PLOT_FORWARD_SELECTION_SUMMARY Generate population summary figure
        %   fig = glm_plotting.plot_forward_selection_summary(results, n_clusters, threshold)
        %
        %   Creates a 3-panel figure showing:
        %     Panel 1: Model complexity distribution (bar chart)
        %     Panel 2: Variable inclusion rates (horizontal bar)
        %     Panel 3: Interaction breakdown (horizontal bar)
        %
        %   Inputs:
        %       results    - Results struct with time.* fields (from GLM analysis)
        %       n_clusters - Total number of clusters analyzed
        %       threshold  - Selection threshold (delta bps)
        %
        %   Output:
        %       fig - Figure handle
        
            colors = glm_plotting.get_main_effect_colors();
            
            fig = figure('Position', [50 50 1800 500], 'Name', 'Forward Selection Summary');
            
            gt_tag = 'time';
            n_total = n_clusters;
            
            % Get forward selection results
            n_vars = results.(gt_tag).n_selected_vars;
            is_spd = results.(gt_tag).is_speed_tuned;
            is_tf = results.(gt_tag).is_tf_tuned;
            is_sf = results.(gt_tag).is_sf_tuned;
            is_or = results.(gt_tag).is_or_tuned;
            has_int = results.(gt_tag).has_interaction;
            
            % Get interaction breakdown
            has_spd_tf = results.(gt_tag).has_speed_x_tf;
            has_spd_sf = results.(gt_tag).has_speed_x_sf;
            has_spd_or = results.(gt_tag).has_speed_x_or;
            has_tf_sf = results.(gt_tag).has_tf_x_sf;
            has_tf_or = results.(gt_tag).has_tf_x_or;
            has_sf_or = results.(gt_tag).has_sf_x_or;
            
            % --- Panel 1: Model Complexity Distribution ---
            subplot(1, 3, 1);
            
            selected_vars_list = results.(gt_tag).selected_vars;
            
            % Build model type counts by complexity
            complexity_counts = zeros(5, 1);  % 0, 1, 2, 3, 4+ vars
            for ni = 1:length(selected_vars_list)
                sv = selected_vars_list{ni};
                if isempty(sv), sv = {}; end
                n_total_vars = length(sv);
                bin_idx = min(n_total_vars + 1, 5);
                complexity_counts(bin_idx) = complexity_counts(bin_idx) + 1;
            end
            
            complexity_labels = {'Null', '1 var', '2 vars', '3 vars', '4+ vars'};
            complexity_colors = [colors.null; colors.speed; colors.tf; colors.sf; colors.fourplus];
            
            b = bar(complexity_counts, 'FaceColor', 'flat', 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5);
            b.CData = complexity_colors;
            
            set(gca, 'XTick', 1:5, 'XTickLabel', complexity_labels, 'XTickLabelRotation', 30, 'FontSize', 8);
            ylabel('# Neurons');
            title('Model Complexity (Forward Selection)', 'FontSize', 10);
            set(gca, 'box', 'off');
            
            % Add count and percentage labels
            for ki = 1:length(complexity_counts)
                if complexity_counts(ki) > 0
                    text(ki, complexity_counts(ki) + max(complexity_counts)*0.03, ...
                        sprintf('%d\n(%.0f%%)', complexity_counts(ki), 100*complexity_counts(ki)/n_total), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
                end
            end
            ylim([0 max(complexity_counts)*1.25]);
            
            % --- Panel 2: Variable Inclusion Rates ---
            subplot(1, 3, 2);
            main_effect_counts = [sum(is_spd); sum(is_tf); sum(is_sf); sum(is_or); sum(has_int)];
            main_effect_labels = {'Speed', 'TF', 'SF', 'OR', 'Any Interact.'};
            main_effect_colors = [colors.speed; colors.tf; colors.sf; colors.or_color; colors.interaction];
            
            bh = barh(main_effect_counts, 'FaceColor', 'flat', 'EdgeColor', 'none');
            bh.CData = main_effect_colors;
            set(gca, 'YTick', 1:5, 'YTickLabel', main_effect_labels, 'FontSize', 9);
            xlabel('# Neurons');
            
            max_main = max(main_effect_counts);
            for ki = 1:length(main_effect_counts)
                text(main_effect_counts(ki) + max_main*0.02, ki, ...
                    sprintf('%d (%.0f%%)', main_effect_counts(ki), 100*main_effect_counts(ki)/n_total), ...
                    'VerticalAlignment', 'middle', 'FontSize', 8);
            end
            xlim([0 max_main*1.35]);
            title(sprintf('Variable Inclusion (n=%d)', n_total), 'FontSize', 10);
            set(gca, 'box', 'off');
            
            % --- Panel 3: Interaction Breakdown ---
            subplot(1, 3, 3);
            int_counts = [sum(has_spd_tf); sum(has_spd_sf); sum(has_spd_or); ...
                          sum(has_tf_sf); sum(has_tf_or); sum(has_sf_or)];
            int_labels = {'Spd x TF', 'Spd x SF', 'Spd x OR', 'TF x SF', 'TF x OR', 'SF x OR'};
            int_colors = [0.9 0.2 0.2; 0.8 0.5 0.2; 0.6 0.2 0.6; 0.5 0.8 0.2; 0.2 0.5 0.8; 0.9 0.6 0.6];
            
            bh_int = barh(int_counts, 'FaceColor', 'flat', 'EdgeColor', 'none');
            bh_int.CData = int_colors;
            set(gca, 'YTick', 1:6, 'YTickLabel', int_labels, 'FontSize', 9);
            xlabel('# Neurons');
            
            max_int = max(max(int_counts), 1);
            for ki = 1:length(int_counts)
                text(int_counts(ki) + max_int*0.02, ki, ...
                    sprintf('%d (%.0f%%)', int_counts(ki), 100*int_counts(ki)/n_total), ...
                    'VerticalAlignment', 'middle', 'FontSize', 8);
            end
            xlim([0 max_int*1.4]);
            title('Interaction Breakdown', 'FontSize', 10);
            set(gca, 'box', 'off');
            
            sgtitle(sprintf('Forward Selection Results (threshold: \\Delta bps > %.3f)', threshold), ...
                'FontSize', 12, 'FontWeight', 'bold');
        end
        
        
        %% ================================================================
        %  Correlation Summary Plot
        %  ================================================================
        
        function fig = plot_correlation_summary(tc_corr_results, model_categories)
        %PLOT_CORRELATION_SUMMARY Generate tuning curve correlation box plots
        %   fig = glm_plotting.plot_correlation_summary(tc_corr_results, model_categories)
        %
        %   Creates a 3x4 figure (Models x Stimulus types) showing boxplots
        %   of Pearson correlations between observed and predicted tuning curves.
        %
        %   Inputs:
        %       tc_corr_results  - Struct with correlation fields (e.g., Selected_Speed_VT)
        %       model_categories - Cell array of unique model category strings
        %
        %   Output:
        %       fig - Figure handle
        
            cond_colors = glm_plotting.get_condition_colors();
            model_cat_markers = glm_plotting.get_model_category_markers();
            
            % Condition labels
            cond_labels = containers.Map();
            cond_labels('T_Vstatic') = 'T';
            cond_labels('V') = 'V';
            cond_labels('VT') = 'VT';
            
            fig = figure('Position', [50 50 1400 900], 'Visible', 'off');
            sgtitle('Tuning Curve Correlations: Observed vs Predicted', ...
                'FontSize', 14, 'FontWeight', 'bold');
            
            stimulus_types = {'Speed', 'TF', 'SF', 'OR'};
            model_names = {'Selected', 'Additive', 'FullInteraction'};
            
            % Valid conditions per stimulus type
            stim_conditions = containers.Map();
            stim_conditions('Speed') = {'T_Vstatic', 'VT'};
            stim_conditions('TF') = {'V', 'VT'};
            stim_conditions('SF') = {'V', 'VT'};
            stim_conditions('OR') = {'V', 'VT'};
            
            for row_i = 1:length(model_names)
                mn = model_names{row_i};
                
                for col_i = 1:length(stimulus_types)
                    stim_name = stimulus_types{col_i};
                    conds = stim_conditions(stim_name);
                    
                    ax = subplot(3, 4, (row_i-1)*4 + col_i);
                    hold on;
                    
                    x_tick_positions = [];
                    x_tick_labels_arr = {};
                    
                    for cond_i = 1:length(conds)
                        cond_name = conds{cond_i};
                        field_name = sprintf('%s_%s_%s', mn, stim_name, cond_name);
                        
                        if ~isfield(tc_corr_results, field_name), continue; end
                        
                        r_vals = tc_corr_results.(field_name);
                        valid_r = r_vals(~isnan(r_vals));
                        
                        if length(valid_r) < 2, continue; end
                        
                        x_pos = cond_i;
                        clr = cond_colors(cond_name);
                        
                        x_tick_positions = [x_tick_positions, x_pos]; %#ok<AGROW>
                        x_tick_labels_arr{end+1} = cond_labels(cond_name); %#ok<AGROW>
                        
                        % Draw boxplot
                        glm_plotting.draw_boxplot(ax, x_pos, r_vals, clr, 0.6);
                        
                        % Jittered scatter points
                        valid_idx = find(~isnan(r_vals));
                        jitter = 0.15 * randn(length(valid_idx), 1);
                        
                        if isfield(tc_corr_results, 'model_category') && ~isempty(model_categories)
                            for mci = 1:length(model_categories)
                                mc = model_categories{mci};
                                mc_mask = strcmp(tc_corr_results.model_category(valid_idx), mc);
                                if ~any(mc_mask), continue; end
                                
                                if model_cat_markers.isKey(mc)
                                    mkr = model_cat_markers(mc);
                                else
                                    mkr = 'o';
                                end
                                
                                mc_idx = valid_idx(mc_mask);
                                scatter(x_pos + jitter(mc_mask), r_vals(mc_idx), ...
                                    12, clr, mkr, 'MarkerFaceColor', 'none', ...
                                    'MarkerEdgeColor', clr, 'LineWidth', 0.5);
                            end
                        else
                            scatter(x_pos + jitter, r_vals(valid_idx), ...
                                12, clr, 'o', 'MarkerFaceColor', 'none', ...
                                'MarkerEdgeColor', clr, 'LineWidth', 0.5);
                        end
                    end
                    
                    % Format axes
                    if ~isempty(x_tick_positions)
                        set(gca, 'XTick', x_tick_positions, 'XTickLabel', x_tick_labels_arr);
                    end
                    xlim([0.3, length(conds)+0.7]);
                    ylim([-1, 1]);
                    plot(xlim, [0, 0], 'k--', 'LineWidth', 0.5);
                    
                    if col_i == 1
                        ylabel(sprintf('%s\nPearson r', mn), 'FontSize', 9);
                    end
                    if row_i == 1
                        title(stim_name, 'FontSize', 10);
                    end
                    set(gca, 'box', 'off', 'FontSize', 8);
                    hold off;
                end
            end
        end
        
        
        %% ================================================================
        %  Selected vs Additive Paired Comparison Plot
        %  ================================================================
        
        function fig = plot_selected_vs_additive(tc_corr_results, model_categories)
        %PLOT_SELECTED_VS_ADDITIVE Generate paired comparison scatter plots
        %   fig = glm_plotting.plot_selected_vs_additive(tc_corr_results, model_categories)
        %
        %   Creates a 2x2 figure comparing Selected vs Additive model performance.
        %   Connected lines show same-cluster comparisons across models.
        %
        %   Layout:
        %     Columns: 1) Pearson r, 2) R²
        %     Rows: 1) Speed tuning, 2) TF tuning
        %
        %   Inputs:
        %       tc_corr_results  - Struct with correlation fields
        %       model_categories - Cell array of unique model category strings
        %
        %   Output:
        %       fig - Figure handle
        
            cond_colors = glm_plotting.get_condition_colors();
            model_cat_markers = glm_plotting.get_model_category_markers();
            
            clr_T = cond_colors('T_Vstatic');
            clr_V = cond_colors('V');
            clr_VT = cond_colors('VT');
            
            fig = figure('Position', [50 50 1000 800], 'Visible', 'off');
            sgtitle('Selected vs Additive Model: Paired Comparison', ...
                'FontSize', 14, 'FontWeight', 'bold');
            
            x_selected = 1;
            x_additive = 2;
            
            % Subplot config: {stim, metric, field1, clr1, label1, field2, clr2, label2}
            subplot_config = {
                {'Speed', 'Pearson r', 'Speed_T_Vstatic', clr_T, 'T', 'Speed_VT', clr_VT, 'VT'};
                {'Speed', 'R²', 'Speed_T_Vstatic', clr_T, 'T', 'Speed_VT', clr_VT, 'VT'};
                {'TF', 'Pearson r', 'TF_V', clr_V, 'V', 'TF_VT', clr_VT, 'VT'};
                {'TF', 'R²', 'TF_V', clr_V, 'V', 'TF_VT', clr_VT, 'VT'}
            };
            
            for sp_i = 1:4
                ax = subplot(2, 2, sp_i);
                hold on;
                
                cfg = subplot_config{sp_i};
                stim_type = cfg{1};
                metric_type = cfg{2};
                field_cond1 = cfg{3};
                clr_cond1 = cfg{4};
                label_cond1 = cfg{5};
                field_cond2 = cfg{6};
                clr_cond2 = cfg{7};
                label_cond2 = cfg{8};
                
                is_r2 = strcmp(metric_type, 'R²');
                
                % Get data
                selected_cond1 = tc_corr_results.(sprintf('Selected_%s', field_cond1));
                additive_cond1 = tc_corr_results.(sprintf('Additive_%s', field_cond1));
                selected_cond2 = tc_corr_results.(sprintf('Selected_%s', field_cond2));
                additive_cond2 = tc_corr_results.(sprintf('Additive_%s', field_cond2));
                
                if is_r2
                    selected_cond1 = selected_cond1.^2;
                    additive_cond1 = additive_cond1.^2;
                    selected_cond2 = selected_cond2.^2;
                    additive_cond2 = additive_cond2.^2;
                end
                
                box_w = 0.18;
                scatter_offset = 0.22;
                
                % --- Condition 1 ---
                valid_c1 = ~isnan(selected_cond1) & ~isnan(additive_cond1);
                valid_idx_c1 = find(valid_c1);
                
                if sum(valid_c1) > 0
                    jitter_c1 = 0.03 * randn(sum(valid_c1), 1);
                    
                    % Connecting lines
                    for ii = 1:length(valid_idx_c1)
                        ci = valid_idx_c1(ii);
                        plot([x_selected - scatter_offset + jitter_c1(ii), x_additive + scatter_offset + jitter_c1(ii)], ...
                            [selected_cond1(ci), additive_cond1(ci)], ...
                            '-', 'Color', [clr_cond1, 0.15], 'LineWidth', 0.3);
                    end
                    
                    % Boxplots
                    glm_plotting.draw_boxplot_with_offset(ax, x_selected, selected_cond1(valid_c1), clr_cond1, box_w);
                    glm_plotting.draw_boxplot_with_offset(ax, x_additive, additive_cond1(valid_c1), clr_cond1, box_w);
                    
                    % Scatter points
                    glm_plotting.draw_model_category_scatter(ax, x_selected - scatter_offset, ...
                        selected_cond1, valid_idx_c1, jitter_c1, clr_cond1, ...
                        tc_corr_results, model_categories, model_cat_markers);
                    glm_plotting.draw_model_category_scatter(ax, x_additive + scatter_offset, ...
                        additive_cond1, valid_idx_c1, jitter_c1, clr_cond1, ...
                        tc_corr_results, model_categories, model_cat_markers);
                end
                
                % --- Condition 2 (VT) ---
                valid_c2 = ~isnan(selected_cond2) & ~isnan(additive_cond2);
                valid_idx_c2 = find(valid_c2);
                
                if sum(valid_c2) > 0
                    jitter_c2 = 0.03 * randn(sum(valid_c2), 1);
                    
                    % Connecting lines
                    for ii = 1:length(valid_idx_c2)
                        ci = valid_idx_c2(ii);
                        plot([x_selected - scatter_offset + jitter_c2(ii), x_additive + scatter_offset + jitter_c2(ii)], ...
                            [selected_cond2(ci), additive_cond2(ci)], ...
                            '-', 'Color', [clr_cond2, 0.15], 'LineWidth', 0.3);
                    end
                    
                    % Boxplots
                    glm_plotting.draw_boxplot_with_offset(ax, x_selected, selected_cond2(valid_c2), clr_cond2, box_w);
                    glm_plotting.draw_boxplot_with_offset(ax, x_additive, additive_cond2(valid_c2), clr_cond2, box_w);
                    
                    % Scatter points
                    glm_plotting.draw_model_category_scatter(ax, x_selected - scatter_offset, ...
                        selected_cond2, valid_idx_c2, jitter_c2, clr_cond2, ...
                        tc_corr_results, model_categories, model_cat_markers);
                    glm_plotting.draw_model_category_scatter(ax, x_additive + scatter_offset, ...
                        additive_cond2, valid_idx_c2, jitter_c2, clr_cond2, ...
                        tc_corr_results, model_categories, model_cat_markers);
                end
                
                % Axes formatting
                set(gca, 'XTick', [x_selected, x_additive], 'XTickLabel', {'Selected', 'Additive'});
                xlim([0.5, 2.5]);
                if is_r2
                    ylabel('R²');
                    ylim([0, 1]);
                else
                    ylabel('Pearson r');
                    ylim([-1, 1]);
                    plot(xlim, [0, 0], 'k--', 'LineWidth', 0.5);
                end
                title(sprintf('%s Tuning: %s', stim_type, metric_type), 'FontSize', 10);
                set(gca, 'box', 'off', 'FontSize', 9);
                hold off;
            end
        end
        
        
        function draw_boxplot_with_offset(ax, x_pos, values, clr, box_w)
        %DRAW_BOXPLOT_WITH_OFFSET Draw boxplot for paired comparison
        %   Internal helper for plot_selected_vs_additive - uses draw_boxplot
            glm_plotting.draw_boxplot(ax, x_pos, values, clr, box_w);
        end
        
        
        function draw_model_category_scatter(ax, x_pos, values, valid_idx, jitter, clr, ...
                tc_corr_results, model_categories, model_cat_markers)
        %DRAW_MODEL_CATEGORY_SCATTER Draw scatter with model-specific markers
        %   Internal helper for plot_selected_vs_additive
            axes(ax); %#ok<LAXES>
            
            if isempty(model_categories) || ~isfield(tc_corr_results, 'model_category')
                scatter(x_pos + jitter, values(valid_idx), 15, clr, 'o', ...
                    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr, 'LineWidth', 0.6);
                return;
            end
            
            for mci = 1:length(model_categories)
                mc = model_categories{mci};
                mc_mask = strcmp(tc_corr_results.model_category(valid_idx), mc);
                if ~any(mc_mask), continue; end
                
                if model_cat_markers.isKey(mc)
                    mkr = model_cat_markers(mc);
                else
                    mkr = 'o';
                end
                
                mc_idx = valid_idx(mc_mask);
                scatter(x_pos + jitter(mc_mask), values(mc_idx), 15, clr, mkr, ...
                    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr, 'LineWidth', 0.6);
            end
        end
        
        
        %% ================================================================
        %  Tuning Curve Error Bar Helper
        %  ================================================================
        
        function plot_tuning_curve_errorbar(ax, bin_centers, means, errors, clr, varargin)
        %PLOT_TUNING_CURVE_ERRORBAR Plot tuning curve with error bars
        %   glm_plotting.plot_tuning_curve_errorbar(ax, bin_centers, means, errors, clr)
        %
        %   Inputs:
        %       ax          - Axes handle
        %       bin_centers - X values (bin centers)
        %       means       - Mean firing rates
        %       errors      - Error values (SEM or STD)
        %       clr         - RGB color [1x3]
        %       varargin    - Optional: 'MarkerSize' (default: 4)
        
            p = inputParser;
            addParameter(p, 'MarkerSize', 4);
            addParameter(p, 'LineWidth', 0.8);
            addParameter(p, 'CapSize', 3);
            parse(p, varargin{:});
            
            axes(ax); %#ok<LAXES>
            hold on;
            
            valid = ~isnan(means);
            errorbar(bin_centers(valid), means(valid), errors(valid), 'o', ...
                'Color', clr, 'MarkerFaceColor', clr, 'MarkerSize', p.Results.MarkerSize, ...
                'LineWidth', p.Results.LineWidth, 'CapSize', p.Results.CapSize);
            plot(bin_centers(valid), means(valid), '-', 'Color', clr, 'LineWidth', 1);
        end
        
        
        %% ================================================================
        %  Scatter Plot Coloring Helpers
        %  ================================================================
        
        function clr_vec = get_scatter_colors_by_speed(speed_vals, speed_range)
        %GET_SCATTER_COLORS_BY_SPEED Get color vector based on speed values
        %   clr_vec = glm_plotting.get_scatter_colors_by_speed(speed_vals, speed_range)
        %
        %   Returns normalized colormap indices for viridis-style coloring.
            speed_norm = (speed_vals - speed_range(1)) / (speed_range(2) - speed_range(1));
            speed_norm = max(0, min(1, speed_norm));
            clr_vec = speed_norm;
        end
        
        
        function clr_mat = get_scatter_colors_by_sf(sf_vals, sf_levels)
        %GET_SCATTER_COLORS_BY_SF Get RGB colors based on SF values
        %   clr_mat = glm_plotting.get_scatter_colors_by_sf(sf_vals, sf_levels)
        %
        %   Inputs:
        %       sf_vals   - Vector of SF values
        %       sf_levels - Reference SF levels [0.003, 0.006, 0.012]
        %
        %   Output:
        %       clr_mat - Nx3 RGB color matrix
            sf_cmap = glm_plotting.get_sf_colormap();
            n = length(sf_vals);
            clr_mat = zeros(n, 3);
            
            for i = 1:n
                [~, idx] = min(abs(sf_levels - sf_vals(i)));
                if idx <= size(sf_cmap, 1)
                    clr_mat(i, :) = sf_cmap(idx, :);
                else
                    clr_mat(i, :) = [0.5 0.5 0.5];  % Gray fallback
                end
            end
        end
        
        
        function clr_mat = get_scatter_colors_by_or(or_vals, or_levels)
        %GET_SCATTER_COLORS_BY_OR Get RGB colors based on orientation values
        %   clr_mat = glm_plotting.get_scatter_colors_by_or(or_vals, or_levels)
        %
        %   Inputs:
        %       or_vals   - Vector of orientation values (radians)
        %       or_levels - Reference OR levels [-pi/4, 0, pi/4, pi/2]
        %
        %   Output:
        %       clr_mat - Nx3 RGB color matrix
            or_cmap = glm_plotting.get_or_colormap();
            n = length(or_vals);
            clr_mat = zeros(n, 3);
            
            for i = 1:n
                if isnan(or_vals(i))
                    clr_mat(i, :) = [0.5 0.5 0.5];  % Gray for NaN
                    continue;
                end
                [~, idx] = min(abs(or_levels - or_vals(i)));
                if idx <= size(or_cmap, 1)
                    clr_mat(i, :) = or_cmap(idx, :);
                else
                    clr_mat(i, :) = [0.5 0.5 0.5];
                end
            end
        end
        
        
        %% ================================================================
        %  Figure Saving Helper
        %  ================================================================
        
        function save_figure(fig, filepath, varargin)
        %SAVE_FIGURE Save figure to file with standard settings
        %   glm_plotting.save_figure(fig, filepath)
        %   glm_plotting.save_figure(fig, filepath, 'Resolution', 300)
        %
        %   Inputs:
        %       fig      - Figure handle
        %       filepath - Full path to save (extension determines format)
        %       varargin - Optional: 'Resolution' (default: 300)
        
            p = inputParser;
            addParameter(p, 'Resolution', 300);
            parse(p, varargin{:});
            
            [~, ~, ext] = fileparts(filepath);
            
            switch lower(ext)
                case '.png'
                    print(fig, filepath, '-dpng', sprintf('-r%d', p.Results.Resolution));
                case '.pdf'
                    print(fig, filepath, '-dpdf', '-bestfit');
                case '.svg'
                    print(fig, filepath, '-dsvg');
                case '.fig'
                    savefig(fig, filepath);
                otherwise
                    print(fig, filepath, '-dpng', sprintf('-r%d', p.Results.Resolution));
            end
        end
        
        
        %% ================================================================
        %  Classification String Builder
        %  ================================================================
        
        function class_str = build_classification_string(results, ci)
        %BUILD_CLASSIFICATION_STRING Build a classification label from results
        %   class_str = glm_plotting.build_classification_string(results, ci)
        %
        %   Builds strings like "Spd+TF +Int" based on forward selection results.
        %
        %   Inputs:
        %       results - Results struct with time.is_speed_tuned, etc.
        %       ci      - Cluster index
        %
        %   Output:
        %       class_str - Formatted classification string
        
            class_str = '';
            
            if results.time.is_speed_tuned(ci)
                class_str = 'Spd';
            end
            if results.time.is_tf_tuned(ci)
                if ~isempty(class_str), class_str = [class_str, '+']; end %#ok<AGROW>
                class_str = [class_str, 'TF']; %#ok<AGROW>
            end
            if results.time.is_sf_tuned(ci)
                if ~isempty(class_str), class_str = [class_str, '+']; end %#ok<AGROW>
                class_str = [class_str, 'SF']; %#ok<AGROW>
            end
            if results.time.is_or_tuned(ci)
                if ~isempty(class_str), class_str = [class_str, '+']; end %#ok<AGROW>
                class_str = [class_str, 'OR']; %#ok<AGROW>
            end
            if isempty(class_str)
                class_str = 'None';
            end
            if results.time.has_interaction(ci)
                class_str = [class_str, ' +Int']; %#ok<AGROW>
            end
        end
        
        
        %% ================================================================
        %  Model Abbreviations
        %  ================================================================
        
        function abbrev = get_model_abbreviations()
        %GET_MODEL_ABBREVIATIONS Return abbreviated model labels for tight spacing
        %   abbrev = glm_plotting.get_model_abbreviations()
        %
        %   Returns a containers.Map with model names as keys
            abbrev = containers.Map();
            abbrev('Null') = 'Null';
            abbrev('Selected') = 'Selected';
            abbrev('Additive') = 'Additive';
            abbrev('FullInteraction') = 'FullInt';
            abbrev('M0') = 'M0';
            abbrev('M0_Speed') = 'M0+S';
            abbrev('M0_Speed_TF') = 'M0+S+TF';
            abbrev('M0_Speed_TF_SF') = 'M0+S+TF+SF';
        end
        
        
        %% ================================================================
        %  Cross-Profile CV Performance Comparison Plot
        %  ================================================================
        
        function fig = plot_cross_profile_cv_comparison(standard_delta_bps, cross_profile_delta_bps)
        %PLOT_CROSS_PROFILE_CV_COMPARISON Generates scatter/boxplot of CV performance
        %   fig = glm_plotting.plot_cross_profile_cv_comparison(standard_bps, cross_profile_bps)
        %
        %   Visualizes the difference in delta bits/spike (Speed vs Null model)
        %   between standard k-fold CV and cross-profile CV (OOD).
        %
        %   Inputs:
        %       standard_delta_bps      - Array of delta bits/spike under standard CV
        %       cross_profile_delta_bps - Array of delta bits/spike under cross-profile CV
        %
        %   Output:
        %       fig - Figure handle
        
            standard_delta_bps = standard_delta_bps(:);
            cross_profile_delta_bps = cross_profile_delta_bps(:);
            valid = isfinite(standard_delta_bps) & isfinite(cross_profile_delta_bps);
            standard_delta_bps = standard_delta_bps(valid);
            cross_profile_delta_bps = cross_profile_delta_bps(valid);

            fig = figure('Name', 'Cross-Profile CV Comparison', 'Color', 'w', ...
                         'Position', [100, 100, 500, 600]);
            hold on;
            
            % Plot styling
            colors = glm_plotting.get_main_effect_colors();
            speed_color = colors.speed;
            
            x_standard = 1;
            x_cross = 2;
            
            % Draw connected lines for each neuron
            n_cells = length(standard_delta_bps);
            for i = 1:n_cells
                plot([x_standard, x_cross], [standard_delta_bps(i), cross_profile_delta_bps(i)], ...
                     '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            end
            
            % Draw points
            scatter(x_standard * ones(n_cells, 1), standard_delta_bps, 30, speed_color, 'filled', 'MarkerEdgeColor', 'w');
            scatter(x_cross * ones(n_cells, 1), cross_profile_delta_bps, 30, speed_color, 'filled', 'MarkerEdgeColor', 'w');
            
            % Draw boxplots via helper
            glm_plotting.draw_boxplot(gca, x_standard, standard_delta_bps, glm_plotting.darken_color(speed_color, 0.5), 0.4);
            glm_plotting.draw_boxplot(gca, x_cross, cross_profile_delta_bps, glm_plotting.darken_color(speed_color, 0.5), 0.4);
            
            % Formatting
            set(gca, 'XTick', [x_standard, x_cross], ...
                     'XTickLabel', {'Standard CV', 'Cross-Profile CV'}, ...
                     'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
            ylabel('\Delta Bits/Spike (Speed - Null)', 'FontSize', 14);
            title('Speed Tuning Generalization', 'FontSize', 14);
            xlim([0.5, 2.5]);
            
            % Significance test
            if n_cells >= 2
                p_val = signrank(standard_delta_bps, cross_profile_delta_bps);
            else
                p_val = NaN;
            end
            
            % Show p-value
            y_max = max([max(standard_delta_bps); max(cross_profile_delta_bps); 0.01]);
            if isnan(p_val)
                p_txt = 'p = n/a';
            else
                p_txt = sprintf('p = %.3e', p_val);
            end
            text(1.5, y_max * 1.05, sprintf('%s\nWilcoxon signed rank', p_txt), ...
                 'HorizontalAlignment', 'center', 'FontSize', 10);
             
            % Threshold line
            yline(0.005, '--k', 'Threshold (0.005)', 'LabelHorizontalAlignment', 'left');
            yline(0, '-k');
        end
        
    end
end
