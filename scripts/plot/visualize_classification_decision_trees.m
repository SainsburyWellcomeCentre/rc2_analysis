% VISUALIZE_CLASSIFICATION_DECISION_TREES
% 
% Creates a single unified tree diagram visualization of the complete
% classification pipeline used in glm_single_cluster_analysis.m:
%   1. Pre-filtering Decision Tree
%   2. Drop-one Classification (as table)
%   3. Additive vs FullInteraction test
%   4. Tuning Shape Classification for Speed
%
% Input: prefilter_decision_tree.csv from GLM analysis
%
% Author: Generated for rc2_analysis
% Date: February 2026

%% ====================================================================
%  Configuration
%  ====================================================================
close all; clearvars;

ctl = RC2Analysis();
ctl.setup_figures({'glm_single_cluster'}, true);

% Path to the results CSV
results_csv = fullfile(ctl.figs.curr_dir, 'prefilter_decision_tree.csv');

% Output directory for figures
output_dir = ctl.figs.curr_dir;

% Figure settings
save_figs = true;

%% ====================================================================
%  Load Data
%  ====================================================================
fprintf('Loading data from: %s\n', results_csv);

if exist(results_csv, 'file')
    data = readtable(results_csv);
    fprintf('Loaded %d clusters\n', height(data));
else
    error('Results CSV not found: %s', results_csv);
end

%% ====================================================================
%  Calculate Statistics
%  ====================================================================

% Pre-filtering statistics
n_total = height(data);
n_not_tuned = sum(strcmp(data.category, 'not_tuned'));
n_speed_tuned_only = sum(contains(data.category, 'speed_tuned'));
n_visually_tuned = sum(strcmp(data.category, 'visually_tuned'));
n_mixed_tuning = sum(strcmp(data.category, 'mixed_tuning'));
n_run_glm = sum(data.should_run_glm == 1);
n_no_glm = n_total - n_run_glm;
n_speed_VT_not_resp = sum(strcmp(data.category, 'speed_tuned_VT_not_responsive'));
n_speed_only = sum(strcmp(data.category, 'speed_tuned_only'));

% Drop-one classification statistics
has_glm = ~isnan(data.glm_delta_bps_drop_speed);
n_glm_analyzed = sum(has_glm);
n_speed_tuned_glm = sum(data.glm_is_speed_tuned == 1, 'omitnan');
n_tf_tuned_glm = sum(data.glm_is_tf_tuned == 1, 'omitnan');
n_sf_tuned_glm = sum(data.glm_is_sf_tuned == 1, 'omitnan');
n_or_tuned_glm = sum(data.glm_is_or_tuned == 1, 'omitnan');

% Interaction test statistics
n_has_interaction = sum(data.glm_has_significant_interaction == 1, 'omitnan');
n_no_interaction = n_glm_analyzed - n_has_interaction;

% Positive control: Visual contribution (Additive vs M0_Speed)
if any(strcmp(data.Properties.VariableNames, 'glm_is_visually_tuned_glm'))
    is_vis_glm = data.glm_is_visually_tuned_glm == 1;
    is_vis_glm(isnan(is_vis_glm)) = false;
    n_visually_tuned_glm = sum(is_vis_glm & has_glm);
else
    is_vis_glm = false(height(data), 1);
    n_visually_tuned_glm = 0;
end

% Speed vs Visual 2-way Venn (positive controls)
% Speed: Additive - Additive_no_Speed > threshold (unique speed contribution)
% Visual: Additive - M0_Speed > threshold (combined visual contribution)
is_spd_unique = data.glm_is_speed_tuned == 1;
is_spd_unique(isnan(is_spd_unique)) = false;

venn2_counts = struct();
venn2_counts.spd_only = sum(has_glm & is_spd_unique & ~is_vis_glm);
venn2_counts.vis_only = sum(has_glm & ~is_spd_unique & is_vis_glm);
venn2_counts.both = sum(has_glm & is_spd_unique & is_vis_glm);
venn2_counts.neither = sum(has_glm & ~is_spd_unique & ~is_vis_glm);

% Cluster IDs for 2-way Venn
get_cluster_ids_local = @(mask) strjoin(arrayfun(@(i) sprintf('%s_%d', data.probe_id{i}, data.cluster_id(i)), ...
    find(mask), 'UniformOutput', false), ', ');

venn2_clusters = struct();
venn2_clusters.spd_only = get_cluster_ids_local(has_glm & is_spd_unique & ~is_vis_glm);
venn2_clusters.vis_only = get_cluster_ids_local(has_glm & ~is_spd_unique & is_vis_glm);
venn2_clusters.both = get_cluster_ids_local(has_glm & is_spd_unique & is_vis_glm);
venn2_clusters.neither = get_cluster_ids_local(has_glm & ~is_spd_unique & ~is_vis_glm);

% SVM pre-filter tuning (from the decision tree)
n_svm_speed_tuned = sum(data.is_speed_tuned_prefilter == 1, 'omitnan');
n_svm_visually_tuned = sum(data.is_visually_tuned == 1, 'omitnan');
n_svm_VT_responsive = sum(data.is_VT_responsive == 1, 'omitnan');

% Calculate Venn diagram intersections for GLM drop-one tuning
% Create boolean vectors for each tuning type (GLM clusters only)
glm_idx = has_glm;
is_spd = data.glm_is_speed_tuned == 1;
is_tf = data.glm_is_tf_tuned == 1;
is_sf = data.glm_is_sf_tuned == 1;
is_or = data.glm_is_or_tuned == 1;

% Replace NaN with false
is_spd(isnan(is_spd)) = false;
is_tf(isnan(is_tf)) = false;
is_sf(isnan(is_sf)) = false;
is_or(isnan(is_or)) = false;

% Count all 16 possible combinations (2^4)
venn_counts = struct();
venn_counts.none = sum(glm_idx & ~is_spd & ~is_tf & ~is_sf & ~is_or);
venn_counts.spd_only = sum(glm_idx & is_spd & ~is_tf & ~is_sf & ~is_or);
venn_counts.tf_only = sum(glm_idx & ~is_spd & is_tf & ~is_sf & ~is_or);
venn_counts.sf_only = sum(glm_idx & ~is_spd & ~is_tf & is_sf & ~is_or);
venn_counts.or_only = sum(glm_idx & ~is_spd & ~is_tf & ~is_sf & is_or);
venn_counts.spd_tf = sum(glm_idx & is_spd & is_tf & ~is_sf & ~is_or);
venn_counts.spd_sf = sum(glm_idx & is_spd & ~is_tf & is_sf & ~is_or);
venn_counts.spd_or = sum(glm_idx & is_spd & ~is_tf & ~is_sf & is_or);
venn_counts.tf_sf = sum(glm_idx & ~is_spd & is_tf & is_sf & ~is_or);
venn_counts.tf_or = sum(glm_idx & ~is_spd & is_tf & ~is_sf & is_or);
venn_counts.sf_or = sum(glm_idx & ~is_spd & ~is_tf & is_sf & is_or);
venn_counts.spd_tf_sf = sum(glm_idx & is_spd & is_tf & is_sf & ~is_or);
venn_counts.spd_tf_or = sum(glm_idx & is_spd & is_tf & ~is_sf & is_or);
venn_counts.spd_sf_or = sum(glm_idx & is_spd & ~is_tf & is_sf & is_or);
venn_counts.tf_sf_or = sum(glm_idx & ~is_spd & is_tf & is_sf & is_or);
venn_counts.all_four = sum(glm_idx & is_spd & is_tf & is_sf & is_or);

% Collect cluster identities for each intersection
% Helper function to get cluster IDs as strings
get_cluster_ids = @(mask) strjoin(arrayfun(@(i) sprintf('%s_%d', data.probe_id{i}, data.cluster_id(i)), ...
    find(mask), 'UniformOutput', false), ', ');

venn_clusters = struct();
venn_clusters.none = get_cluster_ids(glm_idx & ~is_spd & ~is_tf & ~is_sf & ~is_or);
venn_clusters.spd_only = get_cluster_ids(glm_idx & is_spd & ~is_tf & ~is_sf & ~is_or);
venn_clusters.tf_only = get_cluster_ids(glm_idx & ~is_spd & is_tf & ~is_sf & ~is_or);
venn_clusters.sf_only = get_cluster_ids(glm_idx & ~is_spd & ~is_tf & is_sf & ~is_or);
venn_clusters.or_only = get_cluster_ids(glm_idx & ~is_spd & ~is_tf & ~is_sf & is_or);
venn_clusters.spd_tf = get_cluster_ids(glm_idx & is_spd & is_tf & ~is_sf & ~is_or);
venn_clusters.spd_sf = get_cluster_ids(glm_idx & is_spd & ~is_tf & is_sf & ~is_or);
venn_clusters.spd_or = get_cluster_ids(glm_idx & is_spd & ~is_tf & ~is_sf & is_or);
venn_clusters.tf_sf = get_cluster_ids(glm_idx & ~is_spd & is_tf & is_sf & ~is_or);
venn_clusters.tf_or = get_cluster_ids(glm_idx & ~is_spd & is_tf & ~is_sf & is_or);
venn_clusters.sf_or = get_cluster_ids(glm_idx & ~is_spd & ~is_tf & is_sf & is_or);
venn_clusters.spd_tf_sf = get_cluster_ids(glm_idx & is_spd & is_tf & is_sf & ~is_or);
venn_clusters.spd_tf_or = get_cluster_ids(glm_idx & is_spd & is_tf & ~is_sf & is_or);
venn_clusters.spd_sf_or = get_cluster_ids(glm_idx & is_spd & ~is_tf & is_sf & is_or);
venn_clusters.tf_sf_or = get_cluster_ids(glm_idx & ~is_spd & is_tf & is_sf & is_or);
venn_clusters.all_four = get_cluster_ids(glm_idx & is_spd & is_tf & is_sf & is_or);

% Tuning shape statistics (from speed_shape_T column - pre-computed from tuning curves)
if any(strcmp(data.Properties.VariableNames, 'speed_shape_T'))
    shapes = data.speed_shape_T;
    n_mono_inc = sum(strcmp(shapes, 'increasing'));
    n_mono_dec = sum(strcmp(shapes, 'decreasing'));
    n_bandpass = sum(strcmp(shapes, 'bandpass'));
    n_unclassified = sum(strcmp(shapes, 'unclassified') | cellfun(@isempty, shapes));
else
    n_mono_inc = 0; n_mono_dec = 0; n_bandpass = 0; n_unclassified = n_total;
end

% Also count VT shapes for comparison
if any(strcmp(data.Properties.VariableNames, 'speed_shape_VT'))
    shapes_VT = data.speed_shape_VT;
    n_mono_inc_VT = sum(strcmp(shapes_VT, 'increasing'));
    n_mono_dec_VT = sum(strcmp(shapes_VT, 'decreasing'));
    n_bandpass_VT = sum(strcmp(shapes_VT, 'bandpass'));
else
    n_mono_inc_VT = 0; n_mono_dec_VT = 0; n_bandpass_VT = 0;
end

%% ====================================================================
%  Create Single Unified Decision Tree Figure
%  ====================================================================
fprintf('\n--- Creating Unified Classification Decision Tree ---\n');

fig = figure('Name', 'Classification Decision Tree', ...
    'Position', [50 50 2000 1400], 'Color', 'w');

ax = axes('Position', [0.02 0.02 0.96 0.92]);
hold on;

% Define coordinate system
x_range = [0 40];
y_range = [-8 24];
axis([x_range y_range]);
axis off;
set(gca, 'DataAspectRatio', [1.3 1 1]);

% Box dimensions
box_w = 4.0;
box_h = 1.4;
box_h_small = 1.1;

% ========== SECTION 1: PRE-FILTERING DECISION TREE ==========

% Section title
text(9, 23, 'PRE-FILTERING DECISION TREE', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Level 0: Start
draw_box(9, 21.5, box_w, box_h, sprintf('All VISp Clusters\n(n=%d)', n_total));

% Level 1: Speed tuning check
draw_arrow_connected(9, 21.5 - box_h/2, 9, 19 + box_h/2);
draw_box(9, 19, box_w + 0.8, box_h, sprintf('Speed Tuned? (T)\n(stat_T ≠ mot_T)'));

% Branch: Not speed tuned (left)
draw_arrow_connected(9 - box_w/2 - 0.4, 19, 4, 16.5 + box_h/2);
text(5.5, 18, 'No', 'FontSize', 10, 'FontWeight', 'bold');

% Branch: Speed tuned (right)
draw_arrow_connected(9 + box_w/2 + 0.4, 19, 14, 16.5 + box_h/2);
text(12, 18, 'Yes', 'FontSize', 10, 'FontWeight', 'bold');

% Level 2a: Visual tuning (not speed tuned branch)
draw_box(4, 16.5, box_w + 0.5, box_h, sprintf('Visually Tuned? (V)\n(stat_V ≠ mot_V)'));

% Level 2b: Visual tuning (speed tuned branch)
draw_box(14, 16.5, box_w + 0.5, box_h, sprintf('Visually Tuned? (V)\n(stat_V ≠ mot_V)'));

% Outcomes from left branch (not speed tuned)
% No visual -> NOT TUNED
draw_arrow_connected(4 - box_w/2 - 0.25, 16.5, 1.5, 14 + box_h_small/2);
text(2, 15.5, 'No', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(1.5, 14, box_w, box_h_small, sprintf('NOT TUNED\n(n=%d)', n_not_tuned));
draw_terminal_mark(1.5, 14 - box_h_small/2 - 0.4, 'X');

% Yes visual -> VISUALLY TUNED
draw_arrow_connected(4 + box_w/2 + 0.25, 16.5, 6.5, 14 + box_h_small/2);
text(6, 15.5, 'Yes', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(6.5, 14, box_w, box_h_small, sprintf('VISUALLY\nTUNED'));

% Outcomes from right branch (speed tuned)
% Yes visual -> VISUALLY TUNED
draw_arrow_connected(14 + box_w/2 + 0.25, 16.5, 18, 14 + box_h_small/2);
text(16.5, 15.5, 'Yes', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(18, 14, box_w, box_h_small, sprintf('VISUALLY TUNED\n(n=%d)', n_visually_tuned));

% No visual -> VT responsive check
draw_arrow_connected(14 - box_w/2 - 0.25, 16.5, 11, 14 + box_h/2);
text(11.5, 15.5, 'No', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(11, 14, box_w + 0.8, box_h, sprintf('VT Responsive?\n(stat_VT ≠ mot_VT)'));

% VT responsive outcomes
% No -> SPEED ONLY (VT not responsive)
draw_arrow_connected(11 - box_w/2 - 0.4, 14, 7.5, 11.5 + box_h_small/2);
text(8, 13, 'No', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(7.5, 11.5, box_w, box_h_small, sprintf('SPEED ONLY\n(n=%d)', n_speed_VT_not_resp));
draw_terminal_mark(7.5, 11.5 - box_h_small/2 - 0.4, 'X');

% Yes -> Compare T vs VT
draw_arrow_connected(11 + box_w/2 + 0.4, 14, 14.5, 11.5 + box_h/2);
text(13, 13, 'Yes', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(14.5, 11.5, box_w + 0.8, box_h, sprintf('mot_T = mot_VT?\n(Kruskal-Wallis)'));

% T vs VT outcomes
% Same -> SPEED ONLY
draw_arrow_connected(14.5 - box_w/2 - 0.4, 11.5, 11, 9 + box_h_small/2);
text(11.5, 10.5, 'Same', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(11, 9, box_w, box_h_small, sprintf('SPEED ONLY\n(n=%d)', n_speed_only));
draw_terminal_mark(11, 9 - box_h_small/2 - 0.4, 'X');

% Different -> MIXED TUNING -> GLM
draw_arrow_connected(14.5 + box_w/2 + 0.4, 11.5, 18, 9 + box_h_small/2);
text(16.5, 10.5, 'Diff', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(18, 9, box_w, box_h_small, sprintf('MIXED TUNING\n(n=%d)', n_mixed_tuning));

% ========== Connect to GLM ==========
% Draw merge lines from all paths that go to GLM
draw_line_only(6.5, 14 - box_h_small/2, 6.5, 5);
draw_line_only(18, 14 - box_h_small/2, 18, 6);
draw_line_only(18, 9 - box_h_small/2, 18, 6);  % Mixed tuning connects here

% Horizontal merge line
line([6.5, 18], [6, 6], 'Color', 'k', 'LineWidth', 1.2);
draw_arrow_connected(12, 6, 12, 4.5 + box_h/2);

% GLM Analysis box
draw_box(12, 4.5, box_w + 0.8, box_h, sprintf('GLM Analysis\n(n=%d)', n_run_glm));

% ========== SECTION 2: TUNING SHAPE (from Pre-filtering, for ALL speed-tuned) ==========
% Speed tuning shape is computed during pre-filtering for ALL speed-tuned clusters
text(3.5, 7.5, 'TUNING SHAPE (Pre-filter)', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Arrow from speed-tuned branches to shape classification
% This applies to all speed-tuned clusters (including those excluded from GLM)
draw_line_only(14, 16.5 - box_h/2, 14, 8);
draw_line_only(14, 8, 3.5, 8);
draw_arrow_connected(3.5, 8, 3.5, 6.5 + box_h/2);

draw_box(3.5, 6.5, box_w + 0.6, box_h, sprintf('All Speed Tuned\n(n=%d)', n_svm_speed_tuned));

% Shape classification
draw_arrow_connected(3.5, 6.5 - box_h/2, 3.5, 4 + box_h/2);
draw_box(3.5, 4, box_w + 0.6, box_h, sprintf('Asymm Gaussian\nFit to T curve'));

draw_arrow_connected(3.5, 4 - box_h/2, 3.5, 1.5 + box_h_small/2);
draw_box(3.5, 1.5, box_w, box_h_small, 'Peak Location?');

% Shape outcomes
draw_arrow_connected(3.5 - box_w/2, 1.5 - box_h_small/2, 1, -0.5);
text(1, -1.2, sprintf('Decr.\n(n=%d)', n_mono_dec), 'FontSize', 9, 'HorizontalAlignment', 'center');

draw_arrow_connected(3.5, 1.5 - box_h_small/2, 3.5, -0.5);
text(3.5, -1.2, sprintf('Bandpass\n(n=%d)', n_bandpass), 'FontSize', 9, 'HorizontalAlignment', 'center');

draw_arrow_connected(3.5 + box_w/2, 1.5 - box_h_small/2, 6, -0.5);
text(6, -1.2, sprintf('Incr.\n(n=%d)', n_mono_inc), 'FontSize', 9, 'HorizontalAlignment', 'center');

% Note about VT shapes
text(3.5, -2.2, sprintf('(VT: Inc=%d, Dec=%d, Band=%d)', n_mono_inc_VT, n_mono_dec_VT, n_bandpass_VT), ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');

% ========== SECTION 2b: SPEED vs VISUAL 2-WAY VENN (Positive Controls) ==========
text(11, -3.5, 'POSITIVE CONTROLS (Speed vs Visual)', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% 2-circle Venn diagram
venn2_cx = 11;
venn2_cy = -6;
venn2_r = 2.2;  % radius
venn2_offset = 1.4;  % overlap offset

% Draw circles
theta_circ = linspace(0, 2*pi, 100);
% Speed circle (left, red)
x_spd = venn2_r * cos(theta_circ) + venn2_cx - venn2_offset;
y_spd = venn2_r * sin(theta_circ) + venn2_cy;
patch(x_spd, y_spd, [0.8 0.2 0.2], 'FaceAlpha', 0.25, 'EdgeColor', [0.7 0.1 0.1], 'LineWidth', 1.5);

% Visual circle (right, blue)
x_vis = venn2_r * cos(theta_circ) + venn2_cx + venn2_offset;
y_vis = venn2_r * sin(theta_circ) + venn2_cy;
patch(x_vis, y_vis, [0.2 0.2 0.8], 'FaceAlpha', 0.25, 'EdgeColor', [0.1 0.1 0.7], 'LineWidth', 1.5);

% Labels
text(venn2_cx - venn2_offset - venn2_r - 0.5, venn2_cy + 1.8, sprintf('Speed\n(Add - Add_{no\\_Spd})\nn=%d', n_speed_tuned_glm), ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1], 'HorizontalAlignment', 'center');
text(venn2_cx + venn2_offset + venn2_r + 0.5, venn2_cy + 1.8, sprintf('Visual\n(Add - M0_{Spd})\nn=%d', n_visually_tuned_glm), ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.7], 'HorizontalAlignment', 'center');

% Counts in regions
text(venn2_cx - venn2_offset - 1.0, venn2_cy, num2str(venn2_counts.spd_only), 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(venn2_cx, venn2_cy, num2str(venn2_counts.both), 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(venn2_cx + venn2_offset + 1.0, venn2_cy, num2str(venn2_counts.vis_only), 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Neither count (outside)
if venn2_counts.neither > 0
    text(venn2_cx, venn2_cy - venn2_r - 1.2, sprintf('Neither: %d', venn2_counts.neither), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');
end

% Explanation
text(venn2_cx, venn2_cy + venn2_r + 1.0, sprintf('(Δ > 0.01 bps threshold)'), ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');

% ========== SECTION 3: GLM MODEL COMPARISON ==========
text(30, 23, 'GLM MODEL SELECTION', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Arrow from GLM to model comparison
draw_arrow_connected(12 + box_w/2 + 0.4, 4.5, 22, 4.5);

% Additive vs FullInteraction test
draw_box(26, 4.5, box_w + 2.2, box_h, sprintf('Interaction Test\nFull - Additive > 0.01?'));

% Outcomes
draw_arrow_connected(26 - box_w/2 - 1.1, 4.5 - box_h/2, 22, 2 + box_h_small/2);
text(21, 3.3, 'No', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(22, 2, box_w, box_h_small, sprintf('Additive\n(n=%d)', n_no_interaction));

draw_arrow_connected(26 + box_w/2 + 1.1, 4.5 - box_h/2, 30, 2 + box_h_small/2);
text(30.5, 3.3, 'Yes', 'FontSize', 10, 'FontWeight', 'bold');
draw_box(30, 2, box_w, box_h_small, sprintf('FullInteraction\n(n=%d)', n_has_interaction));

% Model formulas
text(22, 0.5, 'β₀ + f(Spd) + g(TF) + SF + OR', 'FontSize', 9, 'HorizontalAlignment', 'center');
text(30, 0.5, '+ all interactions', 'FontSize', 9, 'HorizontalAlignment', 'center');

% ========== SECTION 4: DROP-ONE VENN DIAGRAM ==========
text(30, 21, 'DROP-ONE TUNING (Venn)', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Draw 4-set Venn diagram using overlapping ellipses
% Center of Venn diagram
venn_cx = 30;
venn_cy = 14;

% Ellipse parameters (tilted ellipses for 4-set Venn)
% Using the classic 4-ellipse Venn layout
ellipse_a = 3.8;  % semi-major axis
ellipse_b = 2.2;  % semi-minor axis

% Define ellipse positions and rotations (degrees)
% Speed (top-left), TF (top-right), SF (bottom-left), OR (bottom-right)
ellipse_params = struct();
ellipse_params.spd = struct('cx', venn_cx - 1.2, 'cy', venn_cy + 1.0, 'angle', 45);
ellipse_params.tf = struct('cx', venn_cx + 1.2, 'cy', venn_cy + 1.0, 'angle', -45);
ellipse_params.sf = struct('cx', venn_cx - 1.2, 'cy', venn_cy - 1.0, 'angle', -45);
ellipse_params.or = struct('cx', venn_cx + 1.2, 'cy', venn_cy - 1.0, 'angle', 45);

% Colors for each set (transparent)
colors = struct();
colors.spd = [0.8 0.2 0.2 0.2];  % red
colors.tf = [0.2 0.6 0.2 0.2];   % green
colors.sf = [0.2 0.2 0.8 0.2];   % blue
colors.or = [0.8 0.6 0.0 0.2];   % orange

% Draw ellipses
draw_rotated_ellipse(ellipse_params.spd.cx, ellipse_params.spd.cy, ellipse_a, ellipse_b, ellipse_params.spd.angle, colors.spd);
draw_rotated_ellipse(ellipse_params.tf.cx, ellipse_params.tf.cy, ellipse_a, ellipse_b, ellipse_params.tf.angle, colors.tf);
draw_rotated_ellipse(ellipse_params.sf.cx, ellipse_params.sf.cy, ellipse_a, ellipse_b, ellipse_params.sf.angle, colors.sf);
draw_rotated_ellipse(ellipse_params.or.cx, ellipse_params.or.cy, ellipse_a, ellipse_b, ellipse_params.or.angle, colors.or);

% Labels for each set (outside the ellipses)
text(venn_cx - 5.5, venn_cy + 3, sprintf('Speed\n(n=%d)', n_speed_tuned_glm), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1], 'HorizontalAlignment', 'center');
text(venn_cx + 5.5, venn_cy + 3, sprintf('TF\n(n=%d)', n_tf_tuned_glm), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.1 0.5 0.1], 'HorizontalAlignment', 'center');
text(venn_cx - 5.5, venn_cy - 3, sprintf('SF\n(n=%d)', n_sf_tuned_glm), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.7], 'HorizontalAlignment', 'center');
text(venn_cx + 5.5, venn_cy - 3, sprintf('OR\n(n=%d)', n_or_tuned_glm), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.7 0.5 0.0], 'HorizontalAlignment', 'center');

% Place counts in the Venn regions (approximate positions)
% Single sets only
if venn_counts.spd_only > 0
    text(venn_cx - 3.2, venn_cy + 2.4, num2str(venn_counts.spd_only), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.tf_only > 0
    text(venn_cx + 3.2, venn_cy + 2.4, num2str(venn_counts.tf_only), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.sf_only > 0
    text(venn_cx - 3.2, venn_cy - 2.4, num2str(venn_counts.sf_only), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.or_only > 0
    text(venn_cx + 3.2, venn_cy - 2.4, num2str(venn_counts.or_only), 'FontSize', 9, 'HorizontalAlignment', 'center');
end

% Two-way intersections
if venn_counts.spd_tf > 0
    text(venn_cx, venn_cy + 2.6, num2str(venn_counts.spd_tf), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.spd_sf > 0
    text(venn_cx - 2.6, venn_cy, num2str(venn_counts.spd_sf), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.spd_or > 0
    text(venn_cx - 1.0, venn_cy + 1.0, num2str(venn_counts.spd_or), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.tf_sf > 0
    text(venn_cx + 1.0, venn_cy - 1.0, num2str(venn_counts.tf_sf), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.tf_or > 0
    text(venn_cx + 2.6, venn_cy, num2str(venn_counts.tf_or), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.sf_or > 0
    text(venn_cx, venn_cy - 2.6, num2str(venn_counts.sf_or), 'FontSize', 9, 'HorizontalAlignment', 'center');
end

% Three-way intersections
if venn_counts.spd_tf_sf > 0
    text(venn_cx - 1.2, venn_cy + 0.4, num2str(venn_counts.spd_tf_sf), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.spd_tf_or > 0
    text(venn_cx + 0.4, venn_cy + 1.2, num2str(venn_counts.spd_tf_or), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.spd_sf_or > 0
    text(venn_cx - 0.4, venn_cy - 1.2, num2str(venn_counts.spd_sf_or), 'FontSize', 9, 'HorizontalAlignment', 'center');
end
if venn_counts.tf_sf_or > 0
    text(venn_cx + 1.2, venn_cy - 0.4, num2str(venn_counts.tf_sf_or), 'FontSize', 9, 'HorizontalAlignment', 'center');
end

% Four-way intersection (center)
if venn_counts.all_four > 0
    text(venn_cx, venn_cy, num2str(venn_counts.all_four), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Not tuned to any (outside)
if venn_counts.none > 0
    text(venn_cx, venn_cy - 5.5, sprintf('None: %d', venn_counts.none), 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');
end

% Formula explanation
text(venn_cx, venn_cy + 5.5, sprintf('(n=%d GLM analyzed, Δ > 0.01 bps)', n_glm_analyzed), ...
    'FontSize', 9, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');

% ========== LEGEND ==========
legend_x = 37;
legend_y = 21;
text(legend_x, legend_y + 1.5, 'LEGEND', 'FontSize', 11, 'FontWeight', 'bold');
draw_box(legend_x, legend_y, 2.2, 0.7, 'Node');
text(legend_x + 2, legend_y, '= Decision', 'FontSize', 10);
draw_terminal_mark(legend_x, legend_y - 1.5, 'X');
text(legend_x + 0.6, legend_y - 1.5, '= No GLM', 'FontSize', 10);
text(legend_x - 0.5, legend_y - 3, 'T = Translation (speed)', 'FontSize', 9);
text(legend_x - 0.5, legend_y - 3.9, 'V = Visual only', 'FontSize', 9);
text(legend_x - 0.5, legend_y - 4.8, 'VT = Visual + Translation', 'FontSize', 9);
text(legend_x - 0.5, legend_y - 6, 'stat = stationary', 'FontSize', 9);
text(legend_x - 0.5, legend_y - 6.9, 'mot = motion', 'FontSize', 9);

% ========== SUMMARY BOX ==========
summary_str = sprintf(['SUMMARY\n' ...
    '───────────────────\n' ...
    'Total clusters: %d\n' ...
    'Run GLM: %d (%.1f%%)\n' ...
    'Skip GLM: %d (%.1f%%)\n' ...
    '───────────────────\n' ...
    'Additive wins: %d\n' ...
    'FullInteraction: %d\n' ...
    '───────────────────\n' ...
    'Speed shapes (all):\n' ...
    '  Increasing: %d\n' ...
    '  Decreasing: %d\n' ...
    '  Bandpass: %d'], ...
    n_total, n_run_glm, 100*n_run_glm/n_total, n_no_glm, 100*n_no_glm/n_total, ...
    n_no_interaction, n_has_interaction, ...
    n_mono_inc, n_mono_dec, n_bandpass);

text(37, 11, summary_str, 'FontSize', 10, 'FontName', 'FixedWidth', ...
    'VerticalAlignment', 'top', 'BackgroundColor', [0.95 0.95 0.95], ...
    'EdgeColor', 'k', 'Margin', 6);

% ========== CLUSTER IDENTITY LIST ==========
% Build a compact list of cluster IDs for each non-empty intersection
cluster_list_lines = {'CLUSTER IDs BY TUNING', '─────────────────────────'};

% Speed vs Visual positive controls (2-way)
cluster_list_lines{end+1} = '--- Positive Controls ---';
if venn2_counts.spd_only > 0
    cluster_list_lines{end+1} = sprintf('Spd unique only: %s', truncate_str(venn2_clusters.spd_only, 55));
end
if venn2_counts.both > 0
    cluster_list_lines{end+1} = sprintf('Spd+Vis both: %s', truncate_str(venn2_clusters.both, 55));
end
if venn2_counts.vis_only > 0
    cluster_list_lines{end+1} = sprintf('Vis only: %s', truncate_str(venn2_clusters.vis_only, 55));
end
if venn2_counts.neither > 0
    cluster_list_lines{end+1} = sprintf('Neither: %s', truncate_str(venn2_clusters.neither, 55));
end

cluster_list_lines{end+1} = '--- Drop-one Details ---';

% Single tuning
if venn_counts.spd_only > 0
    cluster_list_lines{end+1} = sprintf('Spd only: %s', truncate_str(venn_clusters.spd_only, 60));
end
if venn_counts.tf_only > 0
    cluster_list_lines{end+1} = sprintf('TF only: %s', truncate_str(venn_clusters.tf_only, 60));
end
if venn_counts.sf_only > 0
    cluster_list_lines{end+1} = sprintf('SF only: %s', truncate_str(venn_clusters.sf_only, 60));
end
if venn_counts.or_only > 0
    cluster_list_lines{end+1} = sprintf('OR only: %s', truncate_str(venn_clusters.or_only, 60));
end

% Two-way combinations
if venn_counts.spd_tf > 0
    cluster_list_lines{end+1} = sprintf('Spd+TF: %s', truncate_str(venn_clusters.spd_tf, 60));
end
if venn_counts.spd_sf > 0
    cluster_list_lines{end+1} = sprintf('Spd+SF: %s', truncate_str(venn_clusters.spd_sf, 60));
end
if venn_counts.spd_or > 0
    cluster_list_lines{end+1} = sprintf('Spd+OR: %s', truncate_str(venn_clusters.spd_or, 60));
end
if venn_counts.tf_sf > 0
    cluster_list_lines{end+1} = sprintf('TF+SF: %s', truncate_str(venn_clusters.tf_sf, 60));
end
if venn_counts.tf_or > 0
    cluster_list_lines{end+1} = sprintf('TF+OR: %s', truncate_str(venn_clusters.tf_or, 60));
end
if venn_counts.sf_or > 0
    cluster_list_lines{end+1} = sprintf('SF+OR: %s', truncate_str(venn_clusters.sf_or, 60));
end

% Three-way combinations
if venn_counts.spd_tf_sf > 0
    cluster_list_lines{end+1} = sprintf('Spd+TF+SF: %s', truncate_str(venn_clusters.spd_tf_sf, 60));
end
if venn_counts.spd_tf_or > 0
    cluster_list_lines{end+1} = sprintf('Spd+TF+OR: %s', truncate_str(venn_clusters.spd_tf_or, 60));
end
if venn_counts.spd_sf_or > 0
    cluster_list_lines{end+1} = sprintf('Spd+SF+OR: %s', truncate_str(venn_clusters.spd_sf_or, 60));
end
if venn_counts.tf_sf_or > 0
    cluster_list_lines{end+1} = sprintf('TF+SF+OR: %s', truncate_str(venn_clusters.tf_sf_or, 60));
end

% Four-way
if venn_counts.all_four > 0
    cluster_list_lines{end+1} = sprintf('All four: %s', truncate_str(venn_clusters.all_four, 60));
end

% None
if venn_counts.none > 0
    cluster_list_lines{end+1} = sprintf('None: %s', truncate_str(venn_clusters.none, 60));
end

cluster_list_str = strjoin(cluster_list_lines, '\n');
text(22, -1, cluster_list_str, 'FontSize', 8, 'FontName', 'FixedWidth', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [0.98 0.98 0.98], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);

% Title
title('Complete GLM Classification Decision Tree', 'FontSize', 15, 'FontWeight', 'bold');

if save_figs
    % Save as PNG only
    print(fig, fullfile(output_dir, 'classification_decision_tree.png'), '-dpng', '-r300');
    fprintf('Saved: classification_decision_tree.png\n');
end

fprintf('\n=== Decision tree figure saved ===\n');

%% ====================================================================
%  Helper Functions
%  ====================================================================

function draw_box(x, y, w, h, txt)
    % Draw a black/white rectangular box with text
    % x, y: center position
    % w, h: width and height
    % txt: text to display
    
    rectangle('Position', [x - w/2, y - h/2, w, h], ...
        'FaceColor', 'w', ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.5, ...
        'Curvature', [0.15, 0.15]);
    
    text(x, y, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 9, ...
        'Interpreter', 'none');
end

function draw_table_cell(x, y, w, h, txt, is_header)
    % Draw a table cell
    if is_header
        facecolor = [0.85 0.85 0.85];
        fontweight = 'bold';
    else
        facecolor = 'w';
        fontweight = 'normal';
    end
    
    rectangle('Position', [x - w/2, y - h/2, w, h], ...
        'FaceColor', facecolor, ...
        'EdgeColor', 'k', ...
        'LineWidth', 1);
    
    text(x, y, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 9, ...
        'FontWeight', fontweight, ...
        'Interpreter', 'none');
end

function draw_arrow_connected(x1, y1, x2, y2)
    % Draw an arrow that connects exactly to the specified points
    
    % Arrow head parameters
    head_len = 0.25;
    head_width = 0.15;
    
    % Calculate direction
    dx = x2 - x1;
    dy = y2 - y1;
    len = sqrt(dx^2 + dy^2);
    
    if len < 0.01
        return;
    end
    
    % Unit vector
    ux = dx / len;
    uy = dy / len;
    
    % Arrow head base
    bx = x2 - head_len * ux;
    by = y2 - head_len * uy;
    
    % Perpendicular for head width
    px = -uy * head_width;
    py = ux * head_width;
    
    % Draw line
    line([x1, bx], [y1, by], 'Color', 'k', 'LineWidth', 1.2);
    
    % Draw arrow head
    patch([x2, bx + px, bx - px], [y2, by + py, by - py], 'k', ...
        'EdgeColor', 'k', 'LineWidth', 1);
end

function draw_line_only(x1, y1, x2, y2)
    % Draw a simple line without arrow
    line([x1, x2], [y1, y2], 'Color', 'k', 'LineWidth', 1.2);
end

function draw_terminal_mark(x, y, symbol)
    % Draw a terminal marker (X for stop)
    if strcmp(symbol, 'X')
        plot(x, y, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
    else
        plot(x, y, 'ko', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'k');
    end
end

function draw_rotated_ellipse(cx, cy, a, b, angle_deg, color)
    % Draw a rotated ellipse with transparency
    % cx, cy: center
    % a, b: semi-major and semi-minor axes
    % angle_deg: rotation angle in degrees
    % color: [R G B alpha] or [R G B]
    
    theta = linspace(0, 2*pi, 100);
    x = a * cos(theta);
    y = b * sin(theta);
    
    % Rotation matrix
    angle_rad = angle_deg * pi / 180;
    R = [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)];
    
    % Rotate points
    xy_rot = R * [x; y];
    x_rot = xy_rot(1, :) + cx;
    y_rot = xy_rot(2, :) + cy;
    
    % Draw with transparency
    if length(color) == 4
        patch(x_rot, y_rot, color(1:3), 'FaceAlpha', color(4), 'EdgeColor', color(1:3), 'LineWidth', 1.5);
    else
        patch(x_rot, y_rot, color, 'FaceAlpha', 0.2, 'EdgeColor', color, 'LineWidth', 1.5);
    end
end

function result = ternary(condition, true_val, false_val)
    % Simple ternary operator
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function result = truncate_str(s, n)
    % Truncate string s to n characters, adding '...' if truncated
    if length(s) > n
        result = [s(1:n) '...'];
    else
        result = s;
    end
end
