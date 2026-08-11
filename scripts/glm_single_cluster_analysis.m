% GLM-based single-cluster analysis of multisensory integration
%
% Version 10: Hierarchical Forward Selection
%
% This script first applies a decision tree to reduce the number of clusters
% that need full GLM analysis, then uses HIERARCHICAL FORWARD SELECTION to 
% determine which variables each neuron encodes.
%
% KEY CHANGE (v10): Two-phase selection ensures proper model hierarchy:
%   Phase 1: Select main effects (Speed, TF, SF, OR)
%   Phase 2: Select interactions only if BOTH parents were selected
%
% ONSET DYNAMICS KERNEL (from v9):
% The onset kernel uses raised cosine bases spanning 0-2s after motion onset,
% capturing onset transients and adaptation effects. The intercept (β₀) is
% anchored to stationary periods where both speed=0 AND onset_kernel=0.
%   - Stationary bins have time_since_onset < 0 → onset kernel = 0
%   - Motion bins have time_since_onset >= 0 → onset kernel captures transient
%
% PRE-FILTERING DECISION TREE:
% ============================
% For each cluster:
%   1. Check T significance: stationary_T != motion_T (Wilcoxon signrank)
%   2. Check V significance: stationary_V != motion_V (Wilcoxon signrank)
%   3. Check VT significance: stationary_VT != motion_VT (Wilcoxon signrank)
%   
%   SPEED TUNED CLASSIFICATION:
%       T significant OR VT significant -> "speed_tuned"
%   
%   GLM SELECTION:
%       VT significant -> run GLM
%       OR any combo of 2 significant (T+V, T+VT, V+VT) -> run GLM
%       Otherwise -> do not run GLM
%
% FORWARD SELECTION PROCEDURE (Hierarchical, inspired by Hardcastle et al. 2017):
% ================================================================================
%   PHASE 1 - Main Effects:
%     1. Start with Null model (intercept + onset kernel)
%     2. Test main effects only: Speed, TF, SF, OR
%     3. Add the best if Δ CV bits/spike > threshold (0.005)
%     4. Repeat until no main effect improves performance
%
%   PHASE 2 - Interactions (hierarchy-respecting):
%     5. Only test interactions where BOTH parent main effects were selected
%        e.g., Speed×TF is only eligible if Speed AND TF were added in Phase 1
%     6. Add the best eligible interaction if Δ CV bps > threshold
%     7. Repeat until no eligible interaction improves performance
%
%   This ensures proper model hierarchy and interpretability.
%
% COMPARISON MODELS (always fitted for visualization):
%   Null:            Intercept + onset kernel only
%   Selected:        Forward-selection winner (neuron-specific)
%   Additive:        All main effects (Speed + TF + SF + OR)
%   FullInteraction: Additive + all pairwise interactions
%
% NEURON CLASSIFICATION (based on forward selection):
%   - is_spprofileeed_tuned: Speed in selected variables
%   - is_tf_tuned:    TF in selected variables
%   - is_sf_tuned:    SF in selected variables
%   - is_or_tuned:    OR in selected variables
%   - has_interaction: Any interaction in selected variables
%
% OUTPUTS:
%   - prefilter_decision_tree.csv: Pre-filtering + GLM selection results
%   - csv/glm_model_comparison.csv: Model comparison with selected variables
%   - csv/glm_coefficients.csv: GLM coefficients for each cluster
%   - csv/glm_selection_history.csv: Forward selection history per cluster
%
% See also: aritmetic_sum, speed_tuning_all_clusters_all_conditions,
%           tf_tuning_all_clusters_all_conditions, create_tables_by_batch

%% ====================================================================
%  Section 1: Configuration
%  ====================================================================
close all; clearvars;

% --- Initialize parallel pool ---
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
    fprintf('Started parallel pool with %d workers\n', pool.NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', pool.NumWorkers);
end

% --- Initialize RC2Analysis controller for path configuration ---
ctl = RC2Analysis();

% --- Open diary to capture all printed output ---
log_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
log_filename = sprintf('glm_analysis_log_%s.txt', log_timestamp);
log_file_tmp = fullfile(tempdir, log_filename);   % temp location until fig dir is known
diary(log_file_tmp);
fprintf('Analysis started: %s\n\n', datestr(now));

experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};
restricted              = false;   % only clusters tuned in V and T_Vstatic

% Binning (matches pre-computed tuning tables: 20 bins at 5% percentile each)
n_speed_bins            = 20;     % speed bins per trial (for observed data plotting)
n_tf_bins               = 20;     % TF bins per trial

% Time-bin GLM settings
time_bin_width          = 0.1;  % 2 ms bins

% Basis functions: raised cosine
n_speed_bases           = 5;      % number of raised cosine bases for speed
n_tf_bases              = 5;      % number of raised cosine bases for TF
speed_range             = [0, 50];  % cm/s (updated from tuning tables later)
tf_range                = [0, 7.3]; % Hz   (updated from tuning tables later)

% Onset dynamics kernel (Park et al. 2014 style)
% Captures transient response at motion onset / offset
n_onset_bases           = 6;      % number of raised cosine bases for onset kernel
onset_range             = [0, 2.0];  % seconds after motion onset (0 = onset time)
% Note: stationary periods have time_since_onset < 0, so onset bases = 0

% Cross-validation
n_cv_folds              = 5;

% Prediction cap (physiological limit for firing rate)
max_predicted_fr        = 400;    % Hz - cap predictions to prevent exp() explosion

% Observed FR smoothing for visualization (Park et al. 2014 style)
% The model operates on raw spike counts, but for visualization we smooth
% the observed spike train to match the naturally smooth model predictions.
obs_fr_smooth_width     = 0.100;  % 100 ms boxcar (Park 2014 used 100ms)

% Figures
save_figs               = true;
overwrite               = true;
figure_dir              = {'glm_single_cluster'};

% CSV output
csv_output              = true;

% --- Pre-computed tuning table paths (via path_config) ---
tuning_curves_dir       = fullfile(ctl.path_config.formatted_data_dir, 'csvs', 'tuning_curves');
tf_tuning_curves_dir    = fullfile(ctl.path_config.formatted_data_dir, 'csvs', 'tf_tuning_curves');

% Motion cloud data paths (via path_config)
mc_sequence_path        = fullfile(ctl.path_config.motion_clouds_root, 'motion_cloud_sequence_250414.mat');
mc_folders_path         = fullfile(ctl.path_config.motion_clouds_root, 'image_folders.mat');

% --- Batch definitions ---
exclude_patterns = {'theta0p000_Btheta3p142_sf00p006_Bsf0p004_VX0p000_BV2p000'};

sf_values_map = containers.Map(...
    {'sf00p003', 'sf00p006', 'sf00p012'}, ...
    {0.003, 0.006, 0.012});
sf_keys = keys(sf_values_map);

or_values_map = containers.Map(...
    {'theta0p000', 'theta0p785', 'theta1p571', 'theta-0p785'}, ...
    {0, pi/4, pi/2, -pi/4});
or_keys = keys(or_values_map);

batch_patterns = { ...
    {'sf00p003_Bsf0p002_VX1p002', 'sf00p006_Bsf0p002_VX0p501', 'sf00p012_Bsf0p002_VX0p250'}, ...
    {'sf00p003_Bsf0p002_VX2p003', 'sf00p006_Bsf0p002_VX1p002', 'sf00p012_Bsf0p002_VX0p501'}, ...
    {'sf00p003_Bsf0p002_VX4p006', 'sf00p006_Bsf0p002_VX2p003', 'sf00p012_Bsf0p002_VX1p002'} };
batch_gains = [1/30, 2/30, 4/30];  % Hz/(cm/s)

% Hierarchical Forward Selection Framework
% ========================================
% Two-phase selection to ensure proper model hierarchy:
%
% PHASE 1 - Main effects:
%   Test: Speed, TF, SF, OR (one at a time)
%   Add best if Δ CV bps > threshold, repeat until none pass
%
% PHASE 2 - Interactions (hierarchy-respecting):
%   Only test interactions where BOTH parents were selected in Phase 1
%   e.g., Speed_x_TF is only eligible if Speed AND TF were added
%   Add best if Δ CV bps > threshold, repeat until none pass

% Candidate variables for forward selection
candidate_main_effects = {'Speed', 'TF', 'SF', 'OR'};
candidate_interactions = {'Speed_x_TF', 'Speed_x_SF', 'Speed_x_OR', ...
                          'TF_x_SF', 'TF_x_OR', 'SF_x_OR'};
all_candidates = [candidate_main_effects, candidate_interactions];
n_candidates = length(all_candidates);

% Selection threshold (delta bits-per-spike for inclusion)
% NOTE: 0.01 was Hardcastle 2017 default, but may be too conservative.
% A proper approach would use a null distribution from shuffled data.
selection_threshold = 0.005;  % Δ bps > 0.005 = variable contributes meaningfully

% Comparison models for visualization (always fitted for reference)
comparison_models = {'Null', 'Additive', 'FullInteraction'};
n_comparison_models = length(comparison_models);

% Keep backward compatibility: old variable name
classification_threshold = selection_threshold;

% Model formula strings (for figure titles)
model_formulas = containers.Map();
model_formulas('Null')            = '\beta_0 + K_{onset}';
model_formulas('M0')              = '\beta_0';
model_formulas('M0_Speed')        = '\beta_0 + f(Spd) + K_{onset}';
model_formulas('M0_Speed_TF')     = '\beta_0 + f(Spd) + g(TF) + K_{onset}';
model_formulas('M0_Speed_TF_SF')  = '\beta_0 + f(Spd) + g(TF) + SF + K_{onset}';
model_formulas('Additive')        = '\beta_0 + f(Spd) + g(TF) + SF + OR + K_{onset}';
model_formulas('FullInteraction') = 'Additive + Spd*TF + Spd*SF + Spd*OR + TF*SF + TF*OR + SF*OR';
model_formulas('Selected')        = '\beta_0 + K_{onset} + [selected vars]';
model_formulas('Additive_no_Speed') = '\beta_0 + g(TF) + SF + OR + K_{onset}';
model_formulas('Additive_no_TF')    = '\beta_0 + f(Spd) + SF + OR + K_{onset}';
model_formulas('Additive_no_SF')    = '\beta_0 + f(Spd) + g(TF) + OR + K_{onset}';
model_formulas('Additive_no_OR')    = '\beta_0 + f(Spd) + g(TF) + SF + K_{onset}';

fprintf('=== GLM Single Cluster Analysis (v10: Hierarchical Forward Selection) ===\n');
fprintf('Configuration:\n');
fprintf('  Time-bin width: %.0f ms\n', time_bin_width*1000);
fprintf('  Speed bases: %d | TF bases: %d | Onset bases: %d\n', n_speed_bases, n_tf_bases, n_onset_bases);
fprintf('  Onset kernel range: [%.1f, %.1f] s\n', onset_range(1), onset_range(2));
fprintf('  CV folds: %d | Restricted: %d\n', n_cv_folds, restricted);
fprintf('  Speed range: [%.1f, %.1f] cm/s\n', speed_range(1), speed_range(2));
fprintf('  TF range: [%.1f, %.1f] Hz\n', tf_range(1), tf_range(2));
fprintf('  Onset dynamics: YES (K_onset captures transient after motion onset)\n');
fprintf('  Stationary baseline: YES (beta_0 anchored to stationary periods where onset=0)\n');
fprintf('  Observed FR smoothing: %.0f ms boxcar (for visualization only)\n', obs_fr_smooth_width*1000);
fprintf('  Model selection: HIERARCHICAL FORWARD SELECTION\n');
fprintf('    Phase 1: Main effects (%s)\n', strjoin(candidate_main_effects, ', '));
fprintf('    Phase 2: Interactions (only if both parents selected)\n');
fprintf('  Selection threshold: delta bps > %.3f\n', selection_threshold);
fprintf('  Comparison models: %s\n', strjoin(comparison_models, ', '));
fprintf('  Pre-filtering: Decision tree using SVM and tuning tables\n\n');

%% ====================================================================
%  Section 2: Load motion cloud sequence and cloud names
%  ====================================================================
fprintf('--- Loading motion cloud metadata ---\n');

P = load(mc_sequence_path);
if isfield(P, 'presentation_sequence')
    mc_sequence = P.presentation_sequence;
else
    fns = fieldnames(P);
    for i = 1:numel(fns)
        v = P.(fns{i});
        if isnumeric(v) && isvector(v)
            mc_sequence = v; break;
        end
    end
end

S = load(mc_folders_path);
fns = fieldnames(S);
for i = 1:numel(fns)
    v = S.(fns{i});
    if isstruct(v) && isfield(v, 'name')
        cloud_names = {v.name};
        break;
    elseif iscell(v)
        cloud_names = v(:)';
        break;
    elseif isstring(v)
        cloud_names = cellstr(v(:))';
        break;
    end
end

fprintf('  Loaded %d sequence entries, %d cloud names\n', length(mc_sequence), length(cloud_names));

%% ====================================================================
%  Section 2b: Pre-filtering Decision Tree
%  ====================================================================
% This section implements the pre-filtering decision tree to reduce the
% number of clusters that need full GLM analysis. The logic is:
%
% For each cluster, test SVM significance (stationary vs motion) for T, V, VT:
%
% SPEED TUNED CLASSIFICATION:
%   T significant OR VT significant -> cluster is speed tuned
%
% GLM SELECTION:
%   VT significant -> run GLM
%   OR any combo of 2 significant conditions:
%     - T + V significant
%     - T + VT significant
%     - V + VT significant
%   -> run GLM
%
%   Otherwise (0 or 1 significant, and that 1 is not VT) -> do not run GLM
%
% The pre-filtering uses:
%   - svm_table: stationary vs motion firing rates per trial
%   - Wilcoxon sign-rank test for stationary vs motion comparisons

fprintf('\n--- Pre-filtering Decision Tree ---\n');

probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

% --- Move diary log into the figures directory ---
diary off;
log_file = fullfile(ctl.figs.curr_dir, log_filename);
movefile(log_file_tmp, log_file);
diary(log_file);
fprintf('Log file: %s\n', log_file);

% Initialize pre-filtering results structure
prefilter_results = struct();
prefilter_results.probe_id = {};
prefilter_results.cluster_id = [];
prefilter_results.is_T_significant = [];       % stationary_T != motion_T (SVM)
prefilter_results.is_V_significant = [];       % stationary_V != motion_V (SVM)
prefilter_results.is_VT_significant = [];      % stationary_VT != motion_VT (SVM)
prefilter_results.is_spprofileeed_tuned = [];         % T OR VT significant
prefilter_results.should_run_glm = [];         % final decision: run GLM or not
prefilter_results.category = {};               % category label for this cluster
prefilter_results.p_T_svm = [];                % p-value for T stationary vs motion
prefilter_results.p_V_svm = [];                % p-value for V stationary vs motion
prefilter_results.p_VT_svm = [];               % p-value for VT stationary vs motion
prefilter_results.direction_T = [];            % direction of T modulation
prefilter_results.direction_V = [];            % direction of V modulation
prefilter_results.direction_VT = [];           % direction of VT modulation
prefilter_results.n_significant = [];          % number of significant conditions (0, 1, 2, or 3)

% Significance threshold for Wilcoxon sign-rank tests
p_threshold = 0.05;

% Counters for summary report
n_prefilter_clusters = 0;
n_T_significant = 0;
n_V_significant = 0;
n_VT_significant = 0;
n_speed_tuned = 0;         % T OR VT significant
n_should_run_glm = 0;
n_glm_by_VT_only = 0;      % VT significant alone qualifies for GLM
n_glm_by_T_and_V = 0;      % T+V both significant (no VT)
n_glm_by_T_and_VT = 0;     % T+VT both significant
n_glm_by_V_and_VT = 0;     % V+VT both significant
n_glm_by_all_three = 0;    % all 3 significant
n_no_significant = 0;      % none significant
n_one_significant_not_VT = 0;  % only T or only V significant (excluded)

for probe_i = 1:length(probe_ids)
    pid = probe_ids{probe_i};
    fprintf('\nPre-filtering Probe %d/%d: %s\n', probe_i, length(probe_ids), pid);
    
    % Load formatted data for this probe
    data = ctl.load_formatted_data(pid);
    visp_clusters = data.VISp_clusters();
    
    if isempty(visp_clusters)
        fprintf('  No VISp clusters, skipping\n');
        continue;
    end
    
    fprintf('  VISp clusters: %d\n', length(visp_clusters));
    
    % Process each VISp cluster
    for ci = 1:length(visp_clusters)
        cid = visp_clusters(ci).id;
        n_prefilter_clusters = n_prefilter_clusters + 1;
        
        prefilter_results.probe_id{n_prefilter_clusters} = pid;
        prefilter_results.cluster_id(n_prefilter_clusters) = cid;
        
        % ====== Step 1: Check T significance (stationary_T vs motion_T) ======
        [sig_T, p_T, dir_T] = data.is_stationary_vs_motion_significant(cid, 'T_Vstatic');
        prefilter_results.p_T_svm(n_prefilter_clusters) = p_T;
        prefilter_results.direction_T(n_prefilter_clusters) = dir_T;
        is_T_significant = sig_T && ~isempty(sig_T);
        prefilter_results.is_T_significant(n_prefilter_clusters) = is_T_significant;
        
        % ====== Step 2: Check V significance (stationary_V vs motion_V) ======
        [sig_V, p_V, dir_V] = data.is_stationary_vs_motion_significant(cid, 'V');
        prefilter_results.p_V_svm(n_prefilter_clusters) = p_V;
        prefilter_results.direction_V(n_prefilter_clusters) = dir_V;
        is_V_significant = sig_V && ~isempty(sig_V);
        prefilter_results.is_V_significant(n_prefilter_clusters) = is_V_significant;
        
        % ====== Step 3: Check VT significance (stationary_VT vs motion_VT) ======
        [sig_VT, p_VT, dir_VT] = data.is_stationary_vs_motion_significant(cid, 'VT');
        prefilter_results.p_VT_svm(n_prefilter_clusters) = p_VT;
        prefilter_results.direction_VT(n_prefilter_clusters) = dir_VT;
        is_VT_significant = sig_VT && ~isempty(sig_VT);
        prefilter_results.is_VT_significant(n_prefilter_clusters) = is_VT_significant;
        
        % ====== Count significant conditions ======
        n_sig = is_T_significant + is_V_significant + is_VT_significant;
        prefilter_results.n_significant(n_prefilter_clusters) = n_sig;
        
        % ====== Update counters for individual conditions ======
        if is_T_significant
            n_T_significant = n_T_significant + 1;
        end
        if is_V_significant
            n_V_significant = n_V_significant + 1;
        end
        if is_VT_significant
            n_VT_significant = n_VT_significant + 1;
        end
        
        % ====== SPEED TUNED CLASSIFICATION: T OR VT significant ======
        is_spprofileeed_tuned = is_T_significant || is_VT_significant;
        prefilter_results.is_spprofileeed_tuned(n_prefilter_clusters) = is_spprofileeed_tuned;
        if is_spprofileeed_tuned
            n_speed_tuned = n_speed_tuned + 1;
        end
        
        % ====== GLM SELECTION: VT significant OR any combo of 2 significant ======
        % VT significant -> run GLM
        % OR T+V both significant -> run GLM  
        % OR T+VT both significant -> run GLM
        % OR V+VT both significant -> run GLM
        should_run_glm = false;
        category = 'not_selected';
        
        if is_VT_significant && is_T_significant && is_V_significant
            % All 3 significant
            should_run_glm = true;
            category = 'all_three_significant';
            n_glm_by_all_three = n_glm_by_all_three + 1;
        elseif is_VT_significant && is_T_significant
            % T + VT significant
            should_run_glm = true;
            category = 'T_and_VT_significant';
            n_glm_by_T_and_VT = n_glm_by_T_and_VT + 1;
        elseif is_VT_significant && is_V_significant
            % V + VT significant
            should_run_glm = true;
            category = 'V_and_VT_significant';
            n_glm_by_V_and_VT = n_glm_by_V_and_VT + 1;
        elseif is_VT_significant
            % VT significant alone
            should_run_glm = true;
            category = 'VT_only_significant';
            n_glm_by_VT_only = n_glm_by_VT_only + 1;
        elseif is_T_significant && is_V_significant
            % T + V significant (but not VT)
            should_run_glm = true;
            category = 'T_and_V_significant';
            n_glm_by_T_and_V = n_glm_by_T_and_V + 1;
        elseif n_sig == 0
            % None significant
            category = 'none_significant';
            n_no_significant = n_no_significant + 1;
        else
            % Only T or only V significant (excluded from GLM)
            if is_T_significant
                category = 'T_only_significant';
            else
                category = 'V_only_significant';
            end
            n_one_significant_not_VT = n_one_significant_not_VT + 1;
        end
        
        if should_run_glm
            n_should_run_glm = n_should_run_glm + 1;
        end
        
        prefilter_results.should_run_glm(n_prefilter_clusters) = should_run_glm;
        prefilter_results.category{n_prefilter_clusters} = category;
    end
end

% Print pre-filtering summary
fprintf('\n====================================================================\n');
fprintf('  PRE-FILTERING SUMMARY\n');
fprintf('====================================================================\n');
fprintf('  Total VISp clusters scanned: %d\n', n_prefilter_clusters);
fprintf('\n  --- SVM Significance (stationary vs motion) ---\n');
fprintf('    T significant (SVM):    %d (%.1f%%)\n', n_T_significant, 100*n_T_significant/n_prefilter_clusters);
fprintf('    V significant (SVM):    %d (%.1f%%)\n', n_V_significant, 100*n_V_significant/n_prefilter_clusters);
fprintf('    VT significant (SVM):   %d (%.1f%%)\n', n_VT_significant, 100*n_VT_significant/n_prefilter_clusters);
fprintf('\n  --- Speed Tuned Classification (T OR VT significant) ---\n');
fprintf('    Speed tuned:            %d (%.1f%%)\n', n_speed_tuned, 100*n_speed_tuned/n_prefilter_clusters);
fprintf('    Not speed tuned:        %d (%.1f%%)\n', n_prefilter_clusters - n_speed_tuned, 100*(n_prefilter_clusters - n_speed_tuned)/n_prefilter_clusters);
fprintf('\n  --- GLM Selection (VT sig, or any 2 of T/V/VT sig) ---\n');
fprintf('    All three significant (T+V+VT):   %d (%.1f%%)\n', n_glm_by_all_three, 100*n_glm_by_all_three/n_prefilter_clusters);
fprintf('    T+VT significant (no V):          %d (%.1f%%)\n', n_glm_by_T_and_VT, 100*n_glm_by_T_and_VT/n_prefilter_clusters);
fprintf('    V+VT significant (no T):          %d (%.1f%%)\n', n_glm_by_V_and_VT, 100*n_glm_by_V_and_VT/n_prefilter_clusters);
fprintf('    VT only significant:              %d (%.1f%%)\n', n_glm_by_VT_only, 100*n_glm_by_VT_only/n_prefilter_clusters);
fprintf('    T+V significant (no VT):          %d (%.1f%%)\n', n_glm_by_T_and_V, 100*n_glm_by_T_and_V/n_prefilter_clusters);
fprintf('    -------------------------------------------\n');
fprintf('    TOTAL for GLM:                    %d (%.1f%%)\n', n_should_run_glm, 100*n_should_run_glm/n_prefilter_clusters);
fprintf('\n  --- Excluded from GLM ---\n');
fprintf('    None significant:                 %d (%.1f%%)\n', n_no_significant, 100*n_no_significant/n_prefilter_clusters);
fprintf('    Only T or V significant (not VT): %d (%.1f%%)\n', n_one_significant_not_VT, 100*n_one_significant_not_VT/n_prefilter_clusters);
fprintf('    -------------------------------------------\n');
fprintf('    TOTAL excluded:                   %d (%.1f%%)\n', n_prefilter_clusters - n_should_run_glm, 100*(n_prefilter_clusters - n_should_run_glm)/n_prefilter_clusters);
fprintf('====================================================================\n\n');
% Save pre-filtering results to CSV
prefilter_table = table(...
    string(prefilter_results.probe_id'), ...
    prefilter_results.cluster_id', ...
    prefilter_results.is_T_significant', ...
    prefilter_results.is_V_significant', ...
    prefilter_results.is_VT_significant', ...
    prefilter_results.n_significant', ...
    prefilter_results.is_spprofileeed_tuned', ...
    prefilter_results.should_run_glm', ...
    string(prefilter_results.category'), ...
    prefilter_results.p_T_svm', ...
    prefilter_results.p_V_svm', ...
    prefilter_results.p_VT_svm', ...
    prefilter_results.direction_T', ...
    prefilter_results.direction_V', ...
    prefilter_results.direction_VT', ...
    'VariableNames', {'probe_id', 'cluster_id', 'is_T_significant', 'is_V_significant', ...
        'is_VT_significant', 'n_significant', 'is_spprofileeed_tuned', 'should_run_glm', 'category', ...
        'p_T_svm', 'p_V_svm', 'p_VT_svm', ...
        'direction_T', 'direction_V', 'direction_VT'});

prefilter_csv_path = fullfile(ctl.figs.curr_dir, 'prefilter_decision_tree.csv');
writetable(prefilter_table, prefilter_csv_path);
fprintf('Pre-filtering results saved to: %s\n\n', prefilter_csv_path);

% Create a set of cluster IDs that should be analyzed
clusters_for_glm = containers.Map('KeyType', 'char', 'ValueType', 'any');
for pfi = 1:n_prefilter_clusters
    if prefilter_results.should_run_glm(pfi)
        key = sprintf('%s_%d', prefilter_results.probe_id{pfi}, prefilter_results.cluster_id(pfi));
        clusters_for_glm(key) = struct(...
            'probe_id', prefilter_results.probe_id{pfi}, ...
            'cluster_id', prefilter_results.cluster_id(pfi), ...
            'category', prefilter_results.category{pfi});
    end
end

fprintf('Clusters passing pre-filter for GLM: %d\n\n', clusters_for_glm.Count);

%% ====================================================================
%  Section 3: Data collection — Speed-bin data for observed tuning curves
%  ====================================================================
% Note: Speed-bin data is still used for plotting observed data (Row 1 in
% tuning curve figures). GLM fitting is done only with time-bin data.
fprintf('\n--- Loading pre-computed tuning tables (for observed data plotting) ---\n');

% Helper: check whether a string contains any substring in a cell array
contains_any = @(str, subs) any(cellfun(@(s) contains(str, s), subs));

% Build stimulus lookup: trial_id -> {sf, batch_gain, or, is_excluded}
trial_stim = struct();
for tid = 1:length(mc_sequence)
    mc_id = mc_sequence(tid);
    if mc_id < 1 || mc_id > length(cloud_names)
        trial_stim(tid).sf = NaN;
        trial_stim(tid).batch_gain = NaN;
        trial_stim(tid).or = NaN;
        trial_stim(tid).excluded = true;
        continue;
    end
    cname = cloud_names{mc_id};
    
    if contains_any(cname, exclude_patterns)
        trial_stim(tid).sf = NaN;
        trial_stim(tid).batch_gain = NaN;
        trial_stim(tid).or = NaN;
        trial_stim(tid).excluded = true;
    else
        [sf_v, ~, or_v] = parse_cloud_name(cname, sf_keys, sf_values_map, or_keys, or_values_map);
        trial_stim(tid).sf = sf_v;
        trial_stim(tid).or = or_v;
        trial_stim(tid).excluded = false;
        
        trial_stim(tid).batch_gain = NaN;
        for bi = 1:length(batch_patterns)
            if contains_any(cname, batch_patterns{bi})
                trial_stim(tid).batch_gain = batch_gains(bi);
                break;
            end
        end
    end
end
fprintf('  Built stimulus lookup for %d trials\n', length(trial_stim));

% --- Load pre-computed bin edges ---
first_spd_file = fullfile(tuning_curves_dir, [probe_ids{1} '.mat']);
first_tf_file  = fullfile(tf_tuning_curves_dir, [probe_ids{1} '.mat']);
D_tmp = load(first_spd_file);
spd_bin_edges   = D_tmp.tuning_curves{1}(1).bin_edges(:)';
spd_bin_centers = D_tmp.tuning_curves{1}(1).bin_centers(:)';
D_tmp = load(first_tf_file);
tf_bin_edges   = D_tmp.tuning_curves{1}(1).bin_edges(:)';
tf_bin_centers = D_tmp.tuning_curves{1}(1).bin_centers(:)';
clear D_tmp;

speed_range = [spd_bin_edges(1), spd_bin_edges(end)];
tf_range    = [tf_bin_edges(1),  tf_bin_edges(end)];
fprintf('  Speed range from tuning table: [%.2f, %.2f] cm/s\n', speed_range(1), speed_range(2));
fprintf('  TF range from tuning table: [%.4f, %.4f] Hz\n', tf_range(1), tf_range(2));

% --- Build speed-bin master table ---
% Each row = one (cluster x trial x speed_bin), spike_count = rate * time
max_rows = 500000;
col_probe_id    = cell(max_rows, 1);
col_cluster_id  = zeros(max_rows, 1);
col_trial_id    = zeros(max_rows, 1);
col_condition   = cell(max_rows, 1);
col_speed       = zeros(max_rows, 1);
col_tf          = zeros(max_rows, 1);
col_sf          = zeros(max_rows, 1);
col_orientation = zeros(max_rows, 1);
col_batch_gain  = zeros(max_rows, 1);
col_spike_count = zeros(max_rows, 1);
col_time_in_bin = zeros(max_rows, 1);
col_bin_idx     = zeros(max_rows, 1);

total_clusters = 0;
total_rows = 0;

% --- Also collect info needed for time-bin GLM ---
% We store: per probe, which cluster IDs are VISp (and restricted),
% and which trial_ids appear in the tuning tables per condition.
probe_info = struct();

for probe_i = 1:length(probe_ids)
    pid = probe_ids{probe_i};
    fprintf('\nProbe %d/%d: %s\n', probe_i, length(probe_ids), pid);
    
    % --- Load FormattedData for VISp filtering ---
    data = ctl.load_formatted_data(pid);
    visp_clusters = data.VISp_clusters();
    visp_ids = arrayfun(@(x) x.id, visp_clusters);
    
    % --- Apply pre-filtering: only keep clusters that passed the decision tree ---
    keep = false(size(visp_ids));
    for ci = 1:length(visp_ids)
        key = sprintf('%s_%d', pid, visp_ids(ci));
        if clusters_for_glm.isKey(key)
            keep(ci) = true;
        end
    end
    visp_ids_filtered = visp_ids(keep);
    
    fprintf('  VISp clusters: %d | After pre-filter: %d\n', length(visp_ids), length(visp_ids_filtered));
    
    % Apply additional restricted filtering if enabled
    if restricted && ~isempty(visp_ids_filtered)
        keep_restr = true(size(visp_ids_filtered));
        for ci = 1:length(visp_ids_filtered)
            [~, ~, dir_T] = data.is_stationary_vs_motion_significant(visp_ids_filtered(ci), 'T_Vstatic');
            [~, ~, dir_V] = data.is_stationary_vs_motion_significant(visp_ids_filtered(ci), 'V');
            if dir_T == 0 || dir_V == 0
                keep_restr(ci) = false;
            end
        end
        visp_ids_filtered = visp_ids_filtered(keep_restr);
        fprintf('  After restricted filter: %d\n', length(visp_ids_filtered));
    end
    
    if isempty(visp_ids_filtered)
        fprintf('  No clusters remaining, skipping\n');
        continue;
    end
    total_clusters = total_clusters + length(visp_ids_filtered);
    
    visp_id_set = visp_ids_filtered(:)';
    
    % Store probe info for time-bin GLM later
    probe_info(probe_i).pid = pid;
    probe_info(probe_i).data = data;
    probe_info(probe_i).visp_ids = visp_id_set;
    
    % --- Load tuning tables (speed and TF) ---
    tuning_file = fullfile(tuning_curves_dir, [pid '.mat']);
    if ~exist(tuning_file, 'file')
        fprintf('  WARNING: tuning file not found: %s, skipping\n', tuning_file);
        continue;
    end
    D_spd = load(tuning_file);
    probe_info(probe_i).D_spd = D_spd;  % Store for Row 1 plotting
    fprintf('  Loaded speed tuning: %s\n', tuning_file);
    
    % Load TF tuning table
    tf_tuning_file = fullfile(tf_tuning_curves_dir, [pid '.mat']);
    if exist(tf_tuning_file, 'file')
        D_tf = load(tf_tuning_file);
        probe_info(probe_i).D_tf = D_tf;  % Store for Row 1 plotting
        fprintf('  Loaded TF tuning: %s\n', tf_tuning_file);
    else
        probe_info(probe_i).D_tf = [];
        fprintf('  WARNING: TF tuning file not found: %s\n', tf_tuning_file);
    end
    
    % --- Process each condition (simplified: rate * time from tables) ---
    for cond_i = 1:length(trial_group_labels)
        condition = trial_group_labels{cond_i};
        
        cond_match = find(strcmp(D_spd.trial_groups, condition));
        if isempty(cond_match)
            fprintf('    Condition %s not found, skipping\n', condition);
            continue;
        end
        
        tc_array = D_spd.tuning_curves{cond_match};
        if isempty(tc_array)
            fprintf('    Condition %s: empty tuning curves, skipping\n', condition);
            continue;
        end
        
        tc_cluster_ids = arrayfun(@(x) double(x.cluster_id), tc_array);
        
        n_cond_rows = 0;
        for ci = 1:length(visp_id_set)
            cid = visp_id_set(ci);
            tc_idx = find(tc_cluster_ids == cid, 1);
            if isempty(tc_idx)
                continue;
            end
            
            tc = tc_array(tc_idx);
            bin_centers = tc.bin_centers(:)';
            n_bins = length(bin_centers);
            trial_ids = tc.trial_ids(:)';
            tuning_matrix = tc.tuning;    % n_bins x n_trials (firing rate, spk/s)
            time_matrix   = tc.timing;    % n_bins x n_trials (seconds)
            
            for ti = 1:length(trial_ids)
                trial_id = trial_ids(ti);
                
                if trial_id < 1 || trial_id > length(trial_stim)
                    continue;
                end
                stim = trial_stim(trial_id);
                if stim.excluded
                    continue;
                end
                
                sf_val = stim.sf;
                gain_val = stim.batch_gain;
                or_val = stim.or;
                
                for bi = 1:n_bins
                    t_in_bin = time_matrix(bi, ti);
                    fr_in_bin = tuning_matrix(bi, ti);
                    
                    if isnan(fr_in_bin) || isnan(t_in_bin) || t_in_bin < 0.01
                        continue;
                    end
                    
                    speed_val = bin_centers(bi);
                    
                    % Estimated spike count = rate * time
                    n_spikes = fr_in_bin * t_in_bin;
                    
                    % Compute TF and features based on condition
                    switch condition
                        case 'T_Vstatic'
                            tf_val = 0;
                            sf_trial = NaN;
                            or_trial = NaN;
                            gain_trial = NaN;
                        case 'V'
                            tf_val = gain_val * speed_val;
                            sf_trial = sf_val;
                            or_trial = or_val;
                            gain_trial = gain_val;
                            speed_val = 0;  % no locomotion in V
                        case 'VT'
                            tf_val = gain_val * speed_val;
                            sf_trial = sf_val;
                            or_trial = or_val;
                            gain_trial = gain_val;
                        otherwise
                            tf_val = 0;
                            sf_trial = NaN;
                            or_trial = NaN;
                            gain_trial = NaN;
                    end
                    
                    total_rows = total_rows + 1;
                    n_cond_rows = n_cond_rows + 1;
                    
                    if total_rows > max_rows
                        max_rows = max_rows * 2;
                        col_probe_id{max_rows}    = [];
                        col_cluster_id(max_rows)  = 0;
                        col_trial_id(max_rows)    = 0;
                        col_condition{max_rows}   = [];
                        col_speed(max_rows)       = 0;
                        col_tf(max_rows)          = 0;
                        col_sf(max_rows)          = 0;
                        col_orientation(max_rows) = 0;
                        col_batch_gain(max_rows)  = 0;
                        col_spike_count(max_rows) = 0;
                        col_time_in_bin(max_rows) = 0;
                        col_bin_idx(max_rows)     = 0;
                    end
                    
                    col_probe_id{total_rows}    = pid;
                    col_cluster_id(total_rows)  = cid;
                    col_trial_id(total_rows)    = trial_id;
                    col_condition{total_rows}   = condition;
                    col_speed(total_rows)       = speed_val;
                    col_tf(total_rows)          = tf_val;
                    col_sf(total_rows)          = sf_trial;
                    col_orientation(total_rows) = or_trial;
                    col_batch_gain(total_rows)  = gain_trial;
                    col_spike_count(total_rows) = n_spikes;
                    col_time_in_bin(total_rows) = t_in_bin;
                    col_bin_idx(total_rows)     = bi;
                end
            end
        end
        fprintf('    %s: %d rows\n', condition, n_cond_rows);
    end
end

fprintf('\n--- Speed-bin data collection complete (for observed data plotting) ---\n');
fprintf('  Total clusters: %d\n', total_clusters);
fprintf('  Total rows (bins): %d\n', total_rows);

% Trim and build table
T_master_spd = table(...
    string(col_probe_id(1:total_rows)), ...
    col_cluster_id(1:total_rows), ...
    col_trial_id(1:total_rows), ...
    string(col_condition(1:total_rows)), ...
    col_speed(1:total_rows), ...
    col_tf(1:total_rows), ...
    col_sf(1:total_rows), ...
    col_orientation(1:total_rows), ...
    col_batch_gain(1:total_rows), ...
    col_spike_count(1:total_rows), ...
    col_time_in_bin(1:total_rows), ...
    col_bin_idx(1:total_rows), ...
    'VariableNames', {'probe_id', 'cluster_id', 'trial_id', 'condition', ...
        'speed', 'tf', 'sf', 'orientation', 'batch_gain', ...
        'spike_count', 'time_in_bin', 'bin_idx'});

fprintf('  Speed-bin table: %d rows x %d columns (used for observed data plotting)\n', height(T_master_spd), width(T_master_spd));

%% ====================================================================
%  Section 3b: Data collection — Time-bin GLM (100 ms bins from raw data)
%  ====================================================================
fprintf('\n--- Building time-bin GLM data (%.0f ms bins) [PARALLEL] ---\n', time_bin_width*1000);

% Collect all tables from each probe, then concatenate
all_probe_tables = cell(length(probe_info), 1);
max_trial_duration = 0;  % track for temporal basis range

for probe_i = 1:length(probe_info)
    if isempty(probe_info(probe_i).pid)
        continue;
    end
    pid = probe_info(probe_i).pid;
    data = probe_info(probe_i).data;
    visp_id_set = probe_info(probe_i).visp_ids;
    
    fprintf('\nProbe %d: %s (time-bin)\n', probe_i, pid);
    
    % Collect all trial IDs for this probe across conditions
    all_trial_ids_map = containers.Map('KeyType', 'int32', 'ValueType', 'char');
    for cond_i = 1:length(trial_group_labels)
        condition = trial_group_labels{cond_i};
        trials_cond = data.get_trials_with_trial_group_label(condition);
        for ti_k = 1:length(trials_cond)
            tid_k = trials_cond{ti_k}.trial_id;
            if tid_k >= 1 && tid_k <= length(trial_stim) && ~trial_stim(tid_k).excluded
                all_trial_ids_map(int32(tid_k)) = condition;
            end
        end
    end
    
    all_trial_ids = cell2mat(keys(all_trial_ids_map));
    all_trial_conditions = values(all_trial_ids_map, num2cell(all_trial_ids));
    n_trials_probe = length(all_trial_ids);
    fprintf('  %d trials to process\n', n_trials_probe);
    
    % Pre-load spike times for VISp clusters into cell array (parfor-compatible)
    n_visp = length(visp_id_set);
    spike_times_all = cell(n_visp, 1);
    valid_cids_mask = false(n_visp, 1);
    for ci_cache = 1:n_visp
        cid_cache = visp_id_set(ci_cache);
        try
            cluster_obj = data.get_cluster_with_id(cid_cache);
            spike_times_all{ci_cache} = cluster_obj.spike_times;
            valid_cids_mask(ci_cache) = true;
        catch ME
            fprintf('    WARNING: cluster %d spike_times: %s\n', cid_cache, ME.message);
        end
    end
    valid_cids = visp_id_set(valid_cids_mask);
    spike_times_valid = spike_times_all(valid_cids_mask);
    n_valid_clusters = length(valid_cids);
    fprintf('  Spike times cached for %d clusters\n', n_valid_clusters);
    
    % Pre-load all trial objects (cannot load in parfor)
    fprintf('  Pre-loading trial objects...\n');
    trial_objs = cell(n_trials_probe, 1);
    for ti_pre = 1:n_trials_probe
        trial_id = all_trial_ids(ti_pre);
        try
            trial_obj = data.get_trials_with_trial_ids(trial_id);
            if iscell(trial_obj)
                trial_obj = trial_obj{1};
            end
            trial_objs{ti_pre} = trial_obj.to_aligned;
        catch
            trial_objs{ti_pre} = [];
        end
    end
    
    % Prepare stimulus info arrays for parfor
    trial_sf_vals = zeros(n_trials_probe, 1);
    trial_gain_vals = zeros(n_trials_probe, 1);
    trial_or_vals = zeros(n_trials_probe, 1);
    for ti_pre = 1:n_trials_probe
        stim = trial_stim(all_trial_ids(ti_pre));
        trial_sf_vals(ti_pre) = stim.sf;
        trial_gain_vals(ti_pre) = stim.batch_gain;
        trial_or_vals(ti_pre) = stim.or;
    end
    
    % === PARALLEL LOOP over trials ===
    trial_results = cell(n_trials_probe, 1);  % Each cell will hold rows for this trial
    trial_max_dur = zeros(n_trials_probe, 1);
    
    t_probe_start = tic;
    parfor ti_main = 1:n_trials_probe
        trial_id = all_trial_ids(ti_main);
        condition = all_trial_conditions{ti_main};
        sf_val = trial_sf_vals(ti_main);
        gain_val = trial_gain_vals(ti_main);
        or_val = trial_or_vals(ti_main);
        
        aligned_obj = trial_objs{ti_main};
        if isempty(aligned_obj)
            continue;
        end
        
        tr_probe_t = aligned_obj.probe_t;
        tr_vel     = aligned_obj.velocity;
        tr_mmask   = aligned_obj.motion_mask;
        tr_smask   = aligned_obj.stationary_mask;
        
        % --- Find continuous motion period ---
        motion_idx = find(tr_mmask);
        if isempty(motion_idx)
            continue;
        end
        motion_start_idx = motion_idx(1);
        motion_end_idx   = motion_idx(end);
        
        t_motion_start = tr_probe_t(motion_start_idx);
        t_motion_end   = tr_probe_t(motion_end_idx);
        motion_dur = t_motion_end - t_motion_start;
        trial_max_dur(ti_main) = motion_dur;
        
        if motion_dur < time_bin_width
            continue;
        end
        
        % Create time bin edges
        bin_edges_t = t_motion_start : time_bin_width : t_motion_end;
        n_tbins = length(bin_edges_t) - 1;
        if n_tbins < 1
            continue;
        end
        
        % --- Vectorised per-bin sample assignment ---
        sample_bin = discretize(tr_probe_t, bin_edges_t);
        sample_bin(isnan(sample_bin)) = 0;
        
        valid_mask   = sample_bin > 0;
        bins_vec     = sample_bin(valid_mask);
        mmask_vec    = tr_mmask(valid_mask);
        absvel_vec   = abs(tr_vel(valid_mask)) .* mmask_vec;
        
        n_samp_per_bin   = accumarray(bins_vec(:), 1,             [n_tbins 1]);
        n_motion_per_bin = accumarray(bins_vec(:), mmask_vec(:),  [n_tbins 1]);
        sum_speed_per_bin = accumarray(bins_vec(:), absvel_vec(:), [n_tbins 1]);
        
        good_bins = find(n_samp_per_bin > 0 & n_motion_per_bin >= n_samp_per_bin * 0.5);
        if isempty(good_bins)
            continue;
        end
        
        mean_speed_per_bin = sum_speed_per_bin(good_bins) ./ n_motion_per_bin(good_bins);
        bin_centres = (bin_edges_t(1:end-1) + bin_edges_t(2:end)) / 2;
        time_in_trial_vec = bin_centres(good_bins)' - t_motion_start;
        
        n_good = length(good_bins);
        switch condition
            case 'T_Vstatic'
                speed_vec = mean_speed_per_bin;
                tf_vec    = zeros(n_good, 1);
                sf_vec    = NaN(n_good, 1);
                or_vec    = NaN(n_good, 1);
                gain_vec  = NaN(n_good, 1);
            case 'V'
                speed_vec = zeros(n_good, 1);
                tf_vec    = gain_val * mean_speed_per_bin;
                sf_vec    = repmat(sf_val, n_good, 1);
                or_vec    = repmat(or_val, n_good, 1);
                gain_vec  = repmat(gain_val, n_good, 1);
            case 'VT'
                speed_vec = mean_speed_per_bin;
                tf_vec    = gain_val * mean_speed_per_bin;
                sf_vec    = repmat(sf_val, n_good, 1);
                or_vec    = repmat(or_val, n_good, 1);
                gain_vec  = repmat(gain_val, n_good, 1);
            otherwise
                speed_vec = mean_speed_per_bin;
                tf_vec    = zeros(n_good, 1);
                sf_vec    = NaN(n_good, 1);
                or_vec    = NaN(n_good, 1);
                gain_vec  = NaN(n_good, 1);
        end
        
        % Collect motion bin data for all clusters
        motion_rows = cell(n_valid_clusters, 1);
        for ci = 1:n_valid_clusters
            st = spike_times_valid{ci};
            spike_counts_all = histcounts(st, bin_edges_t)';
            spike_counts_good = spike_counts_all(good_bins);
            
            motion_rows{ci} = [repmat(valid_cids(ci), n_good, 1), ...
                               repmat(trial_id, n_good, 1), ...
                               speed_vec, tf_vec, sf_vec, or_vec, gain_vec, ...
                               spike_counts_good, time_in_trial_vec, time_in_trial_vec];
        end
        motion_data = vertcat(motion_rows{:});
        motion_conds = repmat({condition}, size(motion_data, 1), 1);
        
        % --- STATIONARY PERIOD DATA ---
        stat_data = [];
        stat_conds = {};
        stationary_before_motion = tr_smask & (1:length(tr_smask))' < motion_start_idx;
        stationary_idx = find(stationary_before_motion);
        
        if ~isempty(stationary_idx)
            stat_end_idx = stationary_idx(end);
            gaps = find(diff(stationary_idx) > 1);
            if ~isempty(gaps)
                stat_start_idx = stationary_idx(gaps(end) + 1);
            else
                stat_start_idx = stationary_idx(1);
            end
            
            t_stat_start = tr_probe_t(stat_start_idx);
            t_stat_end   = tr_probe_t(stat_end_idx);
            stat_dur = t_stat_end - t_stat_start;
            
            if stat_dur >= time_bin_width
                bin_edges_stat = t_stat_start : time_bin_width : t_stat_end;
                n_stat_bins = length(bin_edges_stat) - 1;
                
                if n_stat_bins >= 1
                    n_stat = n_stat_bins;
                    bin_centres_stat = (bin_edges_stat(1:end-1) + bin_edges_stat(2:end)) / 2;
                    time_since_onset_stat = bin_centres_stat' - t_motion_start;
                    
                    stat_rows = cell(n_valid_clusters, 1);
                    for ci = 1:n_valid_clusters
                        st = spike_times_valid{ci};
                        spike_counts_stat = histcounts(st, bin_edges_stat)';
                        
                        stat_rows{ci} = [repmat(valid_cids(ci), n_stat, 1), ...
                                         repmat(trial_id, n_stat, 1), ...
                                         zeros(n_stat, 1), zeros(n_stat, 1), ...
                                         NaN(n_stat, 1), NaN(n_stat, 1), NaN(n_stat, 1), ...
                                         spike_counts_stat, zeros(n_stat, 1), time_since_onset_stat];
                    end
                    stat_data = vertcat(stat_rows{:});
                    stat_conds = repmat({'stationary'}, size(stat_data, 1), 1);
                end
            end
        end
        
        % Combine motion and stationary data
        all_data = [motion_data; stat_data];
        all_conds = [motion_conds; stat_conds];
        
        trial_results{ti_main} = struct('data', all_data, 'conds', {all_conds}, 'pid', pid);
    end
    
    elapsed_probe = toc(t_probe_start);
    fprintf('  Parallel processing done in %.1f s\n', elapsed_probe);
    
    % Track max duration
    max_trial_duration = max(max_trial_duration, max(trial_max_dur));
    
    % Combine results from all trials for this probe
    probe_data = [];
    probe_conds = {};
    probe_pids = {};
    for ti_main = 1:n_trials_probe
        if ~isempty(trial_results{ti_main}) && ~isempty(trial_results{ti_main}.data)
            probe_data = [probe_data; trial_results{ti_main}.data];
            probe_conds = [probe_conds; trial_results{ti_main}.conds];
            probe_pids = [probe_pids; repmat({pid}, size(trial_results{ti_main}.data, 1), 1)];
        end
    end
    
    if ~isempty(probe_data)
        all_probe_tables{probe_i} = table(...
            string(probe_pids), ...
            probe_data(:,1), probe_data(:,2), string(probe_conds), ...
            probe_data(:,3), probe_data(:,4), probe_data(:,5), ...
            probe_data(:,6), probe_data(:,7), probe_data(:,8), ...
            probe_data(:,9), probe_data(:,10), ...
            'VariableNames', {'probe_id', 'cluster_id', 'trial_id', 'condition', ...
                'speed', 'tf', 'sf', 'orientation', 'batch_gain', ...
                'spike_count', 'time_in_trial', 'time_since_onset'});
    end
    
    fprintf('  Probe done: %d trials, %d rows\n', n_trials_probe, size(probe_data, 1));
end

% Concatenate all probe tables
T_master_time = vertcat(all_probe_tables{:});
total_rows_t = height(T_master_time);

fprintf('  Time-bin table: %d rows x %d columns\n', height(T_master_time), width(T_master_time));

% Summarize onset dynamics statistics
onset_stats_motion = T_master_time.time_since_onset(T_master_time.condition ~= "stationary");
onset_stats_stat = T_master_time.time_since_onset(T_master_time.condition == "stationary");
fprintf('  Onset dynamics range:\n');
fprintf('    Motion bins: t_since_onset in [%.2f, %.2f] s\n', min(onset_stats_motion), max(onset_stats_motion));
fprintf('    Stationary bins: t_since_onset in [%.2f, %.2f] s (negative = before onset)\n', min(onset_stats_stat), max(onset_stats_stat));

% Summarize data by condition (including stationary baseline)
cond_counts = groupcounts(T_master_time, 'condition');
fprintf('  Data breakdown by condition:\n');
for ci_cond = 1:height(cond_counts)
    fprintf('    %s: %d bins (%.1f%%)\n', cond_counts.condition(ci_cond), ...
        cond_counts.GroupCount(ci_cond), 100*cond_counts.GroupCount(ci_cond)/height(T_master_time));
end

% Track stationary vs motion for reference
n_stationary = sum(T_master_time.condition == "stationary");
n_motion = height(T_master_time) - n_stationary;
fprintf('  Stationary baseline: %d bins (%.1f%%) | Motion: %d bins (%.1f%%)\n', ...
    n_stationary, 100*n_stationary/height(T_master_time), ...
    n_motion, 100*n_motion/height(T_master_time));

% Temporal basis range (based on max motion duration, not full trial)
trial_duration_range = [0, max_trial_duration];
fprintf('  Max motion duration: %.2f s\n', max_trial_duration);

%% ====================================================================
%  Section 4: Regularisation setup and cluster enumeration
%  ====================================================================

% Get unique clusters (use T_master_time for the time-bin GLM)
unique_clusters = unique(T_master_time(:, {'probe_id', 'cluster_id'}), 'rows');
n_unique_clusters = height(unique_clusters);
fprintf('  Unique clusters: %d\n', n_unique_clusters);

% Ridge regularisation for FullInteraction model (always applied)
lambda_time = 1.0;
fprintf('  Regularisation: lambda=%.1f (applied to FullInteraction model)\n', lambda_time);

%% ====================================================================
%  Section 4b: Basis function visualisation
%  ====================================================================
fprintf('\n--- Plotting basis functions ---\n');

x_speed_fine = linspace(speed_range(1), speed_range(2), 500)';
x_tf_fine    = linspace(tf_range(1),    tf_range(2),    500)';
x_onset_fine = linspace(-0.5, onset_range(2), 500)';  % Include some negative time

B_speed_fine = make_raised_cosine_basis(x_speed_fine, n_speed_bases, speed_range(1), speed_range(2));
B_tf_fine    = make_raised_cosine_basis(x_tf_fine,    n_tf_bases,    tf_range(1),    tf_range(2));
B_onset_fine = make_onset_kernel_basis(x_onset_fine, n_onset_bases, onset_range(2));

fig_basis = figure('Position', [50 50 2200 800], 'Name', 'Basis Functions Overview');

% (a) Speed bases — linear x-axis
subplot(2, 4, 1);
cmap_speed = lines(n_speed_bases);
for bi = 1:n_speed_bases
    plot(x_speed_fine, B_speed_fine(:, bi), 'Color', cmap_speed(bi,:), 'LineWidth', 1.8); hold on;
end
xlabel('Speed (cm/s)'); ylabel('Basis amplitude');
title(sprintf('(a) Speed basis (n=%d)', n_speed_bases));
set(gca, 'box', 'off'); xlim(speed_range);

% (b) Speed bases — log-shifted axis
subplot(2, 4, 2);
epsilon = 0.5;
for bi = 1:n_speed_bases
    plot(log(x_speed_fine + epsilon), B_speed_fine(:, bi), ...
        'Color', cmap_speed(bi,:), 'LineWidth', 1.8); hold on;
end
xlabel('log(Speed + \epsilon)'); ylabel('Basis amplitude');
title('(b) Speed basis (log-shifted)');
set(gca, 'box', 'off');

% (c) TF bases
subplot(2, 4, 3);
cmap_tf = lines(n_tf_bases);
for bi = 1:n_tf_bases
    plot(x_tf_fine, B_tf_fine(:, bi), 'Color', cmap_tf(bi,:), 'LineWidth', 1.8); hold on;
end
xlabel('Temporal frequency (Hz)'); ylabel('Basis amplitude');
title(sprintf('(c) TF basis (n=%d)', n_tf_bases));
set(gca, 'box', 'off'); xlim(tf_range);

% (d) Onset kernel bases (Park et al. 2014 style)
subplot(2, 4, 4);
cmap_onset = lines(n_onset_bases);
for bi = 1:n_onset_bases
    plot(x_onset_fine, B_onset_fine(:, bi), 'Color', cmap_onset(bi,:), 'LineWidth', 1.8); hold on;
end
xline(0, 'k--', 'Motion onset', 'LabelOrientation', 'horizontal', 'LineWidth', 1.5);
xlabel('Time since motion onset (s)'); ylabel('Basis amplitude');
title(sprintf('(d) Onset kernel (n=%d)', n_onset_bases));
set(gca, 'box', 'off'); xlim([-0.5, onset_range(2)]);

% (e) SF dummy coding
subplot(2, 4, 5);
sf_levels = sort(cell2mat(values(sf_values_map)));
n_sf = length(sf_levels);
bar_data_sf = eye(n_sf);
b_sf = bar(1:n_sf, bar_data_sf, 'stacked');
colormap_sf = [0.85 0.85 0.85; lines(n_sf - 1)];
for si = 1:n_sf
    b_sf(si).FaceColor = 'flat';
    b_sf(si).CData = colormap_sf(si,:);
end
set(gca, 'XTick', 1:n_sf, 'XTickLabel', ...
    arrayfun(@(v) sprintf('%.3f', v), sf_levels, 'UniformOutput', false));
xlabel('SF level'); ylabel('Dummy value');
title(sprintf('(e) SF dummy (%d levels)', n_sf));
set(gca, 'box', 'off');

% (f) Orientation dummy coding
subplot(2, 4, 6);
or_levels = sort(cell2mat(values(or_values_map)));
n_or = length(or_levels);
bar_data_or = eye(n_or);
b_or = bar(1:n_or, bar_data_or, 'stacked');
colormap_or = [0.85 0.85 0.85; lines(n_or - 1)];
for oi = 1:n_or
    b_or(oi).FaceColor = 'flat';
    b_or(oi).CData = colormap_or(oi,:);
end
set(gca, 'XTick', 1:n_or, 'XTickLabel', ...
    arrayfun(@(v) sprintf('%.1f°', rad2deg(v)), or_levels, 'UniformOutput', false));
xlabel('Orientation'); ylabel('Dummy value');
title(sprintf('(f) OR dummy (%d levels)', n_or));
set(gca, 'box', 'off');

% (g) Speed x TF interaction (tensor product, 2D contour)
subplot(2, 4, 7);
n_grid = 80;
spd_grid = linspace(speed_range(1), speed_range(2), n_grid)';
tf_grid  = linspace(tf_range(1),    tf_range(2),    n_grid)';
B_spd_g = make_raised_cosine_basis(spd_grid, n_speed_bases, speed_range(1), speed_range(2));
B_tf_g  = make_raised_cosine_basis(tf_grid,  n_tf_bases,    tf_range(1),    tf_range(2));

[SPD_mesh, TF_mesh] = meshgrid(spd_grid, tf_grid);
n_examples = min(4, n_speed_bases * n_tf_bases);
example_pairs = [1 1; n_speed_bases 1; 1 n_tf_bases; n_speed_bases n_tf_bases];
example_pairs = example_pairs(1:n_examples, :);
cmap_inter = lines(n_examples);
hold on;
leg_entries = cell(1, n_examples);
for ei = 1:n_examples
    si = example_pairs(ei, 1);
    ti = example_pairs(ei, 2);
    Z = B_tf_g(:, ti) * B_spd_g(:, si)';
    contour(SPD_mesh, TF_mesh, Z, [0.25, 0.5, 0.75], 'LineColor', cmap_inter(ei,:), 'LineWidth', 1.5);
    leg_entries{ei} = sprintf('B_{spd,%d} x B_{tf,%d}', si, ti);
end
xlabel('Speed (cm/s)'); ylabel('TF (Hz)');
title(sprintf('(g) Speed x TF interaction (%d bases)', n_speed_bases * n_tf_bases));
legend(leg_entries, 'Location', 'eastoutside', 'FontSize', 7);
set(gca, 'box', 'off');

% (h) Onset kernel: example weighted sums (shows flexibility of basis)
subplot(2, 4, 8);
% Show several example kernels formed by different weight combinations
n_examples_onset = 5;
t_kernel = linspace(0, onset_range(2), 200)';
B_kernel = make_onset_kernel_basis(t_kernel, n_onset_bases, onset_range(2));

% Generate diverse weight patterns
rng(42);  % reproducible
example_weights = [
    1, 0.5, 0.2, 0.05, 0.01, 0;        % fast decay (early transient)
    0, 0.2, 0.5, 0.8, 0.5, 0.2;        % delayed peak
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5;      % sustained
    1, -0.3, 0.1, 0, 0, 0;             % overshoot
    0.2, 0.4, 0.6, 0.8, 1, 0.8;        % slow rise
];
example_weights = example_weights(:, 1:n_onset_bases);  % trim to actual basis count

cmap_ex = lines(n_examples_onset);
leg_ex = cell(1, n_examples_onset);
labels_ex = {'Fast decay', 'Delayed peak', 'Sustained', 'Overshoot', 'Slow rise'};
hold on;
for ei = 1:n_examples_onset
    w = example_weights(ei, :)';
    kernel_shape = B_kernel * w;
    kernel_shape = kernel_shape / max(abs(kernel_shape));  % normalize for display
    plot(t_kernel, kernel_shape, 'Color', cmap_ex(ei,:), 'LineWidth', 1.8);
    leg_ex{ei} = labels_ex{ei};
end
hold off;
xline(0, 'k--', 'LineWidth', 1);
xlabel('Time since motion onset (s)'); ylabel('Kernel amplitude (normalized)');
title('(h) Example onset kernels (weighted sums)');
legend(leg_ex, 'Location', 'northeast', 'FontSize', 7);
set(gca, 'box', 'off'); xlim([0, onset_range(2)]);

sgtitle('GLM Basis Functions (v9: with Onset Dynamics)', 'FontSize', 14, 'FontWeight', 'bold');

if save_figs
    print(fig_basis, fullfile(ctl.figs.curr_dir, 'basis_functions_overview.png'), '-dpng', '-r300');
    fprintf('  Saved basis_functions_overview.png\n');
end

%% ====================================================================
%  Section 5: Model fitting and cross-validation (Hierarchical Forward Selection)
%  ====================================================================
fprintf('\n--- Hierarchical Forward Selection GLM ---\n');
fprintf('    Phase 1: Main effects (Speed, TF, SF, OR)\n');
fprintf('    Phase 2: Interactions (only if both parents selected)\n\n');

glm_types = {'time'};
glm_type_labels = {'Time-bin'};

% Models to fit for visualization (comparison)
visualization_models = {'Null', 'Additive', 'FullInteraction'};  % Selected added dynamically
n_vis_models = length(visualization_models);

% Storage for GLM results
for gt = 1:length(glm_types)
    gt_tag = glm_types{gt};
    
    results.(gt_tag) = struct();
    results.(gt_tag).probe_id = cell(n_unique_clusters, 1);
    results.(gt_tag).cluster_id = zeros(n_unique_clusters, 1);
    results.(gt_tag).n_trials = zeros(n_unique_clusters, 1);
    results.(gt_tag).n_bins_total = zeros(n_unique_clusters, 1);
    results.(gt_tag).n_spikes_total = zeros(n_unique_clusters, 1);
    
    % --- Forward selection results ---
    results.(gt_tag).selected_vars = cell(n_unique_clusters, 1);         % Cell array of selected variable names
    results.(gt_tag).selected_vars_str = cell(n_unique_clusters, 1);     % Comma-separated string
    results.(gt_tag).n_selected_vars = zeros(n_unique_clusters, 1);      % Number of selected variables
    results.(gt_tag).selection_rounds = zeros(n_unique_clusters, 1);     % Rounds until convergence
    results.(gt_tag).selection_history = cell(n_unique_clusters, 1);     % Full selection history
    
    % --- CV performance for key models ---
    results.(gt_tag).Null_cv_bps = zeros(n_unique_clusters, 1);
    results.(gt_tag).Selected_cv_bps = zeros(n_unique_clusters, 1);
    results.(gt_tag).Additive_cv_bps = zeros(n_unique_clusters, 1);
    results.(gt_tag).FullInteraction_cv_bps = zeros(n_unique_clusters, 1);
    
    % --- Model fit metrics for all comparison models ---
    for mi = 1:n_vis_models
        ml = visualization_models{mi};
        results.(gt_tag).([ml '_aic']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_bic']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_train_bps']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_deviance']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_n_params']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_dispersion']) = zeros(n_unique_clusters, 1);
    end
    % Also for Selected model
    results.(gt_tag).Selected_aic = zeros(n_unique_clusters, 1);
    results.(gt_tag).Selected_bic = zeros(n_unique_clusters, 1);
    results.(gt_tag).Selected_train_bps = zeros(n_unique_clusters, 1);
    results.(gt_tag).Selected_deviance = zeros(n_unique_clusters, 1);
    results.(gt_tag).Selected_n_params = zeros(n_unique_clusters, 1);
    results.(gt_tag).Selected_dispersion = zeros(n_unique_clusters, 1);
    
    % --- Delta comparisons ---
    results.(gt_tag).delta_bps_selected_vs_null = zeros(n_unique_clusters, 1);   % Selected - Null
    results.(gt_tag).delta_bps_additive_vs_null = zeros(n_unique_clusters, 1);   % Additive - Null
    results.(gt_tag).delta_bps_interaction = zeros(n_unique_clusters, 1);        % FullInteraction - Additive
    results.(gt_tag).delta_bps_selected_vs_additive = zeros(n_unique_clusters, 1); % Selected - Additive (should be ~0 or negative if selection worked)
    
    % --- Classification flags (based on forward selection) ---
    results.(gt_tag).is_spprofileeed_tuned = false(n_unique_clusters, 1);        % Speed in selected_vars
    results.(gt_tag).is_tf_tuned = false(n_unique_clusters, 1);           % TF in selected_vars
    results.(gt_tag).is_sf_tuned = false(n_unique_clusters, 1);           % SF in selected_vars
    results.(gt_tag).is_or_tuned = false(n_unique_clusters, 1);           % OR in selected_vars
    results.(gt_tag).has_interaction = false(n_unique_clusters, 1);       % Any interaction in selected_vars
    results.(gt_tag).has_spprofileeed_x_tf = false(n_unique_clusters, 1);
    results.(gt_tag).has_spprofileeed_x_sf = false(n_unique_clusters, 1);
    results.(gt_tag).has_spprofileeed_x_or = false(n_unique_clusters, 1);
    results.(gt_tag).has_tf_x_sf = false(n_unique_clusters, 1);
    results.(gt_tag).has_tf_x_or = false(n_unique_clusters, 1);
    results.(gt_tag).has_sf_x_or = false(n_unique_clusters, 1);
    
    % Backward compatibility
    results.(gt_tag).winning_model = cell(n_unique_clusters, 1);
    results.(gt_tag).has_significant_interaction = false(n_unique_clusters, 1);
end

% Storage for coefficients and predictions
cluster_predictions = struct();
cluster_predictions.time = cell(n_unique_clusters, 1);

% Storage for selection history (for CSV export)
all_selection_history = cell(n_unique_clusters, 1);  % Cell array per cluster

% Storage for coefficients (cell array per cluster)
all_coefficients_per_cluster = cell(n_unique_clusters, 1);

% Prepare data for parfor: extract pid and cid arrays
unique_pids = unique_clusters.probe_id;
unique_cids = unique_clusters.cluster_id;

fprintf('\n--- Fitting GLMs in parallel (%d clusters) ---\n', n_unique_clusters);

% === PARALLEL GLM FITTING ===
parfor_results = cell(n_unique_clusters, 1);

parfor ci = 1:n_unique_clusters
    pid = unique_pids(ci);
    cid = unique_cids(ci);
    
    gt_tag = 'time';  % Only time-bin GLM now
    
    % Use time-bin table
    T_all = T_master_time;
    lambda_reg = lambda_time;
    
    idx = T_all.probe_id == pid & T_all.cluster_id == cid;
    T_cluster = T_all(idx, :);
    
    % Initialize result struct for this cluster
    cl_result = struct();
    cl_result.probe_id = char(pid);
    cl_result.cluster_id = cid;
    cl_result.skipped = false;
    
    if height(T_cluster) < 10
        cl_result.skipped = true;
        cl_result.winning_model = 'N/A';
        cl_result.selected_vars = {};
        cl_result.selected_vars_str = '';
        parfor_results{ci} = cl_result;
        continue;
    end
    
    cl_result.n_trials = length(unique(T_cluster.trial_id));
    cl_result.n_bins_total = height(T_cluster);
    cl_result.n_spikes_total = sum(double(T_cluster.spike_count));

    % Prepare features
    speed_v = T_cluster.speed;
    tf_v = T_cluster.tf;
    sf_v = T_cluster.sf; sf_v(isnan(sf_v)) = 0;
    or_v = T_cluster.orientation; or_v(isnan(or_v)) = 0;
    time_since_onset_v = T_cluster.time_since_onset;
    y = double(T_cluster.spike_count);
    
    % Offset: constant offset (all bins same width)
    offset = log(time_bin_width) * ones(height(T_cluster), 1);
    
    % --- Create CV folds (trial-level, condition-stratified) ---
    unique_trials = unique(T_cluster(:, {'trial_id', 'condition'}), 'rows');
    n_tr = height(unique_trials);
    
    trial_fold = zeros(n_tr, 1);
    conditions = unique(unique_trials.condition);
    for cond_i = 1:length(conditions)
        c_idx = find(unique_trials.condition == conditions(cond_i));
        n_c = length(c_idx);
        perm = randperm(n_c);
        trial_fold(c_idx(perm)) = mod((1:n_c)' - 1, n_cv_folds) + 1;
    end
    
    % Vectorised fold assignment via unique-group mapping
    [~, ~, row_to_trial] = unique(T_cluster(:, {'trial_id', 'condition'}), 'rows');
    fold_ids = trial_fold(row_to_trial);
    
    % --- Pre-compute basis matrices once per cluster (shared across models) ---
    B_speed_cl = make_raised_cosine_basis(speed_v, n_speed_bases, speed_range(1), speed_range(2));
    B_tf_cl    = make_raised_cosine_basis(tf_v, n_tf_bases, tf_range(1), tf_range(2));
    B_onset_cl = make_onset_kernel_basis(time_since_onset_v, n_onset_bases, onset_range(2));
    
    % Get reference levels for SF and OR (for consistent dummy coding)
    sf_ref_levels = sort(unique(sf_v(sf_v ~= 0)));
    or_ref_levels = sort(unique(or_v(or_v ~= 0)));
    
    % =====================================================================
    % FORWARD SELECTION: Find optimal variable set for this neuron
    % =====================================================================
    [selected_vars, selection_history, selected_cv_bps, null_cv_bps] = forward_select_model(...
        B_speed_cl, B_tf_cl, B_onset_cl, sf_v, or_v, y, offset, fold_ids, ...
        selection_threshold, all_candidates, sf_ref_levels, or_ref_levels);
    
    % Store selection results
    cl_result.selected_vars = selected_vars;
    cl_result.selected_vars_str = strjoin(selected_vars, ', ');
    cl_result.n_selected_vars = length(selected_vars);
    cl_result.selection_rounds = length(selection_history);
    cl_result.selection_history = selection_history;
    cl_result.Null_cv_bps = null_cv_bps;
    cl_result.Selected_cv_bps = selected_cv_bps;
    
    % Store selection history for CSV export
    cl_selection_history = cell(length(selection_history), 1);
    for ri = 1:length(selection_history)
        sh = selection_history(ri);
        cl_selection_history{ri} = {char(pid), cid, ri, sh.best_candidate, ...
            sh.delta_bps, sh.added, sh.cv_bps_after};
    end
    cl_result.selection_history_csv = cl_selection_history;
    
    % --- Classification based on forward selection ---
    cl_result.is_spprofileeed_tuned = ismember('Speed', selected_vars);
    cl_result.is_tf_tuned = ismember('TF', selected_vars);
    cl_result.is_sf_tuned = ismember('SF', selected_vars);
    cl_result.is_or_tuned = ismember('OR', selected_vars);
    cl_result.has_spprofileeed_x_tf = ismember('Speed_x_TF', selected_vars);
    cl_result.has_spprofileeed_x_sf = ismember('Speed_x_SF', selected_vars);
    cl_result.has_spprofileeed_x_or = ismember('Speed_x_OR', selected_vars);
    cl_result.has_tf_x_sf = ismember('TF_x_SF', selected_vars);
    cl_result.has_tf_x_or = ismember('TF_x_OR', selected_vars);
    cl_result.has_sf_x_or = ismember('SF_x_OR', selected_vars);
    cl_result.has_interaction = any(contains(selected_vars, '_x_'));
    
    % =====================================================================
    % FIT COMPARISON MODELS: Null, Selected, Additive, FullInteraction
    % =====================================================================
    models_to_fit = {'Null', 'Selected', 'Additive', 'FullInteraction'};
    cl_result.models = struct();
    cl_coefficients = {};
    
    for mi = 1:length(models_to_fit)
        ml = models_to_fit{mi};
        
        % Ridge for FullInteraction
        if strcmp(ml, 'FullInteraction')
            lambda = lambda_reg;
        else
            lambda = 0;
        end
        
        % Build design matrix
        if strcmp(ml, 'Selected')
            [X, col_names] = assemble_design_matrix_selected(B_speed_cl, B_tf_cl, B_onset_cl, ...
                sf_v, or_v, selected_vars, sf_ref_levels, or_ref_levels);
        else
            [X, col_names] = assemble_design_matrix(B_speed_cl, B_tf_cl, B_onset_cl, ...
                sf_v, or_v, ml);
        end
        
        % Check dimensions
        if size(X, 2) >= height(T_cluster)
            cl_result.models.(ml).aic = Inf;
            cl_result.models.(ml).bic = Inf;
            cl_result.models.(ml).cv_bps = -Inf;
            cl_result.models.(ml).train_bps = -Inf;
            cl_result.models.(ml).deviance = Inf;
            cl_result.models.(ml).n_params = size(X, 2);
            cl_result.models.(ml).dispersion = NaN;
            cl_result.models.(ml).skipped = true;
            continue;
        end
        
        cl_result.models.(ml).skipped = false;
        
        % Fit on full data
        res = fit_poisson_glm(X, y, offset, lambda);
        
        cl_result.models.(ml).aic = res.aic;
        cl_result.models.(ml).bic = res.bic;
        cl_result.models.(ml).deviance = res.deviance;
        cl_result.models.(ml).n_params = res.n_params;
        cl_result.models.(ml).dispersion = res.dispersion;
        
        if cl_result.n_spikes_total > 0
            cl_result.models.(ml).train_bps = ...
                (res.log_likelihood / cl_result.n_spikes_total) / log(2);
        else
            cl_result.models.(ml).train_bps = NaN;
        end
        
        % Cross-validate (except for Selected which we already have from forward selection)
        if ~strcmp(ml, 'Selected')
            [~, cv_bps, ~, cv_predicted] = cross_validate_glm(X, y, offset, fold_ids, lambda);
            cl_result.models.(ml).cv_bps = cv_bps;
        else
            % For Selected, refit to get CV predictions for plotting
            [~, ~, ~, cv_predicted] = cross_validate_glm(X, y, offset, fold_ids, lambda);
            cl_result.models.(ml).cv_bps = selected_cv_bps;
        end
        
        % Store coefficients
        n_beta = length(res.beta);
        for ki = 1:n_beta
            cl_coefficients{end+1} = {char(pid), cid, gt_tag, ml, col_names{ki}, ...
                res.beta(ki), res.se(ki)}; %#ok<AGROW>
        end
        
        % Store predictions
        cl_result.models.(ml).predicted_count = res.predicted_count;
        cl_result.models.(ml).cv_predicted_count = cv_predicted;
        cl_result.models.(ml).predicted_fr = res.predicted_count ./ exp(offset);
        cl_result.models.(ml).cv_predicted_fr = cv_predicted ./ exp(offset);
        cl_result.models.(ml).pearson_residuals = res.pearson_residuals;
        cl_result.models.(ml).beta = res.beta;
        cl_result.models.(ml).se = res.se;
        cl_result.models.(ml).col_names = col_names;
    end
    
    cl_result.coefficients = cl_coefficients;
    
    % --- Compute delta comparisons ---
    cl_result.delta_bps_selected_vs_null = cl_result.Selected_cv_bps - cl_result.Null_cv_bps;
    if isfield(cl_result.models, 'Additive') && ~cl_result.models.Additive.skipped
        cl_result.Additive_cv_bps = cl_result.models.Additive.cv_bps;
        cl_result.delta_bps_additive_vs_null = cl_result.Additive_cv_bps - cl_result.Null_cv_bps;
        cl_result.delta_bps_selected_vs_additive = cl_result.Selected_cv_bps - cl_result.Additive_cv_bps;
    else
        cl_result.Additive_cv_bps = NaN;
        cl_result.delta_bps_additive_vs_null = NaN;
        cl_result.delta_bps_selected_vs_additive = NaN;
    end
    if isfield(cl_result.models, 'FullInteraction') && ~cl_result.models.FullInteraction.skipped
        cl_result.FullInteraction_cv_bps = cl_result.models.FullInteraction.cv_bps;
        cl_result.delta_bps_interaction = cl_result.FullInteraction_cv_bps - cl_result.Additive_cv_bps;
    else
        cl_result.FullInteraction_cv_bps = NaN;
        cl_result.delta_bps_interaction = NaN;
    end
    
    % Backward compatibility
    cl_result.has_significant_interaction = cl_result.delta_bps_interaction > selection_threshold;
    if ~isempty(selected_vars)
        cl_result.winning_model = 'Selected';
    else
        cl_result.winning_model = 'Null';
    end
    
    parfor_results{ci} = cl_result;
end

fprintf('Parallel GLM fitting complete. Collecting results...\n');

% === Collect results from parfor into results struct ===
all_selection_history_flat = {};
all_coefficients = {};

for ci = 1:n_unique_clusters
    cl_result = parfor_results{ci};
    gt_tag = 'time';
    
    results.(gt_tag).probe_id{ci} = cl_result.probe_id;
    results.(gt_tag).cluster_id(ci) = cl_result.cluster_id;
    
    if cl_result.skipped
        results.(gt_tag).winning_model{ci} = cl_result.winning_model;
        results.(gt_tag).selected_vars{ci} = {};
        results.(gt_tag).selected_vars_str{ci} = '';
        continue;
    end
    
    results.(gt_tag).n_trials(ci) = cl_result.n_trials;
    results.(gt_tag).n_bins_total(ci) = cl_result.n_bins_total;
    results.(gt_tag).n_spikes_total(ci) = cl_result.n_spikes_total;
    results.(gt_tag).selected_vars{ci} = cl_result.selected_vars;
    results.(gt_tag).selected_vars_str{ci} = cl_result.selected_vars_str;
    results.(gt_tag).n_selected_vars(ci) = cl_result.n_selected_vars;
    results.(gt_tag).selection_rounds(ci) = cl_result.selection_rounds;
    results.(gt_tag).selection_history{ci} = cl_result.selection_history;
    results.(gt_tag).Null_cv_bps(ci) = cl_result.Null_cv_bps;
    results.(gt_tag).Selected_cv_bps(ci) = cl_result.Selected_cv_bps;
    results.(gt_tag).Additive_cv_bps(ci) = cl_result.Additive_cv_bps;
    results.(gt_tag).FullInteraction_cv_bps(ci) = cl_result.FullInteraction_cv_bps;
    
    results.(gt_tag).is_spprofileeed_tuned(ci) = cl_result.is_spprofileeed_tuned;
    results.(gt_tag).is_tf_tuned(ci) = cl_result.is_tf_tuned;
    results.(gt_tag).is_sf_tuned(ci) = cl_result.is_sf_tuned;
    results.(gt_tag).is_or_tuned(ci) = cl_result.is_or_tuned;
    results.(gt_tag).has_spprofileeed_x_tf(ci) = cl_result.has_spprofileeed_x_tf;
    results.(gt_tag).has_spprofileeed_x_sf(ci) = cl_result.has_spprofileeed_x_sf;
    results.(gt_tag).has_spprofileeed_x_or(ci) = cl_result.has_spprofileeed_x_or;
    results.(gt_tag).has_tf_x_sf(ci) = cl_result.has_tf_x_sf;
    results.(gt_tag).has_tf_x_or(ci) = cl_result.has_tf_x_or;
    results.(gt_tag).has_sf_x_or(ci) = cl_result.has_sf_x_or;
    results.(gt_tag).has_interaction(ci) = cl_result.has_interaction;
    
    results.(gt_tag).delta_bps_selected_vs_null(ci) = cl_result.delta_bps_selected_vs_null;
    results.(gt_tag).delta_bps_additive_vs_null(ci) = cl_result.delta_bps_additive_vs_null;
    results.(gt_tag).delta_bps_interaction(ci) = cl_result.delta_bps_interaction;
    results.(gt_tag).delta_bps_selected_vs_additive(ci) = cl_result.delta_bps_selected_vs_additive;
    results.(gt_tag).has_significant_interaction(ci) = cl_result.has_significant_interaction;
    results.(gt_tag).winning_model{ci} = cl_result.winning_model;
    
    % Model-specific results
    models_to_fit = {'Null', 'Selected', 'Additive', 'FullInteraction'};
    for mi = 1:length(models_to_fit)
        ml = models_to_fit{mi};
        if isfield(cl_result.models, ml)
            results.(gt_tag).([ml '_aic'])(ci) = cl_result.models.(ml).aic;
            results.(gt_tag).([ml '_bic'])(ci) = cl_result.models.(ml).bic;
            results.(gt_tag).([ml '_deviance'])(ci) = cl_result.models.(ml).deviance;
            results.(gt_tag).([ml '_n_params'])(ci) = cl_result.models.(ml).n_params;
            results.(gt_tag).([ml '_dispersion'])(ci) = cl_result.models.(ml).dispersion;
            results.(gt_tag).([ml '_train_bps'])(ci) = cl_result.models.(ml).train_bps;
            if ~strcmp(ml, 'Selected')
                results.(gt_tag).([ml '_cv_bps'])(ci) = cl_result.models.(ml).cv_bps;
            end
            
            if ~cl_result.models.(ml).skipped
                cluster_predictions.(gt_tag){ci}.(ml).predicted_count = cl_result.models.(ml).predicted_count;
                cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_count = cl_result.models.(ml).cv_predicted_count;
                cluster_predictions.(gt_tag){ci}.(ml).predicted_fr = cl_result.models.(ml).predicted_fr;
                cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr = cl_result.models.(ml).cv_predicted_fr;
                cluster_predictions.(gt_tag){ci}.(ml).pearson_residuals = cl_result.models.(ml).pearson_residuals;
                cluster_predictions.(gt_tag){ci}.(ml).beta = cl_result.models.(ml).beta;
                cluster_predictions.(gt_tag){ci}.(ml).se = cl_result.models.(ml).se;
                cluster_predictions.(gt_tag){ci}.(ml).col_names = cl_result.models.(ml).col_names;
            end
        end
    end
    
    % Collect selection history
    all_selection_history_flat = [all_selection_history_flat; cl_result.selection_history_csv]; %#ok<AGROW>
    
    % Collect coefficients
    all_coefficients = [all_coefficients; cl_result.coefficients(:)]; %#ok<AGROW>
end

% Rename for backward compatibility
all_selection_history = all_selection_history_flat;

fprintf('\n--- Forward selection complete ---\n');

%% ====================================================================
%  Section 6: Population summary
%  ====================================================================
fprintf('\n--- Population summary ---\n');

for gt = 1:length(glm_types)
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    fprintf('\n=== %s GLM (Forward Selection) ===\n', gt_label);
    
    % Model complexity distribution
    n_vars = results.(gt_tag).n_selected_vars;
    fprintf('Model complexity distribution:\n');
    fprintf('  0 vars (Null only): %d (%.1f%%)\n', sum(n_vars == 0), 100*sum(n_vars == 0)/n_unique_clusters);
    fprintf('  1 var:  %d (%.1f%%)\n', sum(n_vars == 1), 100*sum(n_vars == 1)/n_unique_clusters);
    fprintf('  2 vars: %d (%.1f%%)\n', sum(n_vars == 2), 100*sum(n_vars == 2)/n_unique_clusters);
    fprintf('  3 vars: %d (%.1f%%)\n', sum(n_vars == 3), 100*sum(n_vars == 3)/n_unique_clusters);
    fprintf('  4+ vars: %d (%.1f%%)\n', sum(n_vars >= 4), 100*sum(n_vars >= 4)/n_unique_clusters);
    
    % Variable inclusion rates (main effects)
    fprintf('Variable inclusion rates (forward selection):\n');
    fprintf('  Speed:    %d (%.1f%%)\n', sum(results.(gt_tag).is_spprofileeed_tuned), ...
        100*sum(results.(gt_tag).is_spprofileeed_tuned)/n_unique_clusters);
    fprintf('  TF:       %d (%.1f%%)\n', sum(results.(gt_tag).is_tf_tuned), ...
        100*sum(results.(gt_tag).is_tf_tuned)/n_unique_clusters);
    fprintf('  SF:       %d (%.1f%%)\n', sum(results.(gt_tag).is_sf_tuned), ...
        100*sum(results.(gt_tag).is_sf_tuned)/n_unique_clusters);
    fprintf('  OR:       %d (%.1f%%)\n', sum(results.(gt_tag).is_or_tuned), ...
        100*sum(results.(gt_tag).is_or_tuned)/n_unique_clusters);
    
    % Interaction inclusion rates
    fprintf('Interaction inclusion rates:\n');
    fprintf('  Speed x TF: %d (%.1f%%)\n', sum(results.(gt_tag).has_spprofileeed_x_tf), ...
        100*sum(results.(gt_tag).has_spprofileeed_x_tf)/n_unique_clusters);
    fprintf('  Speed x SF: %d (%.1f%%)\n', sum(results.(gt_tag).has_spprofileeed_x_sf), ...
        100*sum(results.(gt_tag).has_spprofileeed_x_sf)/n_unique_clusters);
    fprintf('  Speed x OR: %d (%.1f%%)\n', sum(results.(gt_tag).has_spprofileeed_x_or), ...
        100*sum(results.(gt_tag).has_spprofileeed_x_or)/n_unique_clusters);
    fprintf('  TF x SF:    %d (%.1f%%)\n', sum(results.(gt_tag).has_tf_x_sf), ...
        100*sum(results.(gt_tag).has_tf_x_sf)/n_unique_clusters);
    fprintf('  TF x OR:    %d (%.1f%%)\n', sum(results.(gt_tag).has_tf_x_or), ...
        100*sum(results.(gt_tag).has_tf_x_or)/n_unique_clusters);
    fprintf('  SF x OR:    %d (%.1f%%)\n', sum(results.(gt_tag).has_sf_x_or), ...
        100*sum(results.(gt_tag).has_sf_x_or)/n_unique_clusters);
    fprintf('  Any interaction: %d (%.1f%%)\n', sum(results.(gt_tag).has_interaction), ...
        100*sum(results.(gt_tag).has_interaction)/n_unique_clusters);
    
    % CV performance summary
    fprintf('CV performance (bits/spike):\n');
    fprintf('  Null:           median=%.4f\n', median(results.(gt_tag).Null_cv_bps));
    fprintf('  Selected:       median=%.4f\n', median(results.(gt_tag).Selected_cv_bps));
    fprintf('  Additive:       median=%.4f\n', median(results.(gt_tag).Additive_cv_bps));
    fprintf('  FullInteraction: median=%.4f\n', median(results.(gt_tag).FullInteraction_cv_bps));
    
    % Delta comparisons
    fprintf('Delta CV bps:\n');
    fprintf('  Selected - Null:     median=%.4f\n', median(results.(gt_tag).delta_bps_selected_vs_null));
    fprintf('  Additive - Null:     median=%.4f\n', median(results.(gt_tag).delta_bps_additive_vs_null));
    fprintf('  FullInt - Additive:  median=%.4f\n', median(results.(gt_tag).delta_bps_interaction));
end

% Summary statistics for paper
fprintf('\n=== Summary Statistics for Paper ===\n');
n_speed_tuned = sum(results.time.is_spprofileeed_tuned);
n_tf_tuned = sum(results.time.is_tf_tuned);
n_sf_tuned = sum(results.time.is_sf_tuned);
n_or_tuned = sum(results.time.is_or_tuned);
n_has_interaction = sum(results.time.has_interaction);

fprintf('Variable selection rates (forward selection, threshold=%.3f bps):\n', selection_threshold);
fprintf('  Speed tuned: %d/%d (%.1f%%)\n', n_speed_tuned, n_unique_clusters, 100*n_speed_tuned/n_unique_clusters);
fprintf('  TF tuned:    %d/%d (%.1f%%)\n', n_tf_tuned, n_unique_clusters, 100*n_tf_tuned/n_unique_clusters);
fprintf('  SF tuned:    %d/%d (%.1f%%)\n', n_sf_tuned, n_unique_clusters, 100*n_sf_tuned/n_unique_clusters);
fprintf('  OR tuned:    %d/%d (%.1f%%)\n', n_or_tuned, n_unique_clusters, 100*n_or_tuned/n_unique_clusters);
fprintf('  Has interaction: %d/%d (%.1f%%)\n', n_has_interaction, n_unique_clusters, 100*n_has_interaction/n_unique_clusters);

%% ====================================================================
%  Section 7: Population visualisations
%  ====================================================================
fprintf('\n--- Generating population figures ---\n');

% --- Figure 2: Forward Selection Summary (Updated to match printed output) ---
% Panel 1: Model complexity distribution (bar chart showing 0, 1, 2, 3, 4+ vars)
% Panel 2: Variable inclusion rates (main effects from forward selection)
% Panel 3: Interaction inclusion rates (breakdown by interaction type)
% Panel 4: Dominant modulator per neuron (pie chart from component importance)
fig2 = figure('Position', [50 50 1800 500], 'Name', 'Forward Selection Summary');

for gt = 1:1  % Only time-bin GLM now
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    n_total = n_unique_clusters;
    
    % Get forward selection results
    n_vars = results.(gt_tag).n_selected_vars;
    is_spprofiled = results.(gt_tag).is_spprofileeed_tuned;
    is_tf = results.(gt_tag).is_tf_tuned;
    is_sf = results.(gt_tag).is_sf_tuned;
    is_or = results.(gt_tag).is_or_tuned;
    has_int = results.(gt_tag).has_interaction;
    
    % Get interaction breakdown
    has_spprofiled_tf = results.(gt_tag).has_spprofileeed_x_tf;
    has_spprofiled_sf = results.(gt_tag).has_spprofileeed_x_sf;
    has_spprofiled_or = results.(gt_tag).has_spprofileeed_x_or;
    has_tf_sf = results.(gt_tag).has_tf_x_sf;
    has_tf_or = results.(gt_tag).has_tf_x_or;
    has_sf_or = results.(gt_tag).has_sf_x_or;
    
    % --- Panel 1: Model Complexity Distribution (Stacked Bar by Model Type) ---
    subplot(1, 3, 1);
    
    % Define colors for main effects (used consistently throughout this figure)
    % Color scheme: Speed=green, SF=yellow, TF=orange, OR=red, Interactions=blue
    color_speed = [0.17 0.63 0.17];   % Green
    color_tf = [1.0 0.50 0.05];       % Orange
    color_sf = [0.95 0.85 0.10];      % Yellow
    color_or = [0.84 0.15 0.16];      % Red
    color_null = [0.75 0.75 0.75];    % Light gray for Null
    color_4plus = [0.35 0.35 0.35];   % Dark gray for 4+
    color_int = [0.20 0.40 0.80];     % Blue for interactions
    
    % Get selected variables for each neuron
    selected_vars_list = results.(gt_tag).selected_vars;
    
    % Helper function to darken color for interaction terms
    darken = @(c) max(c * 0.65, 0);
    
    % Interaction term names and their colors (blue variations)
    interaction_info = containers.Map();
    interaction_info('Speed_x_TF') = struct('label', 'SxT', 'color', [0.12 0.47 0.71]);    % Medium blue
    interaction_info('Speed_x_SF') = struct('label', 'SxSF', 'color', [0.26 0.58 0.85]);   % Light blue
    interaction_info('Speed_x_OR') = struct('label', 'SxO', 'color', [0.05 0.30 0.55]);    % Dark blue
    interaction_info('TF_x_SF') = struct('label', 'TxSF', 'color', [0.40 0.60 0.80]);      % Sky blue
    interaction_info('TF_x_OR') = struct('label', 'TxO', 'color', [0.15 0.40 0.65]);       % Steel blue
    interaction_info('SF_x_OR') = struct('label', 'SFxO', 'color', [0.30 0.50 0.75]);      % Slate blue
    
    % Main effect info (consistent colors)
    main_effect_info = containers.Map();
    main_effect_info('Speed') = struct('label', 'Spd', 'color', color_speed);
    main_effect_info('TF') = struct('label', 'TF', 'color', color_tf);
    main_effect_info('SF') = struct('label', 'SF', 'color', color_sf);
    main_effect_info('OR') = struct('label', 'OR', 'color', color_or);
    
    all_var_names = {'Speed', 'TF', 'SF', 'OR', 'Speed_x_TF', 'Speed_x_SF', 'Speed_x_OR', 'TF_x_SF', 'TF_x_OR', 'SF_x_OR'};
    
    % Build a catalog of unique model types found in the data
    % Key = sorted string of selected vars, Value = count
    model_type_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
    model_type_complexity = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    for ni = 1:length(selected_vars_list)
        sv = selected_vars_list{ni};
        if isempty(sv)
            sv = {};
        end
        
        % Total complexity = total number of selected variables
        n_total_vars = length(sv);
        
        % Create canonical key (sorted variable names)
        if isempty(sv)
            model_key = 'Null';
        else
            model_key = strjoin(sort(sv), '+');
        end
        
        % Count this model type
        if isKey(model_type_counts, model_key)
            model_type_counts(model_key) = model_type_counts(model_key) + 1;
        else
            model_type_counts(model_key) = 1;
            model_type_complexity(model_key) = n_total_vars;
        end
    end
    
    % Get unique model types and their counts
    model_keys = keys(model_type_counts);
    n_model_types = length(model_keys);
    
    % Assign colors based on model composition
    model_colors = zeros(n_model_types, 3);
    model_labels = cell(n_model_types, 1);
    model_complexity_arr = zeros(n_model_types, 1);
    model_counts = zeros(n_model_types, 1);
    model_main_effects = cell(n_model_types, 1);  % Track constituent main effects for striping
    
    for mi = 1:n_model_types
        mk = model_keys{mi};
        model_counts(mi) = model_type_counts(mk);
        model_complexity_arr(mi) = model_type_complexity(mk);
        
        if strcmp(mk, 'Null')
            model_colors(mi, :) = color_null;
            model_labels{mi} = 'Null';
            model_main_effects{mi} = {};
        else
            % Parse variables in model
            vars_in_model = strsplit(mk, '+');
            main_effects = vars_in_model(ismember(vars_in_model, {'Speed', 'TF', 'SF', 'OR'}));
            interactions = vars_in_model(~ismember(vars_in_model, {'Speed', 'TF', 'SF', 'OR'}));
            model_main_effects{mi} = main_effects;  % Store for striping
            
            % Build short label
            label_parts = {};
            for vi = 1:length(main_effects)
                me = main_effects{vi};
                info = main_effect_info(me);
                label_parts{end+1} = info.label; %#ok<AGROW>
            end
            for vi = 1:length(interactions)
                int_name = interactions{vi};
                if isKey(interaction_info, int_name)
                    info = interaction_info(int_name);
                    label_parts{end+1} = info.label; %#ok<AGROW>
                end
            end
            model_labels{mi} = strjoin(label_parts, '+');
            
            % Assign color based on model composition (explicit lookup)
            % Create a canonical key from sorted main effects
            sorted_main = sort(main_effects);
            main_key = strjoin(sorted_main, '+');
            has_int = ~isempty(interactions);
            
            % Pre-defined distinct colors for each combination
            % Single variables - pure colors
            combo_colors = containers.Map('KeyType', 'char', 'ValueType', 'any');
            combo_colors('Speed') = color_speed;
            combo_colors('TF') = color_tf;
            combo_colors('SF') = color_sf;
            combo_colors('OR') = color_or;
            
            % 2-variable combos - blue variations (will be striped with constituent colors)
            combo_colors('Speed+TF') = [0.12 0.47 0.71];    % Medium blue
            combo_colors('SF+Speed') = [0.26 0.58 0.85];    % Light blue
            combo_colors('OR+Speed') = [0.05 0.30 0.55];    % Dark blue
            combo_colors('SF+TF') = [0.40 0.60 0.80];       % Sky blue
            combo_colors('OR+TF') = [0.15 0.40 0.65];       % Steel blue
            combo_colors('OR+SF') = [0.30 0.50 0.75];       % Slate blue
            
            % 3-variable combos - lighter blue shades (will be striped with constituent colors)
            combo_colors('SF+Speed+TF') = [0.50 0.70 0.90]; % Light blue
            combo_colors('OR+Speed+TF') = [0.35 0.55 0.80]; % Medium-light blue
            combo_colors('OR+SF+Speed') = [0.45 0.65 0.85]; % Soft blue
            combo_colors('OR+SF+TF') = [0.55 0.75 0.92];    % Pale blue
            
            % 4-variable combo
            combo_colors('OR+SF+Speed+TF') = color_4plus;
            
            % Look up color
            if isKey(combo_colors, main_key)
                base_color = combo_colors(main_key);
            else
                % Fallback: use gray
                base_color = [0.5 0.5 0.5];
            end
            
            % If has interactions, make slightly darker and add edge
            if has_int
                model_colors(mi, :) = darken(base_color);
            else
                model_colors(mi, :) = base_color;
            end
        end
    end
    
    % Map complexity to bar position (0->1, 1->2, 2->3, 3->4, 4+->5)
    bar_positions = min(model_complexity_arr + 1, 5);
    
    % Build stacked bar data: rows = complexity levels (1-5), cols = model types
    max_subtypes = 15;  % Allow many subtypes per complexity
    bar_data = zeros(5, max_subtypes);
    bar_colors_mat = zeros(5, max_subtypes, 3);
    
    % Track which model types go where for legend and striping
    legend_info = struct('label', {}, 'color', {}, 'count', {}, 'main_effects', {});
    
    % Track bar segment positions for striping
    bar_segment_info = struct('clevel', {}, 'col_idx', {}, 'main_effects', {}, 'model_idx', {});
    
    % Assign counts to bar_data matrix by complexity level
    subtype_idx = ones(5, 1);  % Track column index for each complexity level
    for mi = 1:n_model_types
        clevel = bar_positions(mi);
        col_idx = subtype_idx(clevel);
        if col_idx <= max_subtypes && model_counts(mi) > 0
            bar_data(clevel, col_idx) = model_counts(mi);
            bar_colors_mat(clevel, col_idx, :) = model_colors(mi, :);
            subtype_idx(clevel) = subtype_idx(clevel) + 1;
            
            % Track segment info for striping
            bar_segment_info(end+1) = struct('clevel', clevel, 'col_idx', col_idx, ...
                'main_effects', {model_main_effects{mi}}, 'model_idx', mi);
            
            % Add to legend
            legend_info(end+1) = struct('label', model_labels{mi}, ...
                'color', model_colors(mi, :), 'count', model_counts(mi), ...
                'main_effects', {model_main_effects{mi}});
        end
    end
    
    % Create stacked bar plot
    complexity_labels = {'Null', '1 var', '2 vars', '3 vars', '4+ vars'};
    
    b_stack = bar(bar_data, 'stacked', 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5);
    
    % Apply colors to each stack segment
    for si = 1:size(bar_data, 2)
        if si <= length(b_stack)
            b_stack(si).FaceColor = 'flat';
            for ci_bar = 1:5
                b_stack(si).CData(ci_bar, :) = squeeze(bar_colors_mat(ci_bar, si, :))';
            end
        end
    end
    
    set(gca, 'XTick', 1:5, 'XTickLabel', complexity_labels, 'XTickLabelRotation', 30, 'FontSize', 8);
    ylabel('# Neurons');
    title('Model Complexity (Forward Selection)', 'FontSize', 10);
    set(gca, 'box', 'off');
    hold on;
    
    % Add diagonal stripes for multi-variable models (>=2 main effects)
    bar_width = 0.8;  % Default MATLAB bar width
    for seg_i = 1:length(bar_segment_info)
        seg = bar_segment_info(seg_i);
        me = seg.main_effects;
        
        if length(me) >= 2  % Only stripe multi-variable models
            % Get bar position and segment boundaries
            x_center = seg.clevel;
            x_left = x_center - bar_width/2;
            x_right = x_center + bar_width/2;
            
            % Calculate cumulative heights for this complexity level
            cum_heights = cumsum(bar_data(seg.clevel, :));
            if seg.col_idx == 1
                y_bottom = 0;
            else
                y_bottom = cum_heights(seg.col_idx - 1);
            end
            y_top = cum_heights(seg.col_idx);
            
            if y_top > y_bottom  % Only draw if segment has height
                % Get colors for each main effect
                me_colors = zeros(length(me), 3);
                for mei = 1:length(me)
                    me_colors(mei, :) = main_effect_info(me{mei}).color;
                end
                
                % Draw alternating horizontal bands within the bar segment
                n_effects = length(me);
                seg_height = y_top - y_bottom;
                band_height = seg_height / (n_effects * 3);  % Multiple thin bands
                
                y_pos = y_bottom;
                band_idx = 1;
                while y_pos < y_top
                    clr_idx = mod(band_idx - 1, n_effects) + 1;
                    band_top = min(y_pos + band_height, y_top);
                    
                    % Draw a colored band
                    patch([x_left, x_right, x_right, x_left], ...
                          [y_pos, y_pos, band_top, band_top], ...
                          me_colors(clr_idx, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
                    
                    y_pos = band_top;
                    band_idx = band_idx + 1;
                end
            end
        end
    end
    
    % Compute totals for labels
    complexity_totals = sum(bar_data, 2);
    
    % Add count and percentage labels on top
    for ki = 1:length(complexity_totals)
        if complexity_totals(ki) > 0
            text(ki, complexity_totals(ki) + max(complexity_totals)*0.03, ...
                sprintf('%d\n(%.0f%%)', complexity_totals(ki), 100*complexity_totals(ki)/n_total), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        end
    end
    ylim([0 max(complexity_totals)*1.25]);
    
    % Build legend with striped icons for multi-variable models
    legend_handles = [];
    legend_labels_final = {};
    for li = 1:length(legend_info)
        lbl = legend_info(li).label;
        cnt = legend_info(li).count;
        me = legend_info(li).main_effects;
        
        if length(me) <= 1
            % Single color square for null or single-variable models
            clr = legend_info(li).color;
            h = scatter(nan, nan, 50, clr, 's', 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3]);
        else
            % For multi-var models, use first main effect color with thick edge of second color
            clr1 = main_effect_info(me{1}).color;
            clr2 = main_effect_info(me{2}).color;
            h = scatter(nan, nan, 50, clr1, 's', 'filled', 'MarkerEdgeColor', clr2, 'LineWidth', 2);
        end
        legend_handles = [legend_handles; h]; %#ok<AGROW>
        legend_labels_final{end+1} = sprintf('%s (%d)', lbl, cnt); %#ok<AGROW>
    end
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels_final, 'Location', 'northeast', 'FontSize', 5);
    end
    hold off;
    
    % --- Panel 2: Variable Inclusion Rates (Main Effects) ---
    subplot(1, 3, 2);
    main_effect_counts = [sum(is_spprofiled); sum(is_tf); sum(is_sf); sum(is_or); sum(has_int)];
    main_effect_labels = {'Speed', 'TF', 'SF', 'OR', 'Any Interact.'};
    % Use consistent colors with Panel 1
    main_effect_colors = [color_speed; color_tf; color_sf; color_or; color_int];
    
    bh = barh(main_effect_counts, 'FaceColor', 'flat', 'EdgeColor', 'none');
    bh.CData = main_effect_colors;
    set(gca, 'YTick', 1:5, 'YTickLabel', main_effect_labels, 'FontSize', 9);
    xlabel('# Neurons');
    
    % Add count and percentage labels
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
    int_counts = [sum(has_spprofiled_tf); sum(has_spprofiled_sf); sum(has_spprofiled_or); ...
                  sum(has_tf_sf); sum(has_tf_or); sum(has_sf_or)];
    int_labels = {'Spd x TF', 'Spd x SF', 'Spd x OR', 'TF x SF', 'TF x OR', 'SF x OR'};
    int_colors = [0.9 0.2 0.2; 0.8 0.5 0.2; 0.6 0.2 0.6; 0.5 0.8 0.2; 0.2 0.5 0.8; 0.9 0.6 0.6];
    
    bh_int = barh(int_counts, 'FaceColor', 'flat', 'EdgeColor', 'none');
    bh_int.CData = int_colors;
    set(gca, 'YTick', 1:6, 'YTickLabel', int_labels, 'FontSize', 9);
    xlabel('# Neurons');
    
    % Add count and percentage labels
    max_int = max(max(int_counts), 1);  % Avoid division by zero
    for ki = 1:length(int_counts)
        text(int_counts(ki) + max_int*0.02, ki, ...
            sprintf('%d (%.0f%%)', int_counts(ki), 100*int_counts(ki)/n_total), ...
            'VerticalAlignment', 'middle', 'FontSize', 8);
    end
    xlim([0 max_int*1.4]);
    title('Interaction Breakdown', 'FontSize', 10);
    set(gca, 'box', 'off');
end

sgtitle(sprintf('Forward Selection Results (threshold: \\Delta bps > %.3f)', classification_threshold), ...
    'FontSize', 12, 'FontWeight', 'bold');

if save_figs
    print(fig2, fullfile(ctl.figs.curr_dir, 'model_selection_classification_summary.png'), '-dpng', '-r300');
    fprintf('  Saved model_selection_classification_summary.png\n');
end

% --- Print classification summary to console ---
fprintf('\n====================================================================\n');
fprintf('  NEURON CLASSIFICATION SUMMARY (Drop-One Variable Selection)\n');
fprintf('====================================================================\n');
fprintf('  Threshold: delta bps > %.3f\n\n', classification_threshold);

for gt = 1:1  % Only time-bin GLM now
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    n_total = n_unique_clusters;
    n_spd = sum(results.(gt_tag).is_spprofileeed_tuned);
    n_tf = sum(results.(gt_tag).is_tf_tuned);
    n_sf = sum(results.(gt_tag).is_sf_tuned);
    n_or = sum(results.(gt_tag).is_or_tuned);
    n_inter = sum(results.(gt_tag).has_significant_interaction);
    
    fprintf('  --- %s GLM ---\n', gt_label);
    fprintf('  Speed-tuned (unique):    %3d / %d (%.1f%%)\n', n_spd, n_total, 100*n_spd/n_total);
    fprintf('  TF-tuned (unique):       %3d / %d (%.1f%%)\n', n_tf, n_total, 100*n_tf/n_total);
    fprintf('  SF-tuned (unique):       %3d / %d (%.1f%%)\n', n_sf, n_total, 100*n_sf/n_total);
    fprintf('  OR-tuned (unique):       %3d / %d (%.1f%%)\n', n_or, n_total, 100*n_or/n_total);
    fprintf('  Significant interaction: %3d / %d (%.1f%%)\n\n', n_inter, n_total, 100*n_inter/n_total);
end

fprintf('====================================================================\n');

%% ====================================================================
%  Section 7b: Speed-Profile Cross-Validation
%  ====================================================================
%  Split data by the two speed profiles (first vs second half of the
%  motion cloud sequence) and perform 2-fold cross-validation.
%  Train the best selected GLM on one speed profile, test on the other.
%  Compare to standard 5-fold CV to assess generalization across profiles.
fprintf('\n--- Speed-Profile Cross-Validation (2-fold: profile 1 vs profile 2) ---\n');

% --- Determine speed profile for each trial ---
% The mc_sequence contains two repetitions of the stimulus set.
% First half of trial IDs = speed profile 1, second half = speed profile 2.
n_total_seq = length(mc_sequence);
sp_midpoint = floor(n_total_seq / 2);
fprintf('  Total trials in sequence: %d | Midpoint: %d\n', n_total_seq, sp_midpoint);
fprintf('  Profile 1: trial_id 1–%d | Profile 2: trial_id %d–%d\n', ...
    sp_midpoint, sp_midpoint + 1, n_total_seq);

% Pre-extract variables for parfor compatibility
sp_all_selected_vars = results.time.selected_vars;

% Output storage
sp_cv_out = cell(n_unique_clusters, 1);

fprintf('  Fitting speed-profile CV in parallel (%d clusters)...\n', n_unique_clusters);
parfor ci = 1:n_unique_clusters
    pid = unique_pids(ci);
    cid = unique_cids(ci);

    idx = T_master_time.probe_id == pid & T_master_time.cluster_id == cid;
    T_cl = T_master_time(idx, :);

    r = struct();
    r.valid = false;
    r.sel_spprofile = NaN;   r.null_spprofile = NaN;   r.delta_spprofile = NaN;
    r.sel_standard = NaN;  r.null_standard = NaN;   r.delta_standard = NaN;
    r.has_spprofileeed = false;
    r.speed_spprofileprofile = NaN;  r.speed_standard = NaN;
    r.n_profile1 = 0;         r.n_profile2 = 0;

    if height(T_cl) < 10
        sp_cv_out{ci} = r;
        continue;
    end

    sel_vars = sp_all_selected_vars{ci};

    % --- Speed-profile fold assignment (2 folds) ---
    sp_fold = ones(height(T_cl), 1);
    sp_fold(T_cl.trial_id > sp_midpoint) = 2;
    r.n_profile1 = sum(sp_fold == 1);
    r.n_profile2 = sum(sp_fold == 2);

    if r.n_profile1 < 5 || r.n_profile2 < 5
        sp_cv_out{ci} = r;
        continue;
    end

    % --- Standard CV folds (5-fold, condition-stratified) for comparison ---
    unique_trials = unique(T_cl(:, {'trial_id', 'condition'}), 'rows');
    n_tr = height(unique_trials);
    trial_fold = zeros(n_tr, 1);
    conds_list = unique(unique_trials.condition);
    for cond_i = 1:length(conds_list)
        c_idx = find(unique_trials.condition == conds_list(cond_i));
        n_c = length(c_idx);
        perm_idx = randperm(n_c);
        trial_fold(c_idx(perm_idx)) = mod((1:n_c)' - 1, n_cv_folds) + 1;
    end
    [~, ~, row_to_trial] = unique(T_cl(:, {'trial_id', 'condition'}), 'rows');
    std_fold = trial_fold(row_to_trial);

    % --- Prepare features ---
    speed_v = T_cl.speed;
    tf_v = T_cl.tf;
    sf_v = T_cl.sf; sf_v(isnan(sf_v)) = 0;
    or_v = T_cl.orientation; or_v(isnan(or_v)) = 0;
    tso_v = T_cl.time_since_onset;
    y = double(T_cl.spike_count);
    offset_v = log(time_bin_width) * ones(height(T_cl), 1);

    B_spd = make_raised_cosine_basis(speed_v, n_speed_bases, speed_range(1), speed_range(2));
    B_tf  = make_raised_cosine_basis(tf_v, n_tf_bases, tf_range(1), tf_range(2));
    B_ons = make_onset_kernel_basis(tso_v, n_onset_bases, onset_range(2));

    sf_ref = sort(unique(sf_v(sf_v ~= 0)));
    or_ref = sort(unique(or_v(or_v ~= 0)));

    % --- Null model (intercept + onset kernel) ---
    [X_null, ~] = assemble_design_matrix_selected(B_spd, B_tf, B_ons, ...
        sf_v, or_v, {}, sf_ref, or_ref);
    [~, r.null_spprofile]  = cross_validate_glm(X_null, y, offset_v, sp_fold,  0);
    [~, r.null_standard] = cross_validate_glm(X_null, y, offset_v, std_fold, 0);

    % --- Selected model ---
    [X_sel, ~] = assemble_design_matrix_selected(B_spd, B_tf, B_ons, ...
        sf_v, or_v, sel_vars, sf_ref, or_ref);
    if size(X_sel, 2) >= height(T_cl)
        sp_cv_out{ci} = r;
        continue;
    end
    [~, r.sel_spprofile]  = cross_validate_glm(X_sel, y, offset_v, sp_fold,  0);
    [~, r.sel_standard] = cross_validate_glm(X_sel, y, offset_v, std_fold, 0);

    r.delta_spprofile  = r.sel_spprofile  - r.null_spprofile;
    r.delta_standard = r.sel_standard - r.null_standard;

    % --- Speed contribution (selected model WITH vs WITHOUT speed) ---
    r.has_spprofileeed = ismember('Speed', sel_vars);
    if r.has_spprofileeed
        no_speed_vars = sel_vars(~contains(sel_vars, 'Speed'));
        [X_no_spd, ~] = assemble_design_matrix_selected(B_spd, B_tf, B_ons, ...
            sf_v, or_v, no_speed_vars, sf_ref, or_ref);
        if size(X_no_spd, 2) < height(T_cl)
            [~, sp_no_spd_bps]  = cross_validate_glm(X_no_spd, y, offset_v, sp_fold,  0);
            [~, std_no_spd_bps] = cross_validate_glm(X_no_spd, y, offset_v, std_fold, 0);
            r.speed_spprofileprofile  = r.sel_spprofile  - sp_no_spd_bps;
            r.speed_standard = r.sel_standard - std_no_spd_bps;
        end
    end

    r.valid = true;
    sp_cv_out{ci} = r;
end

% --- Collect results ---
sp_valid       = false(n_unique_clusters, 1);
delta_spprofile    = NaN(n_unique_clusters, 1);   % speed-profile CV delta (selected - null)
delta_standard   = NaN(n_unique_clusters, 1);   % standard CV delta (selected - null)
spd_spprofileprofile      = NaN(n_unique_clusters, 1);   % speed contribution, speed-profile CV
spd_standard     = NaN(n_unique_clusters, 1);   % speed contribution, standard CV
sp_has_spprofileeed   = false(n_unique_clusters, 1);

for ci = 1:n_unique_clusters
    r = sp_cv_out{ci};
    sp_valid(ci)     = r.valid;
    delta_spprofile(ci)  = r.delta_spprofile;
    delta_standard(ci) = r.delta_standard;
    sp_has_spprofileeed(ci) = r.has_spprofileeed;
    spd_spprofileprofile(ci)    = r.speed_spprofileprofile;
    spd_standard(ci)   = r.speed_standard;
end

n_sp_valid = sum(sp_valid);
n_sp_speed = sum(sp_valid & sp_has_spprofileeed);
fprintf('  Speed-profile CV complete: %d / %d clusters valid\n', n_sp_valid, n_unique_clusters);
fprintf('  Clusters with speed selected: %d / %d (%.0f%%)\n', ...
    n_sp_speed, n_sp_valid, 100*n_sp_speed/n_sp_valid);

% --- Summary statistics ---
% Helper: finite-only median and mean (guards against Inf from zero-spike folds)
fin = @(x) x(isfinite(x));

v  = sp_valid;
vs = sp_valid & sp_has_spprofileeed;

% Build arrays of finite values for each metric
d_standard  = fin(delta_standard(v));   % total delta, standard CV
d_spprofileprofile = fin(delta_spprofile(v)); % total delta, speed-profile CV
change = fin(delta_spprofile(v) - delta_standard(v));  % per-cluster absolute change

s_standard  = fin(spd_standard(vs));   % speed contribution, standard CV
s_spprofileprofile = fin(spd_spprofileprofile(vs)); % speed contribution, speed-profile CV

pct_retained_total = 100 * median(d_spprofile) / median(d_standard);
pct_lost_total     = 100 - pct_retained_total;

fprintf('\n====================================================================\n');
fprintf('  SPEED-PROFILE CROSS-VALIDATION SUMMARY\n');
fprintf('====================================================================\n');
fprintf('\n  How well does the selected GLM generalise across speed profiles?\n');
fprintf('  Standard 5-fold CV shuffles trials randomly (train/test share the same\n');
fprintf('  speed distribution). Speed-profile 2-fold CV trains on one full speed\n');
fprintf('  trajectory and tests on the other — a much stricter generalization test.\n');

% Per-cluster retained fraction: only clusters with positive standard CV delta
% (ratio is only meaningful when the model beats the null in standard CV).
pos_t = d_standard > 0;
pos_s = s_standard > 0;
ratio_t = d_spprofile(pos_t) ./ d_standard(pos_t);   % speed-profile CV / standard CV, per cluster
ratio_s = s_spprofile(pos_s) ./ s_standard(pos_s);
pct_retained_total = 100 * median(ratio_t);
pct_lost_total     = 100 - pct_retained_total;
pct_retained_spprofiled   = 100 * median(ratio_s);
pct_lost_spd       = 100 - pct_retained_spprofiled;

fprintf('\n  --- Total information gain: CV(Selected) − CV(Null)  [bits/spike] ---\n');
fprintf('  Each dot = one cluster. Value = how much the selected model beats\n');
fprintf('  the null (intercept + onset) in cross-validated bits/spike.\n');
fprintf('\n    Standard 5-fold CV:   median = %.4f  (n=%d finite)\n', median(d_standard), numel(d_standard));
fprintf('    Speed-profile 2-fold: median = %.4f  (n=%d finite)\n', median(d_spprofile),  numel(d_spprofile));
fprintf('\n    Per-cluster ratio (speed-profile / standard), n=%d clusters with std CV > 0:\n', sum(pos_t));
fprintf('    Median retained = %.1f%%  |  Median lost = %.1f%%\n', pct_retained_total, pct_lost_total);
fprintf('    (Ratio computed per cluster, then median taken across clusters)\n');
fprintf('\n    Per-cluster absolute change (speed-profile − standard):\n');
fprintf('    median = %.4f  (negative = model loses predictive power across profiles)\n', median(change));
if median(change) < 0
    fprintf('    -> Most clusters lose predictive power when tested on a new speed profile.\n');
else
    fprintf('    -> Model generalises as well or better across speed profiles.\n');
end

if any(vs)
    fprintf('\n  --- Speed unique contribution: CV(Selected) − CV(Selected without Speed)  [bits/spike] ---\n');
    fprintf('  Speed-tuned clusters only (n=%d). Value = how many bits/spike are\n', numel(s_standard));
    fprintf('  uniquely explained by the speed variable, under each CV scheme.\n');
    fprintf('\n    Standard 5-fold CV:   median = %.4f\n', median(s_standard));
    fprintf('    Speed-profile 2-fold: median = %.4f\n', median(s_spprofile));
    fprintf('\n    Per-cluster ratio (speed-profile / standard), n=%d clusters with std CV > 0:\n', sum(pos_s));
    fprintf('    Median retained = %.1f%%  |  Median lost = %.1f%%\n', pct_retained_spprofiled, pct_lost_spd);
    if pct_retained_spprofiled < 50
        fprintf('    -> Less than half of speed''s contribution survives across profiles.\n');
        fprintf('       Speed encoding may be partly driven by within-session trajectory\n');
        fprintf('       shape or speed history, not purely speed magnitude tuning.\n');
    elseif pct_retained_spprofiled < 80
        fprintf('    -> Moderate generalisation: speed magnitude is encoded but some\n');
        fprintf('       within-session structure also contributes.\n');
    else
        fprintf('    -> Strong generalisation: speed tuning is largely profile-independent.\n');
    end
end
fprintf('====================================================================\n');

% --- Figure: Speed-Profile CV — paired plots only, y-axis fixed to [-0.5, 0.5] ---
ylim_sp = [-0.5, 0.5];

% Build paired arrays: both values must be finite; outliers outside ylim are
% still drawn but clipped by the fixed axis (count reported in title).
pair_mask_t = isfinite(delta_standard) & isfinite(delta_spprofile) & v(:);
xp_t_standard = delta_standard(pair_mask_t);
xp_t_spprofile  = delta_spprofile(pair_mask_t);
n_pair_t = sum(pair_mask_t);
n_out_t  = sum(xp_t_standard < ylim_sp(1) | xp_t_standard > ylim_sp(2) | ...
               xp_t_spprofile  < ylim_sp(1) | xp_t_spprofile  > ylim_sp(2));

pair_mask_s = isfinite(spd_standard) & isfinite(spd_spprofileprofile) & vs(:);
xp_s_standard = spd_standard(pair_mask_s);
xp_s_spprofileprofile  = spd_spprofileprofile(pair_mask_s);
n_pair_s = sum(pair_mask_s);
n_out_s  = sum(xp_s_standard < ylim_sp(1) | xp_s_standard > ylim_sp(2) | ...
               xp_s_spprofileprofile  < ylim_sp(1) | xp_s_spprofileprofile  > ylim_sp(2));

% --- Statistical tests: is the decrease from standard CV to speed-profile CV significant? ---
%
% Wilcoxon signed-rank test (one-tailed, H1: standard CV > speed-profile CV).
% Chosen over paired t-test because bps differences are right-skewed with
% outliers — Wilcoxon tests the median shift and is robust to those extremes.
% One-tailed: the hypothesis is directional (speed-profile CV is strictly harder).

% Total delta bps
[p_wilcox_t, ~, stats_wilcox_t] = signrank(xp_t_standard, xp_t_spprofile, 'tail', 'right');

% Speed contribution (speed-tuned clusters only)
if any(pair_mask_s)
    [p_wilcox_s, ~, stats_wilcox_s] = signrank(xp_s_standard, xp_s_spprofileprofile, 'tail', 'right');
else
    p_wilcox_s = NaN;
end

% Helper: format p-value as string for figure annotation
fmt_p = @(p) ternary(p < 0.001, 'p < 0.001', sprintf('p = %.3f', p));

fprintf('\n  --- Statistical tests (Wilcoxon signed-rank, one-tailed: std CV > speed-profile CV) ---\n');
fprintf('  Total model delta bps (n=%d pairs):\n', n_pair_t);
fprintf('    W = %.0f,  %s\n', stats_wilcox_t.signedrank, fmt_p(p_wilcox_t));
if any(pair_mask_s)
    fprintf('  Speed contribution (n=%d speed-tuned pairs):\n', n_pair_s);
    fprintf('    W = %.0f,  %s\n', stats_wilcox_s.signedrank, fmt_p(p_wilcox_s));
end

% Colours
col_std_t = [0.55 0.55 0.55];
col_sp_t  = [0.25 0.45 0.75];
col_std_s = [0.55 0.55 0.55];
col_sp_s  = [0.15 0.58 0.15];

% Compute medians used across all three panels
med_t_standard = median(xp_t_standard);
med_t_spprofile  = median(xp_t_spprofile);
med_s_standard = median(xp_s_standard);
med_s_spprofileprofile  = median(xp_s_spprofileprofile);
% Per-cluster retained fraction for bar chart.
% Only clusters where standard CV delta > 0 (ratio only meaningful there).
% Clipped to [0, 100] so the stacked bar always sums to 100%.
pos_t_fig = xp_t_standard > 0;
pos_s_fig = xp_s_standard > 0;
pct_ret_t = min(100 * median(xp_t_spprofile(pos_t_fig) ./ xp_t_standard(pos_t_fig)), 100);
pct_ret_s = min(100 * median(xp_s_spprofileprofile(pos_s_fig) ./ xp_s_standard(pos_s_fig)), 100);

fig_sp = figure('Position', [50 100 1400 550], 'Name', 'Speed-Profile Cross-Validation');

% ---- Panel 1: Paired plot — total delta bps ----
ax1 = subplot(1, 3, 1);
hold(ax1, 'on');
for ki = 1:n_pair_t
    plot(ax1, [1,2], [xp_t_standard(ki), xp_t_spprofile(ki)], '-', ...
        'Color', [0.75 0.75 0.75 0.35], 'LineWidth', 0.7);
end
scatter(ax1, ones(n_pair_t,1),   xp_t_standard, 30, col_std_t, 'filled', 'MarkerFaceAlpha', 0.65);
scatter(ax1, 2*ones(n_pair_t,1), xp_t_spprofile,  30, col_sp_t,  'filled', 'MarkerFaceAlpha', 0.65);
plot(ax1, [1,2], [med_t_standard, med_t_spprofile], 'k-', 'LineWidth', 2.5);
plot(ax1, 1, med_t_standard, 'ko', 'MarkerFaceColor', col_std_t, 'MarkerSize', 9);
plot(ax1, 2, med_t_spprofile,  'ko', 'MarkerFaceColor', col_sp_t,  'MarkerSize', 9);
yline(ax1, 0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
xlim(ax1, [0.5 2.5]);
ylim(ax1, ylim_sp);
set(ax1, 'XTick', [1 2], 'XTickLabel', {'Standard CV', 'Speed-Profile CV'}, ...
    'box', 'off', 'TickDir', 'out', 'FontSize', 11);
ylabel(ax1, {'\Delta bps = CV(Selected model) − CV(Null model)', 'bits / spike'}, 'FontSize', 10);
title(ax1, sprintf('Total model information gain  (n=%d, %d clipped)', n_pair_t, n_out_t), 'FontSize', 11);
text(ax1, 1.5, ylim_sp(2)*0.88, ...
    sprintf('Wilcoxon: %s', fmt_p(p_wilcox_t)), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
hold(ax1, 'off');

% ---- Panel 2: Paired plot — speed unique contribution ----
if any(pair_mask_s)
    ax2 = subplot(1, 3, 2);
    hold(ax2, 'on');
    for ki = 1:n_pair_s
        plot(ax2, [1,2], [xp_s_standard(ki), xp_s_spprofileprofile(ki)], '-', ...
            'Color', [0.75 0.75 0.75 0.35], 'LineWidth', 0.7);
    end
    scatter(ax2, ones(n_pair_s,1),   xp_s_standard, 30, col_std_s, 'filled', 'MarkerFaceAlpha', 0.65);
    scatter(ax2, 2*ones(n_pair_s,1), xp_s_spprofileprofile,  30, col_sp_s,  'filled', 'MarkerFaceAlpha', 0.65);
    plot(ax2, [1,2], [med_s_standard, med_s_spprofileprofile], 'k-', 'LineWidth', 2.5);
    plot(ax2, 1, med_s_standard, 'ko', 'MarkerFaceColor', col_std_s, 'MarkerSize', 9);
    plot(ax2, 2, med_s_spprofileprofile,  'ko', 'MarkerFaceColor', col_sp_s,  'MarkerSize', 9);
    yline(ax2, 0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
    xlim(ax2, [0.5 2.5]);
    ylim(ax2, ylim_sp);
    set(ax2, 'XTick', [1 2], 'XTickLabel', {'Standard CV', 'Speed-Profile CV'}, ...
        'box', 'off', 'TickDir', 'out', 'FontSize', 11);
    ylabel(ax2, {'\Delta bps = CV(Selected) − CV(Selected without Speed)', 'bits / spike'}, 'FontSize', 10);
    title(ax2, sprintf('Speed unique contribution  (n=%d speed-tuned, %d clipped)', n_pair_s, n_out_s), 'FontSize', 11);
    text(ax2, 1.5, ylim_sp(2)*0.88, ...
        sprintf('Wilcoxon: %s', fmt_p(p_wilcox_s)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
    hold(ax2, 'off');
end

% ---- Panel 3: Stacked bar — retained vs profile-specific (% of standard CV) ----
% Each bar is normalised to 100% of the standard CV median.
% Bottom (coloured): the fraction that generalises across speed profiles.
% Top (light grey):  the fraction that is profile-specific (does not generalise).
ax3 = subplot(1, 3, 3);
hold(ax3, 'on');

bar_x      = [1, 2];
retained   = [pct_ret_t,       pct_ret_s];
not_gen    = [100 - pct_ret_t, 100 - pct_ret_s];

% Stacked bar: column 1 = retained %, column 2 = profile-specific %
bar_data = [retained', not_gen'];  % [n_bars × 2]
bh = bar(ax3, bar_x, bar_data, 0.55, 'stacked', 'EdgeColor', 'none');
bh(1).FaceColor = 'flat';
bh(1).CData(1,:) = col_sp_t;
bh(1).CData(2,:) = col_sp_s;
bh(2).FaceColor  = [0.88 0.88 0.88];

% Percentage labels inside each segment
for bi = 1:2
    ret_val = retained(bi);
    non_val = not_gen(bi);
    % Label in retained segment (centred at ret_val/2)
    text(ax3, bar_x(bi), ret_val / 2, sprintf('%.0f%%', ret_val), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w');
    % Label in non-generalising segment (centred in that segment)
    text(ax3, bar_x(bi), ret_val + non_val / 2, sprintf('%.0f%%', non_val), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.4 0.4 0.4]);
end

ylim(ax3, [0 100]);
set(ax3, 'XTick', bar_x, 'XTickLabel', {'Total model', 'Speed'}, ...
    'box', 'off', 'TickDir', 'out', 'FontSize', 11);
ylabel(ax3, {'Median per-cluster ratio: speed-profile CV / standard CV', '(% retained)'}, 'FontSize', 10);
title(ax3, 'How much generalises across speed profiles?', 'FontSize', 11);
legend(bh, {'Generalises', 'Profile-specific'}, 'Location', 'southeast', 'Box', 'off', 'FontSize', 10);
hold(ax3, 'off');

sgtitle('Speed-Profile Cross-Validation: Generalisation Across Speed Profiles', ...
    'FontSize', 11, 'FontWeight', 'bold');

if save_figs
    print(fig_sp, fullfile(ctl.figs.curr_dir, 'speed_profile_cross_validation.png'), '-dpng', '-r300');
    fprintf('  Saved speed_profile_cross_validation.png\n');
end

%% ====================================================================
%  Section 8: Per-cluster visualisations (Model Overview with Scatter)
%  ====================================================================
%  LAYOUT: 4 rows × 6 columns (Null, Selected, Additive, FullInteraction)
%    For each model row:
%      Col 1: Beta swarm plot (grouped by feature type)
%      Col 2: Mean scatter by condition with SEM error bars
%      Col 3: Scatter coloured by speed
%      Col 4: Scatter coloured by TF
%      Col 5: Scatter coloured by SF
%      Col 6: Scatter coloured by OR
%    Selected model row is highlighted with green border as "winner".
%  ====================================================================
fprintf('\n--- Generating per-cluster figures (Model Overview with Scatter) ---\n');

n_cols_scatter = 6;
n_rows_scatter = 4;  % Null, Selected, Additive, FullInteraction

% Models to visualize
vis_model_labels = {'Null', 'Selected', 'Additive', 'FullInteraction'};

% Condition colours
cond_map = containers.Map({'T_Vstatic', 'V', 'VT'}, ...
    {[0.2 0.7 0.2], [0.9 0.6 0.1], [0 0.4 0.8]});
cond_order = {'T_Vstatic', 'V', 'VT'};

% SF / OR colour maps - use highly distinct colors
sf_unique_all = [0.003, 0.006, 0.012];
sf_cmap = [0.2 0.6 0.2;    % 0.003 cpd = green
           0.2 0.4 0.9;    % 0.006 cpd = blue
           0.9 0.2 0.2];   % 0.012 cpd = red
or_unique_all = [-pi/4, 0, pi/4, pi/2];  % sorted order to match unique() output
% Use highly distinct colors for orientations: red, cyan, orange, purple
% Maximally distinct colors for accessibility
or_cmap = [0.9 0.1 0.1;    % -45° = bright red
           0.1 0.8 0.9;    % 0° = cyan
           0.95 0.55 0.0;  % 45° = orange
           0.6 0.1 0.85];  % 90° = purple

% Abbreviated model labels for tight spacing
model_abbrev = containers.Map();
model_abbrev('Null')            = 'Null';
model_abbrev('Selected')        = 'Selected';
model_abbrev('Additive')        = 'Additive';
model_abbrev('FullInteraction') = 'FullInt';
model_abbrev('M0')              = 'M0';
model_abbrev('M0_Speed')        = 'M0+S';
model_abbrev('M0_Speed_TF')     = 'M0+S+TF';
model_abbrev('M0_Speed_TF_SF')  = 'M0+S+TF+SF';

% Row definitions for scatter plots: {glm_type, model_label, display_prefix}
scatter_row_defs = cell(n_rows_scatter, 3);
for mi = 1:n_rows_scatter
    scatter_row_defs{mi, 1} = 'time';
    scatter_row_defs{mi, 2} = vis_model_labels{mi};
    scatter_row_defs{mi, 3} = sprintf('Time: %s', model_abbrev(vis_model_labels{mi}));
end

for ci = 1:n_unique_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);
    
    wm_time = results.time.winning_model{ci};
    selected_vars_str = results.time.selected_vars_str{ci};
    if isempty(selected_vars_str)
        selected_vars_str = '(none)';
    end
    
    % --- Build classification label string based on forward selection ---
    class_time_str = '';
    if results.time.is_spprofileeed_tuned(ci)
        class_time_str = [class_time_str, 'Spd'];  %#ok<AGROW>
    end
    if results.time.is_tf_tuned(ci)
        if ~isempty(class_time_str), class_time_str = [class_time_str, '+'];  end  %#ok<AGROW>
        class_time_str = [class_time_str, 'TF'];  %#ok<AGROW>
    end
    if results.time.is_sf_tuned(ci)
        if ~isempty(class_time_str), class_time_str = [class_time_str, '+'];  end  %#ok<AGROW>
        class_time_str = [class_time_str, 'SF'];  %#ok<AGROW>
    end
    if results.time.is_or_tuned(ci)
        if ~isempty(class_time_str), class_time_str = [class_time_str, '+'];  end  %#ok<AGROW>
        class_time_str = [class_time_str, 'OR'];  %#ok<AGROW>
    end
    if isempty(class_time_str)
        class_time_str = 'None';
    end
    if results.time.has_interaction(ci)
        class_time_str = [class_time_str, ' +Int'];  %#ok<AGROW>
    end
    
    % Figure for n_rows_scatter rows × 6 columns
    fig_scatter = figure('Position', [10 10 240*n_cols_scatter 120*n_rows_scatter], ...
        'Visible', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A3');
    
    sgtitle(sprintf('Probe: %s  Cluster: %d  |  Selected: %s  |  Class: %s', ...
        char(pid), cid, selected_vars_str, class_time_str), 'FontSize', 10, 'FontWeight', 'bold');
    
    % --- Pre-compute smoothed observed FR for this cluster (used in all scatter plots) ---
    % Extract cluster data once (same for all model rows)
    idx_cl_scatter = T_master_time.probe_id == pid & T_master_time.cluster_id == cid;
    T_cluster_scatter = T_master_time(idx_cl_scatter, :);
    
    % Smoothing kernel (boxcar, matching trial-level plots)
    n_smooth_bins_scatter = round(obs_fr_smooth_width / time_bin_width);
    if n_smooth_bins_scatter < 1, n_smooth_bins_scatter = 1; end
    boxcar_kernel_scatter = ones(n_smooth_bins_scatter, 1) / n_smooth_bins_scatter;
    
    % Compute smoothed observed FR per-trial (preserves temporal structure before binning)
    obs_fr_smoothed_scatter = nan(height(T_cluster_scatter), 1);
    unique_trials_scatter = unique(T_cluster_scatter.trial_id);
    for ti_s = 1:length(unique_trials_scatter)
        tr_id = unique_trials_scatter(ti_s);
        tr_mask_s = T_cluster_scatter.trial_id == tr_id;
        
        % Get data for this trial and sort by time
        T_trial_s = T_cluster_scatter(tr_mask_s, :);
        [~, sort_idx_s] = sort(T_trial_s.time_in_trial);
        spike_counts_s = T_trial_s.spike_count(sort_idx_s);
        
        % Smooth and convert to rate
        obs_fr_raw_s = spike_counts_s / time_bin_width;
        obs_fr_smooth_s = conv(obs_fr_raw_s, boxcar_kernel_scatter, 'same');
        
        % Unsort back to original order and store
        obs_fr_unsort_s = nan(size(obs_fr_smooth_s));
        obs_fr_unsort_s(sort_idx_s) = obs_fr_smooth_s;
        obs_fr_smoothed_scatter(tr_mask_s) = obs_fr_unsort_s;
    end
    
    for ri = 1:n_rows_scatter
        gt_tag = scatter_row_defs{ri, 1};
        ml     = scatter_row_defs{ri, 2};
        row_label = scatter_row_defs{ri, 3};
        
        T_all = T_master_time;
        
        % Highlight the Selected model row with green border
        is_winner = strcmp(ml, 'Selected');
        
        idx = T_all.probe_id == pid & T_all.cluster_id == cid;
        T_cluster = T_all(idx, :);
        vt_idx_cl = T_cluster.condition == "VT";
        
        has_predictions = ~isempty(cluster_predictions.(gt_tag){ci}) && ...
            isfield(cluster_predictions.(gt_tag){ci}, ml);
        
        % Grid sizes for scatter plots
        n_spd_grid = length(spd_bin_centers);
        n_tf_grid  = length(tf_bin_centers);
        
        bin_time = time_bin_width * ones(height(T_cluster), 1);
        
        % Col 1: Beta swarm plot (grouped by feature type)
        ax1 = subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 1);
        if has_predictions
            b_all = cluster_predictions.(gt_tag){ci}.(ml).beta;
            se_all = cluster_predictions.(gt_tag){ci}.(ml).se;
            cn_all = cluster_predictions.(gt_tag){ci}.(ml).col_names;
            n_coef = length(b_all);
            
            [grp_mean, ~, grp_clr, grp_lbl, n_grps, grp_betas, grp_cn] = ...
                compute_grouped_params(b_all, se_all, cn_all);
            
            hold on;
            for gii = 1:n_grps
                bv = grp_betas{gii};
                cn_g = grp_cn{gii};
                n_pts = length(bv);
                % Jitter x positions for visibility
                if n_pts > 1
                    jitter = linspace(-0.25, 0.25, n_pts);
                else
                    jitter = 0;
                end
                xpos = gii + jitter;
                % Plot basis numbers as text labels
                for pi = 1:n_pts
                    basis_lbl = extract_basis_label(cn_g{pi});
                    if isempty(basis_lbl), basis_lbl = '-'; end
                    text(xpos(pi), bv(pi), basis_lbl, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 5, 'Color', grp_clr(gii,:), ...
                        'FontWeight', 'bold');
                end
                % Group mean line
                plot([gii-0.3, gii+0.3], [grp_mean(gii), grp_mean(gii)], '-', ...
                    'Color', grp_clr(gii,:), 'LineWidth', 1.2);
            end
            yline(0, 'k-', 'LineWidth', 0.3);
            hold off;
            xlim([0.3, n_grps + 0.7]);
            set(gca, 'XTick', 1:n_grps, 'XTickLabel', grp_lbl, ...
                'XTickLabelRotation', 45, 'FontSize', 5, ...
                'TickLabelInterpreter', 'tex');
            ylabel('\beta', 'FontSize', 7);
            title_str = sprintf('%s (%d coefs)', row_label, n_coef);
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            title_str = row_label;
        end
        if is_winner
            title(sprintf('\\bf\\color{red}%s', title_str), 'FontSize', 7);
            set(ax1, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 0.8, 'box', 'on');
        else
            title(title_str, 'FontSize', 7);
            set(ax1, 'box', 'off');
        end
        set(ax1, 'FontSize', 6);
        
        % Col 2: Mean scatter by condition with STD error bars (across trials)
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 2);
        if has_predictions
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            obs_fr = obs_fr_smoothed_scatter;  % Use pre-computed smoothed FR
            trial_ids_cl = T_cluster.trial_id;
            
            hold on;
            n_conds = length(cond_order);
            mean_pred_cond = nan(n_conds, 1);
            mean_obs_cond = nan(n_conds, 1);
            std_pred_cond = nan(n_conds, 1);
            std_obs_cond = nan(n_conds, 1);
            
            for coi = 1:n_conds
                c_mask = T_cluster.condition == cond_order{coi};
                if sum(c_mask) > 1
                    % Group by trial: compute mean FR per trial within this condition
                    trials_in_cond = unique(trial_ids_cl(c_mask));
                    n_trials_cond = length(trials_in_cond);
                    obs_per_trial_cond = nan(n_trials_cond, 1);
                    pred_per_trial_cond = nan(n_trials_cond, 1);
                    for ti_tr = 1:n_trials_cond
                        tr_mask = c_mask & (trial_ids_cl == trials_in_cond(ti_tr));
                        if sum(tr_mask) > 0
                            obs_per_trial_cond(ti_tr) = mean(obs_fr(tr_mask));
                            pred_per_trial_cond(ti_tr) = mean(cv_pred_fr(tr_mask));
                        end
                    end
                    valid_trials_cond = ~isnan(obs_per_trial_cond);
                    if sum(valid_trials_cond) >= 2
                        mean_pred_cond(coi) = mean(pred_per_trial_cond(valid_trials_cond));
                        mean_obs_cond(coi) = mean(obs_per_trial_cond(valid_trials_cond));
                        std_pred_cond(coi) = std(pred_per_trial_cond(valid_trials_cond));
                        std_obs_cond(coi) = std(obs_per_trial_cond(valid_trials_cond));
                    end
                end
            end
            
            % Plot error bars and markers for each condition
            h_cond = gobjects(n_conds, 1);  % Store handles for legend
            leg_labels_cond = {};
            for coi = 1:n_conds
                if ~isnan(mean_pred_cond(coi))
                    clr = cond_map(cond_order{coi});
                    h_cond(coi) = errorbar(mean_pred_cond(coi), mean_obs_cond(coi), ...
                        std_obs_cond(coi), std_obs_cond(coi), ...
                        std_pred_cond(coi), std_pred_cond(coi), ...
                        'o', 'Color', clr, 'MarkerFaceColor', clr, ...
                        'MarkerSize', 5, 'LineWidth', 0.8, 'CapSize', 3);
                    leg_labels_cond{end+1} = cond_order{coi};
                end
            end
            
            valid_cond = ~isnan(mean_obs_cond);
            if any(valid_cond)
                max_v = max([mean_obs_cond(valid_cond) + std_obs_cond(valid_cond); ...
                             mean_pred_cond(valid_cond) + std_pred_cond(valid_cond)]);
                if max_v > 0, plot([0 max_v], [0 max_v], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off'); end
                ss_res = sum((mean_obs_cond(valid_cond) - mean_pred_cond(valid_cond)).^2);
                ss_tot = sum((mean_obs_cond(valid_cond) - mean(mean_obs_cond(valid_cond))).^2);
                r2 = ternary(ss_tot > 0, 1 - ss_res/ss_tot, NaN);
            else
                r2 = NaN;
            end
            
            % Add legend (after diagonal line to exclude it)
            valid_h = h_cond(isgraphics(h_cond));
            if ~isempty(valid_h)
                lg = legend(valid_h, leg_labels_cond, 'Location', 'northeastoutside', 'FontSize', 4, 'Box', 'on');
                lg.ItemTokenSize = [5, 5];
            end
            hold off;
            
            xlabel('Pred'); ylabel('Obs');
            title(sprintf('Cond R^2=%.2f', r2), 'FontSize', 7);
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
        end
        if is_winner
            set(gca, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 0.8, 'box', 'on');
        else
            set(gca, 'box', 'off');
        end
        set(gca, 'FontSize', 6);
        
        % Col 3: Speed scatter with STD error bars (across trials)
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 3);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = obs_fr_smoothed_scatter(vt_idx_cl);  % Use pre-computed smoothed FR
            trial_ids_vt = T_vt.trial_id;
            unique_trials_vt = unique(trial_ids_vt);
            n_trials_vt = length(unique_trials_vt);
            
            % For each speed bin, compute mean FR per trial, then mean/std across trials
            obs_s = nan(n_spd_grid, 1); pred_s = nan(n_spd_grid, 1); sv = nan(n_spd_grid, 1);
            std_obs_s = nan(n_spd_grid, 1); std_pred_s = nan(n_spd_grid, 1);
            for si = 1:n_spd_grid
                b_idx = T_vt.speed >= spd_bin_edges(si) & T_vt.speed < spd_bin_edges(si+1);
                if sum(b_idx) > 1
                    % Group by trial: compute mean FR per trial in this speed bin
                    obs_per_trial = nan(n_trials_vt, 1);
                    pred_per_trial = nan(n_trials_vt, 1);
                    for ti_tr = 1:n_trials_vt
                        tr_mask = b_idx & (trial_ids_vt == unique_trials_vt(ti_tr));
                        if sum(tr_mask) > 0
                            obs_per_trial(ti_tr) = mean(obs_fr_vt(tr_mask));
                            pred_per_trial(ti_tr) = mean(cv_pred_vt(tr_mask));
                        end
                    end
                    valid_trials = ~isnan(obs_per_trial);
                    if sum(valid_trials) >= 2
                        obs_s(si) = mean(obs_per_trial(valid_trials));
                        pred_s(si) = mean(pred_per_trial(valid_trials));
                        std_obs_s(si) = std(obs_per_trial(valid_trials));
                        std_pred_s(si) = std(pred_per_trial(valid_trials));
                        sv(si) = spd_bin_centers(si);
                    end
                end
            end
            vs = ~isnan(obs_s);
            if any(vs)
                hold on;
                % Plot error bars with color based on speed
                cmap_spd = parula(n_spd_grid);
                for si = find(vs)'
                    [~, cmap_idx] = min(abs(spd_bin_centers - sv(si)));
                    errorbar(pred_s(si), obs_s(si), ...
                        std_obs_s(si), std_obs_s(si), ...
                        std_pred_s(si), std_pred_s(si), ...
                        'o', 'Color', cmap_spd(cmap_idx,:), ...
                        'MarkerFaceColor', cmap_spd(cmap_idx,:), ...
                        'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 2);
                end
                colormap(gca, parula); cb = colorbar; cb.FontSize = 5;
                caxis([min(spd_bin_centers) max(spd_bin_centers)]);
                mx = max([obs_s(vs) + std_obs_s(vs); pred_s(vs) + std_pred_s(vs)]);
                if mx > 0, plot([0 mx], [0 mx], 'k--', 'LineWidth', 0.5); end
                hold off;
                ss_r = sum((obs_s(vs) - pred_s(vs)).^2);
                ss_t = sum((obs_s(vs) - mean(obs_s(vs))).^2);
                r2_s = ternary(ss_t > 0, 1 - ss_r/ss_t, NaN);
                title(sprintf('Spd R^2=%.2f', r2_s), 'FontSize', 7);
            end
            xlabel('Pred'); ylabel('Obs');
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
        end
        if is_winner
            set(gca, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 0.8, 'box', 'on');
        else
            set(gca, 'box', 'off');
        end
        set(gca, 'FontSize', 6);
        
        % Col 4: TF scatter with STD error bars (across trials)
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 4);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = obs_fr_smoothed_scatter(vt_idx_cl);  % Use pre-computed smoothed FR
            trial_ids_vt_tf = T_vt.trial_id;
            unique_trials_vt_tf = unique(trial_ids_vt_tf);
            n_trials_vt_tf = length(unique_trials_vt_tf);
            
            % For each TF bin, compute mean FR per trial, then mean/std across trials
            obs_tf = nan(n_tf_grid, 1); pred_tf = nan(n_tf_grid, 1); tfv = nan(n_tf_grid, 1);
            std_obs_tf = nan(n_tf_grid, 1); std_pred_tf = nan(n_tf_grid, 1);
            for ti = 1:n_tf_grid
                b_idx = T_vt.tf >= tf_bin_edges(ti) & T_vt.tf < tf_bin_edges(ti+1);
                if sum(b_idx) > 1
                    % Group by trial: compute mean FR per trial in this TF bin
                    obs_per_trial_tf = nan(n_trials_vt_tf, 1);
                    pred_per_trial_tf = nan(n_trials_vt_tf, 1);
                    for ti_tr = 1:n_trials_vt_tf
                        tr_mask = b_idx & (trial_ids_vt_tf == unique_trials_vt_tf(ti_tr));
                        if sum(tr_mask) > 0
                            obs_per_trial_tf(ti_tr) = mean(obs_fr_vt(tr_mask));
                            pred_per_trial_tf(ti_tr) = mean(cv_pred_vt(tr_mask));
                        end
                    end
                    valid_trials_tf = ~isnan(obs_per_trial_tf);
                    if sum(valid_trials_tf) >= 2
                        obs_tf(ti) = mean(obs_per_trial_tf(valid_trials_tf));
                        pred_tf(ti) = mean(pred_per_trial_tf(valid_trials_tf));
                        std_obs_tf(ti) = std(obs_per_trial_tf(valid_trials_tf));
                        std_pred_tf(ti) = std(pred_per_trial_tf(valid_trials_tf));
                        tfv(ti) = tf_bin_centers(ti);
                    end
                end
            end
            vt = ~isnan(obs_tf);
            if any(vt)
                hold on;
                % Plot error bars with color based on TF
                cmap_tf = hot(n_tf_grid);
                for ti = find(vt)'
                    [~, cmap_idx] = min(abs(tf_bin_centers - tfv(ti)));
                    errorbar(pred_tf(ti), obs_tf(ti), ...
                        std_obs_tf(ti), std_obs_tf(ti), ...
                        std_pred_tf(ti), std_pred_tf(ti), ...
                        'o', 'Color', cmap_tf(cmap_idx,:), ...
                        'MarkerFaceColor', cmap_tf(cmap_idx,:), ...
                        'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 2);
                end
                colormap(gca, hot); cb = colorbar; cb.FontSize = 5;
                caxis([min(tf_bin_centers) max(tf_bin_centers)]);
                mx = max([obs_tf(vt) + std_obs_tf(vt); pred_tf(vt) + std_pred_tf(vt)]);
                if mx > 0, plot([0 mx], [0 mx], 'k--', 'LineWidth', 0.5); end
                hold off;
                ss_r = sum((obs_tf(vt) - pred_tf(vt)).^2);
                ss_t = sum((obs_tf(vt) - mean(obs_tf(vt))).^2);
                r2_tf = ternary(ss_t > 0, 1 - ss_r/ss_t, NaN);
                title(sprintf('TF R^2=%.2f', r2_tf), 'FontSize', 7);
            end
            xlabel('Pred'); ylabel('Obs');
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
        end
        if is_winner
            set(gca, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 0.8, 'box', 'on');
        else
            set(gca, 'box', 'off');
        end
        set(gca, 'FontSize', 6);
        
        % Col 5: SF scatter with STD error bars (across trials)
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 5);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = obs_fr_smoothed_scatter(vt_idx_cl);  % Use pre-computed smoothed FR
            sf_vt = T_vt.sf;
            trial_ids_vt_sf = T_vt.trial_id;
            unique_trials_vt_sf = unique(trial_ids_vt_sf);
            n_trials_vt_sf = length(unique_trials_vt_sf);
            sf_u = unique(sf_vt(~isnan(sf_vt)));
            
            if ~isempty(sf_u)
                obs_sf = nan(length(sf_u), 1);
                pred_sf = nan(length(sf_u), 1);
                std_obs_sf = nan(length(sf_u), 1);
                std_pred_sf = nan(length(sf_u), 1);
                for sfi = 1:length(sf_u)
                    b_idx = sf_vt == sf_u(sfi);
                    if sum(b_idx) > 1
                        % Group by trial: compute mean FR per trial for this SF level
                        obs_per_trial_sf = nan(n_trials_vt_sf, 1);
                        pred_per_trial_sf = nan(n_trials_vt_sf, 1);
                        for ti_tr = 1:n_trials_vt_sf
                            tr_mask = b_idx & (trial_ids_vt_sf == unique_trials_vt_sf(ti_tr));
                            if sum(tr_mask) > 0
                                obs_per_trial_sf(ti_tr) = mean(obs_fr_vt(tr_mask));
                                pred_per_trial_sf(ti_tr) = mean(cv_pred_vt(tr_mask));
                            end
                        end
                        valid_trials_sf = ~isnan(obs_per_trial_sf);
                        if sum(valid_trials_sf) >= 2
                            obs_sf(sfi) = mean(obs_per_trial_sf(valid_trials_sf));
                            pred_sf(sfi) = mean(pred_per_trial_sf(valid_trials_sf));
                            std_obs_sf(sfi) = std(obs_per_trial_sf(valid_trials_sf));
                            std_pred_sf(sfi) = std(pred_per_trial_sf(valid_trials_sf));
                        end
                    end
                end
                hold on;
                h_sf = gobjects(length(sf_u), 1);  % Store handles for legend
                leg_labels_sf = {};
                
                % DEBUG: Print SF color mapping for first cluster, first row
                if ci == 1 && ri == 1
                    fprintf('\n=== DEBUG: SF scatter color mapping (cluster %d) ===\n', cid);
                    fprintf('  sf_unique_all: [%s]\n', sprintf('%.4f ', sf_unique_all));
                    fprintf('  sf_u from data: [%s]\n', sprintf('%.4f ', sf_u'));
                end
                
                for sfi = 1:length(sf_u)
                    if ~isnan(obs_sf(sfi))
                        [~, sc_idx] = min(abs(sf_unique_all - sf_u(sfi)));
                        sf_color = sf_cmap(sc_idx,:);  % Get the color for this SF level
                        
                        % DEBUG: Print color index assignment
                        if ci == 1 && ri == 1
                            fprintf('    sf_u(%d)=%.4f -> sc_idx=%d, color=[%.2f %.2f %.2f]\n', ...
                                sfi, sf_u(sfi), sc_idx, sf_color(1), sf_color(2), sf_color(3));
                        end
                        
                        h_sf(sfi) = errorbar(pred_sf(sfi), obs_sf(sfi), ...
                            std_obs_sf(sfi), std_obs_sf(sfi), ...
                            std_pred_sf(sfi), std_pred_sf(sfi), ...
                            'o', 'MarkerSize', 5, 'LineWidth', 0.8, 'CapSize', 3);
                        % Explicitly set colors on the handle (more reliable across MATLAB versions)
                        set(h_sf(sfi), 'Color', sf_color, ...
                            'MarkerFaceColor', sf_color, ...
                            'MarkerEdgeColor', sf_color);
                        leg_labels_sf{end+1} = sprintf('%.3f', sf_u(sfi));
                    end
                end
                vsf = ~isnan(obs_sf);
                mx = max([obs_sf(vsf) + std_obs_sf(vsf); pred_sf(vsf) + std_pred_sf(vsf)]);
                if ~isempty(mx) && mx > 0, plot([0 mx], [0 mx], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off'); end
                % Add legend (after diagonal line to exclude it)
                valid_h_sf = h_sf(isgraphics(h_sf));
                if ~isempty(valid_h_sf)
                    lg_sf = legend(valid_h_sf, leg_labels_sf, 'Location', 'northeastoutside', 'FontSize', 4, 'Box', 'on');
                    lg_sf.ItemTokenSize = [5, 5];
                end
                hold off;
                if sum(vsf) > 1
                    ss_r = sum((obs_sf(vsf) - pred_sf(vsf)).^2);
                    ss_t = sum((obs_sf(vsf) - mean(obs_sf(vsf))).^2);
                    r2_sf = ternary(ss_t > 0, 1 - ss_r/ss_t, NaN);
                    title(sprintf('SF R^2=%.2f', r2_sf), 'FontSize', 7);
                end
            end
            xlabel('Pred'); ylabel('Obs');
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
        end
        if is_winner
            set(gca, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 0.8, 'box', 'on');
        else
            set(gca, 'box', 'off');
        end
        set(gca, 'FontSize', 6);
        
        % Col 6: OR scatter with STD error bars (across trials)
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 6);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = obs_fr_smoothed_scatter(vt_idx_cl);  % Use pre-computed smoothed FR
            or_vt = T_vt.orientation;
            trial_ids_vt_or = T_vt.trial_id;
            unique_trials_vt_or = unique(trial_ids_vt_or);
            n_trials_vt_or = length(unique_trials_vt_or);
            or_u = unique(or_vt(~isnan(or_vt)));
            
            if ~isempty(or_u)
                obs_or = nan(length(or_u), 1);
                pred_or = nan(length(or_u), 1);
                std_obs_or = nan(length(or_u), 1);
                std_pred_or = nan(length(or_u), 1);
                for ori = 1:length(or_u)
                    b_idx = or_vt == or_u(ori);
                    if sum(b_idx) > 1
                        % Group by trial: compute mean FR per trial for this orientation
                        obs_per_trial_or = nan(n_trials_vt_or, 1);
                        pred_per_trial_or = nan(n_trials_vt_or, 1);
                        for ti_tr = 1:n_trials_vt_or
                            tr_mask = b_idx & (trial_ids_vt_or == unique_trials_vt_or(ti_tr));
                            if sum(tr_mask) > 0
                                obs_per_trial_or(ti_tr) = mean(obs_fr_vt(tr_mask));
                                pred_per_trial_or(ti_tr) = mean(cv_pred_vt(tr_mask));
                            end
                        end
                        valid_trials_or = ~isnan(obs_per_trial_or);
                        if sum(valid_trials_or) >= 2
                            obs_or(ori) = mean(obs_per_trial_or(valid_trials_or));
                            pred_or(ori) = mean(pred_per_trial_or(valid_trials_or));
                            std_obs_or(ori) = std(obs_per_trial_or(valid_trials_or));
                            std_pred_or(ori) = std(pred_per_trial_or(valid_trials_or));
                        end
                    end
                end
                hold on;
                h_or = gobjects(length(or_u), 1);  % Store handles for legend
                leg_labels_or = {};
                
                % DEBUG: Print orientation color mapping for first cluster, first row
                if ci == 1 && ri == 1
                    fprintf('\n=== DEBUG: OR scatter color mapping (cluster %d) ===\n', cid);
                    fprintf('  or_unique_all: [%s]\n', sprintf('%.4f ', or_unique_all));
                    fprintf('  or_u from data: [%s]\n', sprintf('%.4f ', or_u'));
                end
                
                for ori = 1:length(or_u)
                    if ~isnan(obs_or(ori))
                        [~, oc_idx] = min(abs(or_unique_all - or_u(ori)));
                        or_color = or_cmap(oc_idx,:);  % Get the color for this orientation
                        
                        % DEBUG: Print color index assignment
                        if ci == 1 && ri == 1
                            fprintf('    or_u(%d)=%.4f -> oc_idx=%d, color=[%.2f %.2f %.2f]\n', ...
                                ori, or_u(ori), oc_idx, or_color(1), or_color(2), or_color(3));
                        end
                        
                        h_or(ori) = errorbar(pred_or(ori), obs_or(ori), ...
                            std_obs_or(ori), std_obs_or(ori), ...
                            std_pred_or(ori), std_pred_or(ori), ...
                            'o', 'MarkerSize', 5, 'LineWidth', 0.8, 'CapSize', 3);
                        % Explicitly set colors on the handle (more reliable across MATLAB versions)
                        set(h_or(ori), 'Color', or_color, ...
                            'MarkerFaceColor', or_color, ...
                            'MarkerEdgeColor', or_color);
                        % Convert radians to degrees for legend
                        or_deg = round(rad2deg(or_u(ori)));
                        leg_labels_or{end+1} = sprintf('%d°', or_deg);
                    end
                end
                vor = ~isnan(obs_or);
                mx = max([obs_or(vor) + std_obs_or(vor); pred_or(vor) + std_pred_or(vor)]);
                if ~isempty(mx) && mx > 0, plot([0 mx], [0 mx], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off'); end
                % Add legend (after diagonal line to exclude it)
                valid_h_or = h_or(isgraphics(h_or));
                if ~isempty(valid_h_or)
                    lg_or = legend(valid_h_or, leg_labels_or, 'Location', 'northeastoutside', 'FontSize', 4, 'Box', 'on');
                    lg_or.ItemTokenSize = [5, 5];
                end
                hold off;
                if sum(vor) > 1
                    ss_r = sum((obs_or(vor) - pred_or(vor)).^2);
                    ss_t = sum((obs_or(vor) - mean(obs_or(vor))).^2);
                    r2_or = ternary(ss_t > 0, 1 - ss_r/ss_t, NaN);
                    title(sprintf('OR R^2=%.2f', r2_or), 'FontSize', 7);
                end
            end
            xlabel('Pred'); ylabel('Obs');
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
        end
        if is_winner
            set(gca, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 0.8, 'box', 'on');
        else
            set(gca, 'box', 'off');
        end
        set(gca, 'FontSize', 6);
    end
    
    % Save figure
    if save_figs
        fname_cl = sprintf('cluster_%s_%d.png', strrep(char(pid), filesep, '_'), cid);
        drawnow limitrate;
        exportgraphics(fig_scatter, fullfile(ctl.figs.curr_dir, fname_cl), 'Resolution', 200);
    end
    close(fig_scatter);
    
    if mod(ci, 5) == 0
        fprintf('  Saved %d/%d cluster figures\n', ci, n_unique_clusters);
    end
end

fprintf('  All %d cluster figures saved (Model Overview with Scatter)\n', n_unique_clusters);

%% ====================================================================
%  Section 8b: Per-cluster GLM-reconstructed tuning curves
%  ====================================================================
%  For every cluster, use the fitted GLM coefficients to reconstruct
%  marginal tuning curves for each covariate (Speed, TF, SF, OR).
%
%  Method:
%    For a given covariate (e.g. Speed), evaluate GLM predictions at the
%    pre-calculated bin centers from the tuning tables, while holding all
%    OTHER covariates at their mean (continuous) or mode (categorical) values.
%    The predicted firing rate is:
%       FR(x) = exp( intercept + eta_sweep(x) + eta_others_at_mean )
%
%  Figure layout per cluster (5 rows x 4 columns):
%    Row 1 — Original data tuning curves (observed)
%    Row 2 — Null model predictions (intercept + onset only)
%    Row 3 — Selected model predictions (forward-selected variables)
%    Row 4 — Additive model predictions (all main effects)
%    Row 5 — FullInteraction model predictions
%    Col 1 — Speed tuning | Col 2 — TF tuning
%    Col 3 — SF tuning    | Col 4 — OR tuning
%
%  Each panel shows:
%    • Dots with errorbars at pre-calculated bin centers
%    • Colour coding by condition (T_Vstatic, V, VT)
%  ====================================================================
fprintf('\n--- Generating per-cluster reconstructed tuning curves ---\n');

% Use pre-calculated bin centers for plotting (from tuning tables)
% spd_bin_centers and tf_bin_centers were loaded in Section 3
n_spd_plot_bins = length(spd_bin_centers);
n_tf_plot_bins  = length(tf_bin_centers);

% Pre-compute basis expansions at the bin centers for predictions
B_speed_bincenters_rc = make_raised_cosine_basis(spd_bin_centers(:), n_speed_bases, speed_range(1), speed_range(2));
B_tf_bincenters_rc    = make_raised_cosine_basis(tf_bin_centers(:),  n_tf_bases,    tf_range(1),    tf_range(2));

% Condition colours (reuse from Section 8)
cond_colors_rc = containers.Map({'T_Vstatic', 'V', 'VT'}, ...
    {[0.2 0.7 0.2], [0.9 0.6 0.1], [0 0.4 0.8]});

% Row labels for the figure (5 rows: Observed + 4 models)
row_labels_base = {'Original Data', 'GLM: Null', 'GLM: Selected', 'GLM: Additive', 'GLM: FullInteraction'};

for ci = 1:n_unique_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);

    fig_tc = figure('Position', [30 30 1800 1400], 'Visible', 'off', ...
        'PaperOrientation', 'portrait', 'PaperType', 'A3');

    sgtitle(sprintf('GLM Tuning Curves — Probe: %s  Cluster: %d', ...
        char(pid), cid), 'FontSize', 14, 'FontWeight', 'bold');

    % --- Gather common data from speed-bin data (used for observed data plotting) ---
    T_all_spd_rc = T_master_spd;
    idx_spd_rc = T_all_spd_rc.probe_id == pid & T_all_spd_rc.cluster_id == cid;
    T_cl_spd_rc = T_all_spd_rc(idx_spd_rc, :);
    
    % DEBUG: Print data summary for first cluster
    if ci == 1
        fprintf('\n======== DEBUG: Plotting data summary for Cluster %d ========\n', cid);
        fprintf('  Total rows in T_cl_spd_rc: %d\n', height(T_cl_spd_rc));
        if height(T_cl_spd_rc) > 0
            fprintf('  Conditions present: %s\n', strjoin(unique(string(T_cl_spd_rc.condition)), ', '));
            fprintf('  Columns: %s\n', strjoin(T_cl_spd_rc.Properties.VariableNames, ', '));
            % Show first few rows
            fprintf('  Sample data (first 5 rows):\n');
            for row_i = 1:min(5, height(T_cl_spd_rc))
                fprintf('    Row %d: cond=%s, speed=%.2f, tf=%.2f, sf=%.4f, or=%.4f, spk=%.1f, time=%.3f\n', ...
                    row_i, char(T_cl_spd_rc.condition(row_i)), T_cl_spd_rc.speed(row_i), ...
                    T_cl_spd_rc.tf(row_i), T_cl_spd_rc.sf(row_i), T_cl_spd_rc.orientation(row_i), ...
                    T_cl_spd_rc.spike_count(row_i), T_cl_spd_rc.time_in_bin(row_i));
            end
        else
            fprintf('  WARNING: No data found for this cluster!\n');
        end
    end
    
    % Collect combined data across conditions for covariate ranges
    speed_data = T_cl_spd_rc.speed;
    tf_data    = T_cl_spd_rc.tf;
    sf_data    = T_cl_spd_rc.sf;  sf_data(isnan(sf_data)) = 0;
    or_data    = T_cl_spd_rc.orientation;  % keep NaN for T_Vstatic (don't replace with 0, since 0 is a valid orientation)
    bin_t_spd  = T_cl_spd_rc.time_in_bin;
    
    % Mean / mode of covariates (used as "held constant" values for predictions)
    mean_speed = mean(speed_data);
    mean_tf    = mean(tf_data);
    mode_sf    = mode(sf_data(sf_data ~= 0));
    if isempty(mode_sf) || isnan(mode_sf), mode_sf = 0.006; end
    mode_or    = mode(or_data(~isnan(or_data)));  % exclude NaN (T_Vstatic) when computing mode
    
    % Basis evaluations at the mean/mode values (single-row matrices)
    B_spd_mean = make_raised_cosine_basis(mean_speed, n_speed_bases, speed_range(1), speed_range(2));
    B_tf_mean  = make_raised_cosine_basis(mean_tf,    n_tf_bases,    tf_range(1),    tf_range(2));
    
    % Onset kernel: for tuning curve predictions, we evaluate at "steady state"
    % (late in trial, after onset transient has decayed). Use t = 1.5s post-onset.
    steady_state_time = 1.5;  % seconds after motion onset (in plateau region)
    B_onset_mean = make_onset_kernel_basis(steady_state_time, n_onset_bases, onset_range(2));
    
    % SF / OR unique levels (for categorical tuning curves)
    % Use GLOBAL reference levels for GLM predictions (must match training)
    % Note: SF=0 is included because training converts NaN (T_Vstatic) to 0
    sf_levels_rc = sort([0; 0.003; 0.006; 0.012]);  % global SF levels used during training (includes 0 for T_Vstatic)
    sf_levels_plot = [0.003; 0.006; 0.012];  % actual SF values for plotting (exclude 0)
    or_levels_rc = sort([-pi/4; 0; pi/4; pi/2]); % global OR levels used during training (4 orientations)
    n_sf_rc = length(sf_levels_rc);
    n_sf_plot = length(sf_levels_plot);  % 3 actual SF levels for predictions
    n_or_rc = length(or_levels_rc);
    
    % Track y-limits for each column to synchronize across rows
    ylim_col = cell(1, 4);  % {speed, TF, SF, OR}
    for col_idx = 1:4
        ylim_col{col_idx} = [Inf, -Inf];  % [min, max]
    end
    ax_handles = gobjects(5, 4);  % Store subplot handles for y-limit adjustment (5 rows x 4 cols)
    
    % ================================================================
    %  ROW 1: ORIGINAL DATA TUNING CURVES (from pre-computed tuning tables)
    % ================================================================
    
    % --- Find the probe_info entry for this cluster's probe ---
    probe_idx_rc = [];
    for pi_rc = 1:length(probe_info)
        if ~isempty(probe_info(pi_rc).pid) && strcmp(probe_info(pi_rc).pid, char(pid))
            probe_idx_rc = pi_rc;
            break;
        end
    end
    
    % ---- Column 1: Speed tuning (from speed tuning table) ----
    ax_handles(1, 1) = subplot(5, 4, 1);
    hold on;
    if ~isempty(probe_idx_rc) && isfield(probe_info(probe_idx_rc), 'D_spd') && ~isempty(probe_info(probe_idx_rc).D_spd)
        D_spd_rc = probe_info(probe_idx_rc).D_spd;
        for coi = 1:length(cond_order)
            cond_name = cond_order{coi};
            % Skip V condition - no locomotion speed in visual-only
            if strcmp(cond_name, 'V'), continue; end
            % Find condition in tuning table
            cond_idx_spd = find(strcmp(D_spd_rc.trial_groups, cond_name), 1);
            if isempty(cond_idx_spd) || isempty(D_spd_rc.tuning_curves{cond_idx_spd})
                continue;
            end
            tc_array_spd = D_spd_rc.tuning_curves{cond_idx_spd};
            tc_cluster_ids_spprofiled = arrayfun(@(x) double(x.cluster_id), tc_array_spd);
            tc_idx_spd = find(tc_cluster_ids_spprofiled == cid, 1);
            if isempty(tc_idx_spd), continue; end
            
            tc_spd = tc_array_spd(tc_idx_spd);
            % tuning matrix: n_bins x n_trials (firing rate in Hz)
            tuning_mat_spd = tc_spd.tuning;
            bc_spd = tc_spd.bin_centers(:);
            
            % Compute mean and STD across trials
            fr_mean = nanmean(tuning_mat_spd, 2);
            n_trials_spprofiled = sum(~isnan(tuning_mat_spd), 2);
            fr_std = nanstd(tuning_mat_spd, 0, 2);
            
            vm = ~isnan(fr_mean) & n_trials_spprofiled > 1;
            clr = cond_colors_rc(cond_name);
            errorbar(bc_spd(vm), fr_mean(vm), fr_std(vm), 'o', ...
                'Color', clr, 'MarkerFaceColor', clr, 'MarkerSize', 4, ...
                'LineWidth', 0.8, 'CapSize', 3);
            plot(bc_spd(vm), fr_mean(vm), '-', 'Color', clr, 'LineWidth', 1);
        end
    end
    % Add stationary firing rate (black dot at speed=0)
    idx_stat_rc = T_master_time.probe_id == pid & T_master_time.cluster_id == cid & ...
                  T_master_time.condition == "stationary";
    if sum(idx_stat_rc) > 0
        stat_spikes = T_master_time.spike_count(idx_stat_rc);
        fr_stat_mean = mean(stat_spikes) / time_bin_width;
        fr_stat_std  = std(stat_spikes / time_bin_width);
        errorbar(0, fr_stat_mean, fr_stat_std, 'ko', 'MarkerFaceColor', 'k', ...
            'MarkerSize', 6, 'LineWidth', 1.2, 'CapSize', 4);
    end
    hold off;
    xlabel('Speed (cm/s)'); ylabel('FR (Hz)');
    title('Tuning Table: Speed', 'FontSize', 9);
    set(gca, 'box', 'off', 'FontSize', 8);
    yl = ylim; ylim_col{1} = [min(ylim_col{1}(1), yl(1)), max(ylim_col{1}(2), yl(2))];
    
    % ---- Column 2: TF tuning (from TF tuning table) ----
    ax_handles(1, 2) = subplot(5, 4, 2);
    hold on;
    if ~isempty(probe_idx_rc) && isfield(probe_info(probe_idx_rc), 'D_tf') && ~isempty(probe_info(probe_idx_rc).D_tf)
        D_tf_rc = probe_info(probe_idx_rc).D_tf;
        for coi = 1:length(cond_order)
            cond_name = cond_order{coi};
            % Skip T_Vstatic for TF (no visual flow)
            if strcmp(cond_name, 'T_Vstatic'), continue; end
            
            % Find condition in TF tuning table
            cond_idx_tf = find(strcmp(D_tf_rc.trial_groups, cond_name), 1);
            if isempty(cond_idx_tf) || isempty(D_tf_rc.tuning_curves{cond_idx_tf})
                continue;
            end
            tc_array_tf = D_tf_rc.tuning_curves{cond_idx_tf};
            tc_cluster_ids_tf = arrayfun(@(x) double(x.cluster_id), tc_array_tf);
            tc_idx_tf = find(tc_cluster_ids_tf == cid, 1);
            if isempty(tc_idx_tf), continue; end
            
            tc_tf = tc_array_tf(tc_idx_tf);
            % tuning matrix: n_bins x n_trials (firing rate in Hz)
            tuning_mat_tf = tc_tf.tuning;
            bc_tf = tc_tf.bin_centers(:);
            
            % Compute mean and STD across trials
            fr_mean_tf = nanmean(tuning_mat_tf, 2);
            n_trials_tf = sum(~isnan(tuning_mat_tf), 2);
            fr_std_tf = nanstd(tuning_mat_tf, 0, 2);
            
            vm = ~isnan(fr_mean_tf) & n_trials_tf > 1;
            clr = cond_colors_rc(cond_name);
            errorbar(bc_tf(vm), fr_mean_tf(vm), fr_std_tf(vm), 'o', ...
                'Color', clr, 'MarkerFaceColor', clr, 'MarkerSize', 4, ...
                'LineWidth', 0.8, 'CapSize', 3);
            plot(bc_tf(vm), fr_mean_tf(vm), '-', 'Color', clr, 'LineWidth', 1);
        end
    end
    hold off;
    xlabel('TF (Hz)'); ylabel('FR (Hz)');
    title('Tuning Table: TF', 'FontSize', 9);
    set(gca, 'box', 'off', 'FontSize', 8);
    yl = ylim; ylim_col{2} = [min(ylim_col{2}(1), yl(1)), max(ylim_col{2}(2), yl(2))];
    
    % ---- Column 3: SF tuning (observed) ----
    ax_handles(1, 3) = subplot(5, 4, 3);
    hold on;
    
    % DEBUG: Print SF data diagnostics for first cluster only
    if ci == 1
        fprintf('\n=== DEBUG: SF data for cluster %d ===\n', cid);
        fprintf('  Expected sf_levels_rc: [%s]\n', sprintf('%.4f ', sf_levels_rc));
        fprintf('  Unique sf_data values: [%s]\n', sprintf('%.4f ', unique(sf_data(~isnan(sf_data)))));
        fprintf('  Total sf_data points: %d (non-NaN: %d, non-zero: %d)\n', ...
            length(sf_data), sum(~isnan(sf_data)), sum(sf_data ~= 0 & ~isnan(sf_data)));
    end
    
    for coi = 1:length(cond_order)
        cond_mask = T_cl_spd_rc.condition == cond_order{coi};
        if ~any(cond_mask) || strcmp(cond_order{coi}, 'T_Vstatic'), continue; end
        sf_c  = sf_data(cond_mask);
        spk_c = T_cl_spd_rc.spike_count(cond_mask);
        bt_c  = bin_t_spd(cond_mask);
        
        % DEBUG: Print per-condition SF diagnostics for first cluster
        if ci == 1
            fprintf('  Condition %s: %d rows, unique SF: [%s]\n', ...
                cond_order{coi}, sum(cond_mask), sprintf('%.4f ', unique(sf_c(~isnan(sf_c)))));
        end
        
        fr_obs_sf = nan(n_sf_rc, 1);
        fr_std_sf = nan(n_sf_rc, 1);
        for sfi = 1:n_sf_rc
            % Use tolerance for floating-point comparison (0.001 is safe since SF values are at least 0.003 apart)
            bm = abs(sf_c - sf_levels_rc(sfi)) < 0.001;
            
            % DEBUG: Print match count for first cluster
            if ci == 1
                fprintf('    SF level %.4f: %d matches\n', sf_levels_rc(sfi), sum(bm));
            end
            
            if sum(bm) > 1
                trial_fr = spk_c(bm) ./ bt_c(bm);
                fr_obs_sf(sfi) = mean(trial_fr);
                fr_std_sf(sfi) = std(trial_fr);
            end
        end
        vm = ~isnan(fr_obs_sf);
        clr = cond_colors_rc(cond_order{coi});
        errorbar(find(vm), fr_obs_sf(vm), fr_std_sf(vm), 'o', ...
            'Color', clr, 'MarkerFaceColor', clr, 'MarkerSize', 5, ...
            'LineWidth', 0.8, 'CapSize', 4);
        plot(find(vm), fr_obs_sf(vm), '-', 'Color', clr, 'LineWidth', 1);
    end
    hold off;
    set(gca, 'XTick', 1:n_sf_rc, 'XTickLabel', ...
        arrayfun(@(v) sprintf('%.3f', v), sf_levels_rc, 'UniformOutput', false));
    xlabel('SF (cpd)'); ylabel('FR (Hz)');
    title('Original Data: SF', 'FontSize', 9);
    set(gca, 'box', 'off', 'FontSize', 8);
    yl = ylim; ylim_col{3} = [min(ylim_col{3}(1), yl(1)), max(ylim_col{3}(2), yl(2))];
    title('Original Data: SF', 'FontSize', 9);
    set(gca, 'box', 'off', 'FontSize', 8);
    
    % ---- Column 4: Orientation tuning (observed) ----
    ax_handles(1, 4) = subplot(5, 4, 4);
    hold on;
    
    % DEBUG: Print orientation data diagnostics for first cluster only
    if ci == 1
        fprintf('\n=== DEBUG: Orientation data for cluster %d ===\n', cid);
        fprintf('  Expected or_levels_rc: [%s]\n', sprintf('%.4f ', or_levels_rc));
        fprintf('  Unique or_data values: [%s]\n', sprintf('%.4f ', unique(or_data(~isnan(or_data)))));
        fprintf('  Total or_data points: %d (non-NaN: %d)\n', length(or_data), sum(~isnan(or_data)));
    end
    
    for coi = 1:length(cond_order)
        cond_mask = T_cl_spd_rc.condition == cond_order{coi};
        if ~any(cond_mask) || strcmp(cond_order{coi}, 'T_Vstatic'), continue; end
        or_c  = or_data(cond_mask);
        spk_c = T_cl_spd_rc.spike_count(cond_mask);
        bt_c  = bin_t_spd(cond_mask);
        
        % DEBUG: Print per-condition diagnostics for first cluster
        if ci == 1
            fprintf('  Condition %s: %d rows, unique OR: [%s]\n', ...
                cond_order{coi}, sum(cond_mask), sprintf('%.4f ', unique(or_c(~isnan(or_c)))));
        end
        
        fr_obs_or = nan(n_or_rc, 1);
        fr_std_or = nan(n_or_rc, 1);
        for ori = 1:n_or_rc
            % Use tolerance for floating-point comparison (0.01 is safe since orientations are ~0.785 apart)
            bm = abs(or_c - or_levels_rc(ori)) < 0.01;
            
            % DEBUG: Print match count for first cluster
            if ci == 1
                fprintf('    OR level %.4f: %d matches\n', or_levels_rc(ori), sum(bm));
            end
            
            if sum(bm) >= 1
                trial_fr = spk_c(bm) ./ bt_c(bm);
                fr_obs_or(ori) = mean(trial_fr);
                fr_std_or(ori) = std(trial_fr);
            end
        end
        vm = ~isnan(fr_obs_or);
        clr = cond_colors_rc(cond_order{coi});
        errorbar(find(vm), fr_obs_or(vm), fr_std_or(vm), 'o', ...
            'Color', clr, 'MarkerFaceColor', clr, 'MarkerSize', 5, ...
            'LineWidth', 0.8, 'CapSize', 4);
        plot(find(vm), fr_obs_or(vm), '-', 'Color', clr, 'LineWidth', 1);
    end
    hold off;
    set(gca, 'XTick', 1:n_or_rc, 'XTickLabel', ...
        arrayfun(@(v) format_radians(v), or_levels_rc, 'UniformOutput', false));
    xlabel('Orientation (rad)'); ylabel('FR (Hz)');
    title('Original Data: OR', 'FontSize', 9);
    set(gca, 'box', 'off', 'FontSize', 8);
    yl = ylim; ylim_col{4} = [min(ylim_col{4}(1), yl(1)), max(ylim_col{4}(2), yl(2))];
    
    % ================================================================
    %  ROWS 2-5: GLM PREDICTIONS (Null, Selected, Additive, FullInteraction)
    % ================================================================
    
    % Get selected variables for this cluster (for labeling)
    sel_vars_rc = results.time.selected_vars{ci};
    sel_label = ternary(isempty(sel_vars_rc), 'none', strjoin(sel_vars_rc, '+'));
    
    % Define row mapping: {glm_type, model_label, row_number, display_name}
    % Row 2: Null, Row 3: Selected, Row 4: Additive, Row 5: FullInteraction
    pred_rows = {
        'time', 'Null',            2, 'Null (intercept+onset)';
        'time', 'Selected',        3, sprintf('Selected (%s)', sel_label);
        'time', 'Additive',        4, 'Additive';
        'time', 'FullInteraction', 5, 'FullInteraction'
    };
    
    for pr_i = 1:size(pred_rows, 1)
        gt_tag   = pred_rows{pr_i, 1};
        ml_name  = pred_rows{pr_i, 2};
        row_num  = pred_rows{pr_i, 3};
        disp_name = pred_rows{pr_i, 4};
        
        % Check if this model was fitted
        has_model = ~isempty(cluster_predictions.(gt_tag){ci}) && ...
                    isfield(cluster_predictions.(gt_tag){ci}, ml_name);
        if ~has_model
            % Leave empty subplots with a message
            for col = 1:4
                subplot(5, 4, (row_num-1)*4 + col);
                text(0.5, 0.5, 'Not fitted', 'HorizontalAlignment', 'center', ...
                    'FontSize', 10, 'Color', [0.5 0.5 0.5]);
                axis off;
            end
            continue;
        end
        
        % Flag for Selected model (needs special design matrix assembly)
        is_selected_model = strcmp(ml_name, 'Selected');
        
        % Get coefficients and column names
        beta_rc = cluster_predictions.(gt_tag){ci}.(ml_name).beta;
        col_names_rc = cluster_predictions.(gt_tag){ci}.(ml_name).col_names;
        
        % Prepare onset kernel for predictions (steady-state: t > onset transient)
        % We use the same steady_state_time defined earlier (~1.5s)
        B_onset_use = B_onset_mean;  % Single row, evaluated at steady state
        
        % Row title uses display name
        row_title_suffix = disp_name;
        
        % ---- Column 1: Speed tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 1) = subplot(5, 4, (row_num-1)*4 + 1);
        hold on;
        
        n_bins_spprofiled = n_spd_plot_bins;
        B_onset_rep_spd = repmat(B_onset_use, n_bins_spprofiled, 1);  % Onset kernel at steady state
        
        % T_Vstatic: Speed sweep with TF=0, SF=reference (first level), OR=reference
        B_tf_zero = make_raised_cosine_basis(0, n_tf_bases, tf_range(1), tf_range(2));
        B_tf_rep_T = repmat(B_tf_zero, n_bins_spprofiled, 1);
        sf_sweep_T = repmat(sf_levels_rc(1), n_bins_spprofiled, 1);  % reference SF
        or_sweep_T = repmat(or_levels_rc(1), n_bins_spprofiled, 1);  % reference OR
        if is_selected_model
            X_sweep_T = build_design_matrix_from_colnames(col_names_rc, B_speed_bincenters_rc, B_tf_rep_T, ...
                B_onset_rep_spd, sf_sweep_T, or_sweep_T, sf_levels_rc, or_levels_rc);
        else
            [X_sweep_T, ~] = assemble_design_matrix(B_speed_bincenters_rc, B_tf_rep_T, ...
                B_onset_rep_spd, sf_sweep_T, or_sweep_T, ml_name, sf_levels_rc, or_levels_rc);
        end
        if size(X_sweep_T, 2) == length(beta_rc)
            fr_pred_T = min(exp(X_sweep_T * beta_rc), max_predicted_fr);
            clr_T = cond_colors_rc('T_Vstatic');
            plot(spd_bin_centers, fr_pred_T, 'o-', 'Color', clr_T, ...
                'MarkerFaceColor', clr_T, 'MarkerSize', 4, 'LineWidth', 1.2);
        end
        
        % V: Skip - no locomotion speed in visual-only condition
        
        % VT: Speed sweep with mean TF, mode SF, mode OR (both stimuli on)
        B_tf_rep_VT = repmat(B_tf_mean, n_bins_spprofiled, 1);
        sf_sweep_VT = repmat(mode_sf, n_bins_spprofiled, 1);
        or_sweep_VT = repmat(mode_or, n_bins_spprofiled, 1);
        if is_selected_model
            X_sweep_VT = build_design_matrix_from_colnames(col_names_rc, B_speed_bincenters_rc, B_tf_rep_VT, ...
                B_onset_rep_spd, sf_sweep_VT, or_sweep_VT, sf_levels_rc, or_levels_rc);
        else
            [X_sweep_VT, ~] = assemble_design_matrix(B_speed_bincenters_rc, B_tf_rep_VT, ...
                B_onset_rep_spd, sf_sweep_VT, or_sweep_VT, ml_name, sf_levels_rc, or_levels_rc);
        end
        if size(X_sweep_VT, 2) == length(beta_rc)
            fr_pred_VT = min(exp(X_sweep_VT * beta_rc), max_predicted_fr);
            clr_VT = cond_colors_rc('VT');
            plot(spd_bin_centers, fr_pred_VT, 'o-', 'Color', clr_VT, ...
                'MarkerFaceColor', clr_VT, 'MarkerSize', 4, 'LineWidth', 1.2);
        end
        
        % Stationary: Speed=0, onset kernel at baseline (before motion onset → zeros)
        B_spd_zero_stat = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        B_tf_zero_stat = make_raised_cosine_basis(0, n_tf_bases, tf_range(1), tf_range(2));
        B_onset_stat = zeros(1, n_onset_bases);  % Before motion onset → all zeros
        if is_selected_model
            X_stat = build_design_matrix_from_colnames(col_names_rc, B_spd_zero_stat, B_tf_zero_stat, ...
                B_onset_stat, sf_levels_rc(1), or_levels_rc(1), sf_levels_rc, or_levels_rc);
        else
            [X_stat, ~] = assemble_design_matrix(B_spd_zero_stat, B_tf_zero_stat, ...
                B_onset_stat, sf_levels_rc(1), or_levels_rc(1), ml_name, sf_levels_rc, or_levels_rc);
        end
        if size(X_stat, 2) == length(beta_rc)
            fr_pred_stat = min(exp(X_stat * beta_rc), max_predicted_fr);
            plot(0, fr_pred_stat, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.2);
        end
        
        hold off;
        xlabel('Speed (cm/s)'); ylabel('FR (Hz)');
        title(sprintf('%s: Speed', row_title_suffix), 'FontSize', 9);
        set(gca, 'box', 'off', 'FontSize', 8);
        yl = ylim; ylim_col{1} = [min(ylim_col{1}(1), yl(1)), max(ylim_col{1}(2), yl(2))];
        
        % ---- Column 2: TF tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 2) = subplot(5, 4, (row_num-1)*4 + 2);
        hold on;
        
        n_bins_tf = n_tf_plot_bins;
        B_onset_rep_tf = repmat(B_onset_use, n_bins_tf, 1);  % Onset kernel at steady state
        
        % T_Vstatic: No TF tuning (TF=0 always) - skip
        
        % V: TF sweep with Speed=0 (visual only)
        B_spd_zero_tf = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        B_spd_rep_V_tf = repmat(B_spd_zero_tf, n_bins_tf, 1);
        sf_sweep_V_tf = repmat(mode_sf, n_bins_tf, 1);
        or_sweep_V_tf = repmat(mode_or, n_bins_tf, 1);
        if is_selected_model
            X_sweep_V_tf = build_design_matrix_from_colnames(col_names_rc, B_spd_rep_V_tf, B_tf_bincenters_rc, ...
                B_onset_rep_tf, sf_sweep_V_tf, or_sweep_V_tf, sf_levels_rc, or_levels_rc);
        else
            [X_sweep_V_tf, ~] = assemble_design_matrix(B_spd_rep_V_tf, B_tf_bincenters_rc, ...
                B_onset_rep_tf, sf_sweep_V_tf, or_sweep_V_tf, ml_name, sf_levels_rc, or_levels_rc);
        end
        if size(X_sweep_V_tf, 2) == length(beta_rc)
            fr_pred_V_tf = min(exp(X_sweep_V_tf * beta_rc), max_predicted_fr);
            clr_V = cond_colors_rc('V');
            plot(tf_bin_centers, fr_pred_V_tf, 'o-', 'Color', clr_V, ...
                'MarkerFaceColor', clr_V, 'MarkerSize', 4, 'LineWidth', 1.2);
        end
        
        % VT: TF sweep with mean Speed (both stimuli on)
        B_spd_rep_VT_tf = repmat(B_spd_mean, n_bins_tf, 1);
        sf_sweep_VT_tf = repmat(mode_sf, n_bins_tf, 1);
        or_sweep_VT_tf = repmat(mode_or, n_bins_tf, 1);
        if is_selected_model
            X_sweep_VT_tf = build_design_matrix_from_colnames(col_names_rc, B_spd_rep_VT_tf, B_tf_bincenters_rc, ...
                B_onset_rep_tf, sf_sweep_VT_tf, or_sweep_VT_tf, sf_levels_rc, or_levels_rc);
        else
            [X_sweep_VT_tf, ~] = assemble_design_matrix(B_spd_rep_VT_tf, B_tf_bincenters_rc, ...
                B_onset_rep_tf, sf_sweep_VT_tf, or_sweep_VT_tf, ml_name, sf_levels_rc, or_levels_rc);
        end
        if size(X_sweep_VT_tf, 2) == length(beta_rc)
            fr_pred_VT_tf = min(exp(X_sweep_VT_tf * beta_rc), max_predicted_fr);
            clr_VT = cond_colors_rc('VT');
            plot(tf_bin_centers, fr_pred_VT_tf, 'o-', 'Color', clr_VT, ...
                'MarkerFaceColor', clr_VT, 'MarkerSize', 4, 'LineWidth', 1.2);
        end
        
        hold off;
        xlabel('TF (Hz)'); ylabel('FR (Hz)');
        title(sprintf('%s: TF', row_title_suffix), 'FontSize', 9);
        set(gca, 'box', 'off', 'FontSize', 8);
        yl = ylim; ylim_col{2} = [min(ylim_col{2}(1), yl(1)), max(ylim_col{2}(2), yl(2))];
        
        % ---- Column 3: SF tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 3) = subplot(5, 4, (row_num-1)*4 + 3);
        hold on;
        
        B_onset_one_sf = B_onset_use;  % Onset kernel at steady state (single row)
        
        % T_Vstatic: No SF tuning (no visual stimulus) - skip
        
        % V: SF sweep with Speed=0 (visual only) - use only actual SF values (exclude 0)
        B_spd_zero_sf = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        fr_pred_V_sf = nan(n_sf_plot, 1);
        for sfi = 1:n_sf_plot
            if is_selected_model
                X_one = build_design_matrix_from_colnames(col_names_rc, B_spd_zero_sf, B_tf_mean, ...
                    B_onset_one_sf, sf_levels_plot(sfi), mode_or, sf_levels_rc, or_levels_rc);
            else
                [X_one, ~] = assemble_design_matrix(B_spd_zero_sf, B_tf_mean, ...
                    B_onset_one_sf, sf_levels_plot(sfi), mode_or, ml_name, sf_levels_rc, or_levels_rc);
            end
            if size(X_one, 2) == length(beta_rc)
                fr_pred_V_sf(sfi) = min(exp(X_one * beta_rc), max_predicted_fr);
            end
        end
        vm = ~isnan(fr_pred_V_sf);
        if any(vm)
            clr_V = cond_colors_rc('V');
            plot(find(vm), fr_pred_V_sf(vm), 'o-', 'Color', clr_V, ...
                'MarkerFaceColor', clr_V, 'MarkerSize', 5, 'LineWidth', 1.2);
        end
        
        % VT: SF sweep with mean Speed (both stimuli on) - use only actual SF values (exclude 0)
        fr_pred_VT_sf = nan(n_sf_plot, 1);
        for sfi = 1:n_sf_plot
            if is_selected_model
                X_one = build_design_matrix_from_colnames(col_names_rc, B_spd_mean, B_tf_mean, ...
                    B_onset_one_sf, sf_levels_plot(sfi), mode_or, sf_levels_rc, or_levels_rc);
            else
                [X_one, ~] = assemble_design_matrix(B_spd_mean, B_tf_mean, ...
                    B_onset_one_sf, sf_levels_plot(sfi), mode_or, ml_name, sf_levels_rc, or_levels_rc);
            end
            if size(X_one, 2) == length(beta_rc)
                fr_pred_VT_sf(sfi) = min(exp(X_one * beta_rc), max_predicted_fr);
            end
        end
        vm = ~isnan(fr_pred_VT_sf);
        if any(vm)
            clr_VT = cond_colors_rc('VT');
            plot(find(vm), fr_pred_VT_sf(vm), 'o-', 'Color', clr_VT, ...
                'MarkerFaceColor', clr_VT, 'MarkerSize', 5, 'LineWidth', 1.2);
        end
        
        hold off;
        set(gca, 'XTick', 1:n_sf_plot, 'XTickLabel', ...
            arrayfun(@(v) sprintf('%.3f', v), sf_levels_plot, 'UniformOutput', false));
        xlabel('SF (cpd)'); ylabel('FR (Hz)');
        title(sprintf('%s: SF', row_title_suffix), 'FontSize', 9);
        set(gca, 'box', 'off', 'FontSize', 8);
        yl = ylim; ylim_col{3} = [min(ylim_col{3}(1), yl(1)), max(ylim_col{3}(2), yl(2))];
        
        % ---- Column 4: Orientation tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 4) = subplot(5, 4, (row_num-1)*4 + 4);
        hold on;
        
        B_onset_one_or = B_onset_use;  % Onset kernel at steady state (single row)
        
        % T_Vstatic: No orientation tuning (no visual stimulus) - skip
        
        % V: OR sweep with Speed=0 (visual only)
        B_spd_zero_or = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        fr_pred_V_or = nan(n_or_rc, 1);
        for ori = 1:n_or_rc
            if is_selected_model
                X_one = build_design_matrix_from_colnames(col_names_rc, B_spd_zero_or, B_tf_mean, ...
                    B_onset_one_or, mode_sf, or_levels_rc(ori), sf_levels_rc, or_levels_rc);
            else
                [X_one, ~] = assemble_design_matrix(B_spd_zero_or, B_tf_mean, ...
                    B_onset_one_or, mode_sf, or_levels_rc(ori), ml_name, sf_levels_rc, or_levels_rc);
            end
            if size(X_one, 2) == length(beta_rc)
                fr_pred_V_or(ori) = min(exp(X_one * beta_rc), max_predicted_fr);
            end
        end
        vm = ~isnan(fr_pred_V_or);
        if any(vm)
            clr_V = cond_colors_rc('V');
            plot(find(vm), fr_pred_V_or(vm), 'o-', 'Color', clr_V, ...
                'MarkerFaceColor', clr_V, 'MarkerSize', 5, 'LineWidth', 1.2);
        end
        
        % VT: OR sweep with mean Speed (both stimuli on)
        fr_pred_VT_or = nan(n_or_rc, 1);
        for ori = 1:n_or_rc
            if is_selected_model
                X_one = build_design_matrix_from_colnames(col_names_rc, B_spd_mean, B_tf_mean, ...
                    B_onset_one_or, mode_sf, or_levels_rc(ori), sf_levels_rc, or_levels_rc);
            else
                [X_one, ~] = assemble_design_matrix(B_spd_mean, B_tf_mean, ...
                    B_onset_one_or, mode_sf, or_levels_rc(ori), ml_name, sf_levels_rc, or_levels_rc);
            end
            if size(X_one, 2) == length(beta_rc)
                fr_pred_VT_or(ori) = min(exp(X_one * beta_rc), max_predicted_fr);
            end
        end
        vm = ~isnan(fr_pred_VT_or);
        if any(vm)
            clr_VT = cond_colors_rc('VT');
            plot(find(vm), fr_pred_VT_or(vm), 'o-', 'Color', clr_VT, ...
                'MarkerFaceColor', clr_VT, 'MarkerSize', 5, 'LineWidth', 1.2);
        end
        
        hold off;
        set(gca, 'XTick', 1:n_or_rc, 'XTickLabel', ...
            arrayfun(@(v) format_radians(v), or_levels_rc, 'UniformOutput', false));
        xlabel('Orientation (rad)'); ylabel('FR (Hz)');
        title(sprintf('%s: OR', row_title_suffix), 'FontSize', 9);
        set(gca, 'box', 'off', 'FontSize', 8);
        yl = ylim; ylim_col{4} = [min(ylim_col{4}(1), yl(1)), max(ylim_col{4}(2), yl(2))];
    end  % pred_rows loop

    % --- Synchronize y-axis limits across all rows for each column ---
    for col_idx = 1:4
        % Ensure valid limits (handle Inf if no data plotted)
        if isinf(ylim_col{col_idx}(1)) || isinf(ylim_col{col_idx}(2))
            continue;  % skip if no valid data in this column
        end
        % Add small padding to upper limit; always start from 0
        yl_range = ylim_col{col_idx}(2);  % range from 0 to max
        yl_pad = yl_range * 0.05;
        yl_final = [0, ylim_col{col_idx}(2) + yl_pad];
        for row_idx = 1:5
            if isvalid(ax_handles(row_idx, col_idx)) && isgraphics(ax_handles(row_idx, col_idx))
                ylim(ax_handles(row_idx, col_idx), yl_final);
            end
        end
    end

    % --- Legend (shared) ---
    % Build a legend at the bottom using invisible axes
    ax_leg = axes(fig_tc, 'Position', [0.05 0.01 0.9 0.02], 'Visible', 'off');
    hold(ax_leg, 'on');
    h_leg = gobjects(6, 1);
    h_leg(1) = plot(ax_leg, NaN, NaN, 'o', 'Color', [0.2 0.7 0.2], 'MarkerFaceColor', [0.2 0.7 0.2], 'MarkerSize', 5);
    h_leg(2) = plot(ax_leg, NaN, NaN, 'o', 'Color', [0.9 0.6 0.1], 'MarkerFaceColor', [0.9 0.6 0.1], 'MarkerSize', 5);
    h_leg(3) = plot(ax_leg, NaN, NaN, 'o', 'Color', [0 0.4 0.8], 'MarkerFaceColor', [0 0.4 0.8], 'MarkerSize', 5);
    h_leg(4) = plot(ax_leg, NaN, NaN, 'o-', 'Color', [0.2 0.7 0.2], 'MarkerFaceColor', [0.2 0.7 0.2], 'MarkerSize', 4, 'LineWidth', 1.2);
    h_leg(5) = plot(ax_leg, NaN, NaN, 'o-', 'Color', [0.9 0.6 0.1], 'MarkerFaceColor', [0.9 0.6 0.1], 'MarkerSize', 4, 'LineWidth', 1.2);
    h_leg(6) = plot(ax_leg, NaN, NaN, 'o-', 'Color', [0 0.4 0.8], 'MarkerFaceColor', [0 0.4 0.8], 'MarkerSize', 4, 'LineWidth', 1.2);
    legend(h_leg, {'T\_Vstatic (obs)', 'V (obs)', 'VT (obs)', 'T\_Vstatic (pred)', 'V (pred)', 'VT (pred)'}, ...
        'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 8);
    hold(ax_leg, 'off');

    % Save
    if save_figs
        fname_tc = sprintf('tuning_curves_%s_%d.png', strrep(char(pid), filesep, '_'), cid);
        drawnow limitrate;
        exportgraphics(fig_tc, fullfile(ctl.figs.curr_dir, fname_tc), 'Resolution', 200);
    end
    close(fig_tc);

    if mod(ci, 5) == 0
        fprintf('  Saved %d/%d tuning curve figures\n', ci, n_unique_clusters);
    end
end

fprintf('  All %d tuning curve figures saved\n', n_unique_clusters);

%% ====================================================================
%  Section 8c: Tuning Curve Pearson Correlations per Model
%  ====================================================================
%  For each cluster, compute Pearson correlations between observed and 
%  predicted tuning curves for each condition (T_Vstatic, V, VT) and each
%  stimulus type (Speed, TF, SF, OR).
%
%  Models compared:
%    - Selected (forward-selection winning model)
%    - Additive (all main effects)
%    - FullInteraction (additive + all pairwise interactions)
%
%  Correlations are computed for:
%    - Speed tuning: T_Vstatic, VT (V has no locomotion)
%    - TF tuning: V, VT (T_Vstatic has no visual flow)
%    - SF tuning: V, VT (T_Vstatic has no visual stimulus)
%    - OR tuning: V, VT (T_Vstatic has no visual stimulus)
%
%  Output: Summary figure with violin plots of correlation distributions
%  ====================================================================
fprintf('\n--- Computing tuning curve Pearson correlations ---\n');

% Initialize storage for correlations
tc_corr_results = struct();
model_names_corr = {'Selected', 'Additive', 'FullInteraction'};
stimulus_types = {'Speed', 'TF', 'SF', 'OR'};
condition_names = {'T_Vstatic', 'V', 'VT'};

% Pre-allocate correlation matrices: n_clusters x (models x stimuli x conditions)
for mi = 1:length(model_names_corr)
    mn = model_names_corr{mi};
    for si = 1:length(stimulus_types)
        st = stimulus_types{si};
        for cdi = 1:length(condition_names)
            cn = condition_names{cdi};
            field_name = sprintf('%s_%s_%s', mn, st, cn);
            tc_corr_results.(field_name) = nan(n_unique_clusters, 1);
        end
    end
end

% Store prefilter category and selected model type for each cluster
tc_corr_results.category = cell(n_unique_clusters, 1);
tc_corr_results.probe_id = cell(n_unique_clusters, 1);
tc_corr_results.cluster_id = nan(n_unique_clusters, 1);
tc_corr_results.selected_vars_str = cell(n_unique_clusters, 1);  % Winner model type
tc_corr_results.model_category = cell(n_unique_clusters, 1);     % Simplified category for plotting

for ci = 1:n_unique_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);
    
    tc_corr_results.probe_id{ci} = char(pid);
    tc_corr_results.cluster_id(ci) = cid;
    
    % Get prefilter category
    pf_idx = find(strcmp(prefilter_results.probe_id, char(pid)) & ...
                  prefilter_results.cluster_id == cid, 1);
    if ~isempty(pf_idx)
        tc_corr_results.category{ci} = prefilter_results.category{pf_idx};
    else
        tc_corr_results.category{ci} = 'unknown';
    end
    
    % Get selected model type from forward selection results
    tc_corr_results.selected_vars_str{ci} = results.time.selected_vars_str{ci};
    
    % Categorize for plotting: extract main effects pattern
    sv = results.time.selected_vars_str{ci};
    has_spprofiled = results.time.is_spprofileeed_tuned(ci);
    has_tf = results.time.is_tf_tuned(ci);
    has_sf = results.time.is_sf_tuned(ci);
    has_or = results.time.is_or_tuned(ci);
    has_int = results.time.has_interaction(ci);
    
    % Create simplified category based on main effects
    if isempty(sv) || strcmp(sv, '')
        tc_corr_results.model_category{ci} = 'None';
    elseif has_spprofiled && ~has_tf && ~has_sf && ~has_or
        tc_corr_results.model_category{ci} = 'Speed';
    elseif ~has_spprofiled && has_tf && ~has_sf && ~has_or
        tc_corr_results.model_category{ci} = 'TF';
    elseif has_spprofiled && has_tf && ~has_sf && ~has_or
        if has_int
            tc_corr_results.model_category{ci} = 'Speed+TF+Int';
        else
            tc_corr_results.model_category{ci} = 'Speed+TF';
        end
    elseif ~has_spprofiled && (has_sf || has_or) && ~has_tf
        tc_corr_results.model_category{ci} = 'Visual(SF/OR)';
    elseif has_spprofiled && (has_sf || has_or)
        tc_corr_results.model_category{ci} = 'Speed+Visual';
    elseif has_tf && (has_sf || has_or)
        tc_corr_results.model_category{ci} = 'TF+Visual';
    else
        tc_corr_results.model_category{ci} = 'Other';
    end
    
    % --- Find probe_info entry for this cluster ---
    probe_idx_corr = [];
    for pi_corr = 1:length(probe_info)
        if ~isempty(probe_info(pi_corr).pid) && strcmp(probe_info(pi_corr).pid, char(pid))
            probe_idx_corr = pi_corr;
            break;
        end
    end
    
    % --- Get data for observed tuning curves ---
    T_all_corr = T_master_spd;
    idx_corr = T_all_corr.probe_id == pid & T_all_corr.cluster_id == cid;
    T_cl_corr = T_all_corr(idx_corr, :);
    
    % Covariates for this cluster
    speed_data_c = T_cl_corr.speed;
    tf_data_c = T_cl_corr.tf;
    sf_data_c = T_cl_corr.sf; sf_data_c(isnan(sf_data_c)) = 0;
    or_data_c = T_cl_corr.orientation;
    bin_t_c = T_cl_corr.time_in_bin;
    
    % Mean/mode values for held-constant covariates
    mean_speed_c = mean(speed_data_c);
    mean_tf_c = mean(tf_data_c);
    mode_sf_c = mode(sf_data_c(sf_data_c ~= 0));
    if isempty(mode_sf_c) || isnan(mode_sf_c), mode_sf_c = 0.006; end
    mode_or_c = mode(or_data_c(~isnan(or_data_c)));
    
    % Basis evaluations at mean/mode
    B_spd_mean_c = make_raised_cosine_basis(mean_speed_c, n_speed_bases, speed_range(1), speed_range(2));
    B_tf_mean_c = make_raised_cosine_basis(mean_tf_c, n_tf_bases, tf_range(1), tf_range(2));
    B_onset_mean_c = make_onset_kernel_basis(1.5, n_onset_bases, onset_range(2));
    
    % Reference levels
    sf_levels_c = sort([0; 0.003; 0.006; 0.012]);
    sf_levels_plot_c = [0.003; 0.006; 0.012];
    or_levels_c = sort([-pi/4; 0; pi/4; pi/2]);
    
    % Loop over models
    for mi = 1:length(model_names_corr)
        mn = model_names_corr{mi};
        
        % Check if model exists
        if isempty(cluster_predictions.time{ci}) || ~isfield(cluster_predictions.time{ci}, mn)
            continue;
        end
        
        beta_c = cluster_predictions.time{ci}.(mn).beta;
        col_names_c = cluster_predictions.time{ci}.(mn).col_names;
        is_selected = strcmp(mn, 'Selected');
        
        % ============ SPEED TUNING ============
        % T_Vstatic: Speed sweep with TF=0
        B_tf_zero_c = make_raised_cosine_basis(0, n_tf_bases, tf_range(1), tf_range(2));
        B_tf_rep_T_c = repmat(B_tf_zero_c, n_spd_plot_bins, 1);
        B_onset_rep_spd_c = repmat(B_onset_mean_c, n_spd_plot_bins, 1);
        sf_sweep_T_c = repmat(sf_levels_c(1), n_spd_plot_bins, 1);
        or_sweep_T_c = repmat(or_levels_c(1), n_spd_plot_bins, 1);
        
        if is_selected
            X_spd_T = build_design_matrix_from_colnames(col_names_c, B_speed_bincenters_rc, B_tf_rep_T_c, ...
                B_onset_rep_spd_c, sf_sweep_T_c, or_sweep_T_c, sf_levels_c, or_levels_c);
        else
            [X_spd_T, ~] = assemble_design_matrix(B_speed_bincenters_rc, B_tf_rep_T_c, ...
                B_onset_rep_spd_c, sf_sweep_T_c, or_sweep_T_c, mn, sf_levels_c, or_levels_c);
        end
        
        if size(X_spd_T, 2) == length(beta_c)
            fr_pred_spprofiled_T = min(exp(X_spd_T * beta_c), max_predicted_fr);
            
            % Get observed speed tuning for T_Vstatic from tuning table
            if ~isempty(probe_idx_corr) && isfield(probe_info(probe_idx_corr), 'D_spd') && ~isempty(probe_info(probe_idx_corr).D_spd)
                D_spd_c = probe_info(probe_idx_corr).D_spd;
                cond_idx_T = find(strcmp(D_spd_c.trial_groups, 'T_Vstatic'), 1);
                if ~isempty(cond_idx_T) && ~isempty(D_spd_c.tuning_curves{cond_idx_T})
                    tc_array_T = D_spd_c.tuning_curves{cond_idx_T};
                    tc_cids_T = arrayfun(@(x) double(x.cluster_id), tc_array_T);
                    tc_idx_T = find(tc_cids_T == cid, 1);
                    if ~isempty(tc_idx_T)
                        tc_T = tc_array_T(tc_idx_T);
                        fr_obs_spprofiled_T = nanmean(tc_T.tuning, 2);
                        % Compute correlation
                        vm = ~isnan(fr_obs_spprofiled_T) & ~isnan(fr_pred_spprofiled_T);
                        if sum(vm) >= 3
                            r = corr(fr_obs_spprofiled_T(vm), fr_pred_spprofiled_T(vm), 'type', 'Pearson');
                            tc_corr_results.(sprintf('%s_Speed_T_Vstatic', mn))(ci) = r;
                        end
                    end
                end
            end
        end
        
        % VT: Speed sweep with mean TF
        B_tf_rep_VT_c = repmat(B_tf_mean_c, n_spd_plot_bins, 1);
        sf_sweep_VT_c = repmat(mode_sf_c, n_spd_plot_bins, 1);
        or_sweep_VT_c = repmat(mode_or_c, n_spd_plot_bins, 1);
        
        if is_selected
            X_spd_VT = build_design_matrix_from_colnames(col_names_c, B_speed_bincenters_rc, B_tf_rep_VT_c, ...
                B_onset_rep_spd_c, sf_sweep_VT_c, or_sweep_VT_c, sf_levels_c, or_levels_c);
        else
            [X_spd_VT, ~] = assemble_design_matrix(B_speed_bincenters_rc, B_tf_rep_VT_c, ...
                B_onset_rep_spd_c, sf_sweep_VT_c, or_sweep_VT_c, mn, sf_levels_c, or_levels_c);
        end
        
        if size(X_spd_VT, 2) == length(beta_c)
            fr_pred_spprofiled_VT = min(exp(X_spd_VT * beta_c), max_predicted_fr);
            
            % Get observed speed tuning for VT from tuning table
            if ~isempty(probe_idx_corr) && isfield(probe_info(probe_idx_corr), 'D_spd') && ~isempty(probe_info(probe_idx_corr).D_spd)
                D_spd_c = probe_info(probe_idx_corr).D_spd;
                cond_idx_VT = find(strcmp(D_spd_c.trial_groups, 'VT'), 1);
                if ~isempty(cond_idx_VT) && ~isempty(D_spd_c.tuning_curves{cond_idx_VT})
                    tc_array_VT = D_spd_c.tuning_curves{cond_idx_VT};
                    tc_cids_VT = arrayfun(@(x) double(x.cluster_id), tc_array_VT);
                    tc_idx_VT = find(tc_cids_VT == cid, 1);
                    if ~isempty(tc_idx_VT)
                        tc_VT = tc_array_VT(tc_idx_VT);
                        fr_obs_spprofiled_VT = nanmean(tc_VT.tuning, 2);
                        vm = ~isnan(fr_obs_spprofiled_VT) & ~isnan(fr_pred_spprofiled_VT);
                        if sum(vm) >= 3
                            r = corr(fr_obs_spprofiled_VT(vm), fr_pred_spprofiled_VT(vm), 'type', 'Pearson');
                            tc_corr_results.(sprintf('%s_Speed_VT', mn))(ci) = r;
                        end
                    end
                end
            end
        end
        
        % ============ TF TUNING ============
        B_spd_zero_c = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        B_spd_rep_V_tf_c = repmat(B_spd_zero_c, n_tf_plot_bins, 1);
        B_onset_rep_tf_c = repmat(B_onset_mean_c, n_tf_plot_bins, 1);
        sf_sweep_V_tf_c = repmat(mode_sf_c, n_tf_plot_bins, 1);
        or_sweep_V_tf_c = repmat(mode_or_c, n_tf_plot_bins, 1);
        
        % V condition: TF sweep with Speed=0
        if is_selected
            X_tf_V = build_design_matrix_from_colnames(col_names_c, B_spd_rep_V_tf_c, B_tf_bincenters_rc, ...
                B_onset_rep_tf_c, sf_sweep_V_tf_c, or_sweep_V_tf_c, sf_levels_c, or_levels_c);
        else
            [X_tf_V, ~] = assemble_design_matrix(B_spd_rep_V_tf_c, B_tf_bincenters_rc, ...
                B_onset_rep_tf_c, sf_sweep_V_tf_c, or_sweep_V_tf_c, mn, sf_levels_c, or_levels_c);
        end
        
        if size(X_tf_V, 2) == length(beta_c)
            fr_pred_tf_V = min(exp(X_tf_V * beta_c), max_predicted_fr);
            
            % Get observed TF tuning for V from tuning table
            if ~isempty(probe_idx_corr) && isfield(probe_info(probe_idx_corr), 'D_tf') && ~isempty(probe_info(probe_idx_corr).D_tf)
                D_tf_c = probe_info(probe_idx_corr).D_tf;
                cond_idx_V_tf = find(strcmp(D_tf_c.trial_groups, 'V'), 1);
                if ~isempty(cond_idx_V_tf) && ~isempty(D_tf_c.tuning_curves{cond_idx_V_tf})
                    tc_array_V_tf = D_tf_c.tuning_curves{cond_idx_V_tf};
                    tc_cids_V_tf = arrayfun(@(x) double(x.cluster_id), tc_array_V_tf);
                    tc_idx_V_tf = find(tc_cids_V_tf == cid, 1);
                    if ~isempty(tc_idx_V_tf)
                        tc_V_tf = tc_array_V_tf(tc_idx_V_tf);
                        fr_obs_tf_V = nanmean(tc_V_tf.tuning, 2);
                        vm = ~isnan(fr_obs_tf_V) & ~isnan(fr_pred_tf_V);
                        if sum(vm) >= 3
                            r = corr(fr_obs_tf_V(vm), fr_pred_tf_V(vm), 'type', 'Pearson');
                            tc_corr_results.(sprintf('%s_TF_V', mn))(ci) = r;
                        end
                    end
                end
            end
        end
        
        % VT condition: TF sweep with mean Speed
        B_spd_rep_VT_tf_c = repmat(B_spd_mean_c, n_tf_plot_bins, 1);
        sf_sweep_VT_tf_c = repmat(mode_sf_c, n_tf_plot_bins, 1);
        or_sweep_VT_tf_c = repmat(mode_or_c, n_tf_plot_bins, 1);
        
        if is_selected
            X_tf_VT = build_design_matrix_from_colnames(col_names_c, B_spd_rep_VT_tf_c, B_tf_bincenters_rc, ...
                B_onset_rep_tf_c, sf_sweep_VT_tf_c, or_sweep_VT_tf_c, sf_levels_c, or_levels_c);
        else
            [X_tf_VT, ~] = assemble_design_matrix(B_spd_rep_VT_tf_c, B_tf_bincenters_rc, ...
                B_onset_rep_tf_c, sf_sweep_VT_tf_c, or_sweep_VT_tf_c, mn, sf_levels_c, or_levels_c);
        end
        
        if size(X_tf_VT, 2) == length(beta_c)
            fr_pred_tf_VT = min(exp(X_tf_VT * beta_c), max_predicted_fr);
            
            % Get observed TF tuning for VT
            if ~isempty(probe_idx_corr) && isfield(probe_info(probe_idx_corr), 'D_tf') && ~isempty(probe_info(probe_idx_corr).D_tf)
                D_tf_c = probe_info(probe_idx_corr).D_tf;
                cond_idx_VT_tf = find(strcmp(D_tf_c.trial_groups, 'VT'), 1);
                if ~isempty(cond_idx_VT_tf) && ~isempty(D_tf_c.tuning_curves{cond_idx_VT_tf})
                    tc_array_VT_tf = D_tf_c.tuning_curves{cond_idx_VT_tf};
                    tc_cids_VT_tf = arrayfun(@(x) double(x.cluster_id), tc_array_VT_tf);
                    tc_idx_VT_tf = find(tc_cids_VT_tf == cid, 1);
                    if ~isempty(tc_idx_VT_tf)
                        tc_VT_tf = tc_array_VT_tf(tc_idx_VT_tf);
                        fr_obs_tf_VT = nanmean(tc_VT_tf.tuning, 2);
                        vm = ~isnan(fr_obs_tf_VT) & ~isnan(fr_pred_tf_VT);
                        if sum(vm) >= 3
                            r = corr(fr_obs_tf_VT(vm), fr_pred_tf_VT(vm), 'type', 'Pearson');
                            tc_corr_results.(sprintf('%s_TF_VT', mn))(ci) = r;
                        end
                    end
                end
            end
        end
        
        % ============ SF TUNING ============
        % Compute SF tuning for V and VT conditions
        for cond_sf = {'V', 'VT'}
            cond_name_sf = cond_sf{1};
            if strcmp(cond_name_sf, 'V')
                B_spd_sf = B_spd_zero_c;
            else
                B_spd_sf = B_spd_mean_c;
            end
            
            fr_pred_sf = nan(length(sf_levels_plot_c), 1);
            for sfi = 1:length(sf_levels_plot_c)
                if is_selected
                    X_one = build_design_matrix_from_colnames(col_names_c, B_spd_sf, B_tf_mean_c, ...
                        B_onset_mean_c, sf_levels_plot_c(sfi), mode_or_c, sf_levels_c, or_levels_c);
                else
                    [X_one, ~] = assemble_design_matrix(B_spd_sf, B_tf_mean_c, ...
                        B_onset_mean_c, sf_levels_plot_c(sfi), mode_or_c, mn, sf_levels_c, or_levels_c);
                end
                if size(X_one, 2) == length(beta_c)
                    fr_pred_sf(sfi) = min(exp(X_one * beta_c), max_predicted_fr);
                end
            end
            
            % Get observed SF tuning from time-bin data
            cond_mask_sf = T_cl_corr.condition == cond_name_sf;
            if any(cond_mask_sf)
                sf_c_obs = sf_data_c(cond_mask_sf);
                spk_c_obs = T_cl_corr.spike_count(cond_mask_sf);
                bt_c_obs = bin_t_c(cond_mask_sf);
                
                fr_obs_sf = nan(length(sf_levels_plot_c), 1);
                for sfi = 1:length(sf_levels_plot_c)
                    bm = abs(sf_c_obs - sf_levels_plot_c(sfi)) < 0.001;
                    if sum(bm) >= 1
                        fr_obs_sf(sfi) = mean(spk_c_obs(bm) ./ bt_c_obs(bm));
                    end
                end
                
                vm = ~isnan(fr_obs_sf) & ~isnan(fr_pred_sf);
                if sum(vm) >= 3
                    r = corr(fr_obs_sf(vm), fr_pred_sf(vm), 'type', 'Pearson');
                    tc_corr_results.(sprintf('%s_SF_%s', mn, cond_name_sf))(ci) = r;
                end
            end
        end
        
        % ============ OR TUNING ============
        % Compute OR tuning for V and VT conditions
        for cond_or = {'V', 'VT'}
            cond_name_or = cond_or{1};
            if strcmp(cond_name_or, 'V')
                B_spd_or = B_spd_zero_c;
            else
                B_spd_or = B_spd_mean_c;
            end
            
            fr_pred_or = nan(length(or_levels_c), 1);
            for ori = 1:length(or_levels_c)
                if is_selected
                    X_one = build_design_matrix_from_colnames(col_names_c, B_spd_or, B_tf_mean_c, ...
                        B_onset_mean_c, mode_sf_c, or_levels_c(ori), sf_levels_c, or_levels_c);
                else
                    [X_one, ~] = assemble_design_matrix(B_spd_or, B_tf_mean_c, ...
                        B_onset_mean_c, mode_sf_c, or_levels_c(ori), mn, sf_levels_c, or_levels_c);
                end
                if size(X_one, 2) == length(beta_c)
                    fr_pred_or(ori) = min(exp(X_one * beta_c), max_predicted_fr);
                end
            end
            
            % Get observed OR tuning from time-bin data
            cond_mask_or = T_cl_corr.condition == cond_name_or;
            if any(cond_mask_or)
                or_c_obs = or_data_c(cond_mask_or);
                spk_c_obs = T_cl_corr.spike_count(cond_mask_or);
                bt_c_obs = bin_t_c(cond_mask_or);
                
                fr_obs_or = nan(length(or_levels_c), 1);
                for ori = 1:length(or_levels_c)
                    bm = abs(or_c_obs - or_levels_c(ori)) < 0.01;
                    if sum(bm) >= 1
                        fr_obs_or(ori) = mean(spk_c_obs(bm) ./ bt_c_obs(bm));
                    end
                end
                
                vm = ~isnan(fr_obs_or) & ~isnan(fr_pred_or);
                if sum(vm) >= 3
                    r = corr(fr_obs_or(vm), fr_pred_or(vm), 'type', 'Pearson');
                    tc_corr_results.(sprintf('%s_OR_%s', mn, cond_name_or))(ci) = r;
                end
            end
        end
    end
    
    if mod(ci, 20) == 0
        fprintf('  Computed correlations for %d/%d clusters\n', ci, n_unique_clusters);
    end
end

fprintf('  Completed correlation computation for %d clusters\n', n_unique_clusters);

%% ====================================================================
%  Section 8d: Tuning Curve Correlation Summary Figure
%  ====================================================================
%  Generate box plots showing the distribution of Pearson correlations
%  between OBSERVED and PREDICTED tuning curves (from the plots).
%  
%  Layout:
%    - Rows: Models (Selected, Additive, FullInteraction)
%    - Columns: Stimulus types (Speed, TF, SF, OR)
%    - Within each subplot: 2 box plots for the applicable conditions
%        Speed: T_Vstatic, VT (V has no locomotion)
%        TF/SF/OR: V, VT (T_Vstatic has no visual flow)
%  
%  Each correlation is computed within a single condition - e.g.,
%  observed Speed tuning in T vs predicted Speed tuning in T.
%  ====================================================================
fprintf('\n--- Generating tuning curve correlation summary figure ---\n');

% Colors for conditions (matching tuning curve plots)
cond_colors = containers.Map();
cond_colors('T_Vstatic') = [0.2 0.7 0.2];  % Green - T condition
cond_colors('V') = [0.9 0.6 0.1];          % Orange-yellow - V condition (VF: TF, SF, OR)
cond_colors('VT') = [0 0.4 0.8];           % Blue - VT condition

% Condition labels for x-axis
cond_labels = containers.Map();
cond_labels('T_Vstatic') = 'T';
cond_labels('V') = 'V';
cond_labels('VT') = 'VT';

% --- Model category markers (for winner model type visualization) ---
% Define unique markers for each model category from forward selection
model_cat_markers = containers.Map();
model_cat_markers('None') = 'o';           % Circle - no variables selected
model_cat_markers('Speed') = '^';          % Triangle up - Speed only
model_cat_markers('TF') = 'v';             % Triangle down - TF only  
model_cat_markers('Speed+TF') = 's';       % Square - Speed+TF without interaction
model_cat_markers('Speed+TF+Int') = 'd';   % Diamond - Speed+TF with interaction
model_cat_markers('Visual(SF/OR)') = 'p';  % Pentagon - SF/OR only
model_cat_markers('Speed+Visual') = 'h';   % Hexagram - Speed + SF/OR
model_cat_markers('TF+Visual') = '+';      % Plus - TF + SF/OR
model_cat_markers('Other') = 'x';          % Cross - other combinations

% Get unique categories actually present in data
all_model_cats = unique(tc_corr_results.model_category);
all_model_cats = all_model_cats(~cellfun(@isempty, all_model_cats));  % Remove empty

% Figure: 3 rows (models) x 4 columns (stimulus types)
fig_corr = figure('Position', [50 50 1400 900], 'Visible', 'off');
sgtitle('Tuning Curve Correlations: Observed vs Predicted', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Define layout
stimulus_types_plot = {'Speed', 'TF', 'SF', 'OR'};
model_names_plot = {'Selected', 'Additive', 'FullInteraction'};
model_labels = {'Selected (Winner)', 'Additive', 'FullInteraction'};

% For each stimulus, which conditions are valid:
% Speed: T_Vstatic, VT (V has no locomotion)
% TF, SF, OR: V, VT (T_Vstatic has no visual flow)
stim_conditions = containers.Map();
stim_conditions('Speed') = {'T_Vstatic', 'VT'};
stim_conditions('TF') = {'V', 'VT'};
stim_conditions('SF') = {'V', 'VT'};
stim_conditions('OR') = {'V', 'VT'};

for row_i = 1:length(model_names_plot)
    mn = model_names_plot{row_i};
    
    for col_i = 1:length(stimulus_types_plot)
        stim_name = stimulus_types_plot{col_i};
        conds = stim_conditions(stim_name);
        
        subplot(3, 4, (row_i-1)*4 + col_i);
        hold on;
        
        x_tick_positions = [];
        x_tick_labels = {};
        
        % Plot each condition (2 boxes per subplot)
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
            x_tick_labels{end+1} = cond_labels(cond_name); %#ok<AGROW>
            
            % Box plot
            q = quantile(valid_r, [0.25, 0.5, 0.75]);
            iqr_val = q(3) - q(1);
            whisker_lo = max(min(valid_r), q(1) - 1.5*iqr_val);
            whisker_hi = min(max(valid_r), q(3) + 1.5*iqr_val);
            
            box_width = 0.6;
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
            
            % --- Jittered points with model-category-specific markers ---
            % Find valid indices (non-NaN) to get model categories
            valid_idx = find(~isnan(r_vals));
            valid_r = r_vals(valid_idx);
            jitter = 0.15 * randn(size(valid_r));
            
            % Plot each model category with its own marker
            for mci = 1:length(all_model_cats)
                mc = all_model_cats{mci};
                mc_mask = strcmp(tc_corr_results.model_category(valid_idx), mc);
                if ~any(mc_mask), continue; end
                
                if model_cat_markers.isKey(mc)
                    mkr = model_cat_markers(mc);
                else
                    mkr = 'o';  % fallback
                end
                
                % Plot with marker style
                scatter(x_pos + jitter(mc_mask), valid_r(mc_mask), 30, clr, mkr, ...
                    'MarkerFaceColor', clr, 'MarkerEdgeColor', clr, ...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7, 'LineWidth', 0.8);
            end
            
            % Mean marker (diamond)
            mean_val = mean(valid_r);
            plot(x_pos, mean_val, 'd', 'MarkerSize', 8, ...
                'MarkerFaceColor', 'w', 'MarkerEdgeColor', clr, 'LineWidth', 1.5);
            
            % Print stats
            fprintf('  %s %s %s: n=%d, mean=%.3f, median=%.3f\n', ...
                mn, stim_name, cond_name, length(valid_r), mean_val, q(2));
        end
        
        % Reference line at r=0
        plot([0.5, 2.5], [0, 0], 'k--', 'LineWidth', 0.5);
        
        hold off;
        xlim([0.5, 2.5]);
        ylim([-1, 1]);
        set(gca, 'XTick', x_tick_positions, 'XTickLabel', x_tick_labels);
        
        % Y-axis label only on first column
        if col_i == 1
            ylabel('Pearson r');
        end
        
        % Title: stimulus name on first row
        if row_i == 1
            title(sprintf('%s', stim_name), 'FontSize', 11, 'FontWeight', 'bold');
        end
        
        % Row label on first column
        if col_i == 1
            yl = ylabel({model_labels{row_i}, 'Pearson r'}, 'FontSize', 10);
            yl.Position(1) = yl.Position(1) - 0.3;
        end
        
        set(gca, 'box', 'off', 'FontSize', 9);
    end
end

% Add legend at bottom: condition colors (left) and model type markers (right)
% --- Condition colors legend (left) ---
ax_leg_cond = axes(fig_corr, 'Position', [0.08 0.01 0.25 0.03], 'Visible', 'off');
hold(ax_leg_cond, 'on');
h_leg_cond = gobjects(3, 1);
h_leg_cond(1) = plot(ax_leg_cond, NaN, NaN, 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', cond_colors('T_Vstatic'), 'MarkerEdgeColor', cond_colors('T_Vstatic'));
h_leg_cond(2) = plot(ax_leg_cond, NaN, NaN, 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', cond_colors('V'), 'MarkerEdgeColor', cond_colors('V'));
h_leg_cond(3) = plot(ax_leg_cond, NaN, NaN, 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', cond_colors('VT'), 'MarkerEdgeColor', cond_colors('VT'));
legend(h_leg_cond, {'T', 'V', 'VT'}, ...
    'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 8);
title(ax_leg_cond, 'Condition:', 'FontSize', 8);
hold(ax_leg_cond, 'off');

% --- Model type markers legend (right) ---
ax_leg_model = axes(fig_corr, 'Position', [0.35 0.01 0.60 0.03], 'Visible', 'off');
hold(ax_leg_model, 'on');

% Create legend entries for model categories present in data
n_cats = length(all_model_cats);
h_leg_model = gobjects(n_cats, 1);
model_cat_labels = cell(n_cats, 1);
neutral_clr = [0.3 0.3 0.3];  % Gray for legend markers

for mci = 1:n_cats
    mc = all_model_cats{mci};
    if model_cat_markers.isKey(mc)
        mkr = model_cat_markers(mc);
    else
        mkr = 'o';
    end
    h_leg_model(mci) = plot(ax_leg_model, NaN, NaN, mkr, 'MarkerSize', 8, ...
        'MarkerFaceColor', neutral_clr, 'MarkerEdgeColor', neutral_clr, 'LineWidth', 1);
    model_cat_labels{mci} = mc;
end

legend(h_leg_model, model_cat_labels, ...
    'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 7, 'NumColumns', min(n_cats, 5));
title(ax_leg_model, 'Winner Model:', 'FontSize', 8);
hold(ax_leg_model, 'off');

% Save figure
if save_figs
    fname_corr = 'tuning_curve_correlations_summary.png';
    drawnow;
    exportgraphics(fig_corr, fullfile(ctl.figs.curr_dir, fname_corr), 'Resolution', 200);
    fprintf('  Saved: %s\n', fname_corr);
end
close(fig_corr);

% Store correlations in a CSV for later analysis
if csv_output
    fprintf('\n--- Exporting tuning curve correlations CSV ---\n');
    
    % Build table
    T_corr = table();
    T_corr.probe_id = tc_corr_results.probe_id(:);
    T_corr.cluster_id = tc_corr_results.cluster_id(:);
    T_corr.category = tc_corr_results.category(:);
    T_corr.selected_vars_str = tc_corr_results.selected_vars_str(:);
    T_corr.model_category = tc_corr_results.model_category(:);
    
    % Add all correlation fields
    corr_fields = fieldnames(tc_corr_results);
    for fi = 1:length(corr_fields)
        fn = corr_fields{fi};
        if startsWith(fn, 'Selected_') || startsWith(fn, 'Additive_') || startsWith(fn, 'FullInteraction_')
            T_corr.(fn) = tc_corr_results.(fn)(:);
        end
    end
    
    csv_dir = fullfile(ctl.figs.curr_dir, 'csv');
    if ~exist(csv_dir, 'dir'), mkdir(csv_dir); end
    csv_corr_path = fullfile(csv_dir, 'tuning_curve_correlations.csv');
    writetable(T_corr, csv_corr_path);
    fprintf('  Saved: tuning_curve_correlations.csv (%d rows, %d columns)\n', ...
        height(T_corr), width(T_corr));
end

%% ====================================================================
%  Section 8g: Selected vs Additive Model Comparison (Paired)
%  ====================================================================
%  Connected scatter plot comparing Selected (winner) vs Additive models.
%  Lines connect the same cluster's performance across both models.
%  
%  Layout: 2x2
%    Columns: 1) Pearson correlation, 2) R² (variance explained)
%    Rows: 1) Speed tuning, 2) TF tuning
%  
%  X-axis: Two categories (Selected, Additive)
%  Y-axis: Correlation or R²
%  
%  Markers: Winner model category on Selected side, plain circles on Additive
%  Colors: Condition (T=red, V=green, VT=blue)
%  ====================================================================
fprintf('\n--- Generating Selected vs Additive paired comparison ---\n');

fig_paired = figure('Position', [50 50 1000 800], 'Visible', 'off');
sgtitle('Selected vs Additive Model: Paired Comparison', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Colors for each condition (matching other plots)
clr_T = cond_colors('T_Vstatic');  % Red
clr_V = cond_colors('V');          % Green
clr_VT = cond_colors('VT');        % Blue

% X positions for the two model types
x_selected = 1;
x_additive = 2;

% Define subplot layout: [stim_type, metric_type, field_cond1, color1, label1, field_cond2, color2, label2]
% Row 1 (Speed): T_Vstatic and VT (both have locomotion)
% Row 2 (TF): V and VT (both have visual flow)
subplot_config = {
    {'Speed', 'Pearson r', 'Speed_T_Vstatic', clr_T, 'T', 'Speed_VT', clr_VT, 'VT'};      % (1,1)
    {'Speed', 'R²', 'Speed_T_Vstatic', clr_T, 'T', 'Speed_VT', clr_VT, 'VT'};             % (1,2)
    {'TF', 'Pearson r', 'TF_V', clr_V, 'V', 'TF_VT', clr_VT, 'VT'};                        % (2,1)
    {'TF', 'R²', 'TF_V', clr_V, 'V', 'TF_VT', clr_VT, 'VT'}                                 % (2,2)
};

for sp_i = 1:4
    subplot(2, 2, sp_i);
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
    
    % Get data for Selected and Additive models
    selected_cond1 = tc_corr_results.(sprintf('Selected_%s', field_cond1));
    additive_cond1 = tc_corr_results.(sprintf('Additive_%s', field_cond1));
    selected_cond2 = tc_corr_results.(sprintf('Selected_%s', field_cond2));
    additive_cond2 = tc_corr_results.(sprintf('Additive_%s', field_cond2));
    
    % Convert to R² if needed
    if is_r2
        selected_cond1 = selected_cond1.^2;
        additive_cond1 = additive_cond1.^2;
        selected_cond2 = selected_cond2.^2;
        additive_cond2 = additive_cond2.^2;
    end
    
    % ====== Condition 1 (T or V depending on row) ======
    valid_c1 = ~isnan(selected_cond1) & ~isnan(additive_cond1);
    valid_idx_c1 = find(valid_c1);
    
    % Box width and scatter offset
    box_w = 0.18;
    scatter_offset = 0.22;  % offset for scatter points from box center
    
    if sum(valid_c1) > 0
        % Draw connecting lines (thin, between scatter points)
        jitter_c1 = 0.03 * randn(sum(valid_c1), 1);
        for ii = 1:length(valid_idx_c1)
            ci = valid_idx_c1(ii);
            plot([x_selected - scatter_offset + jitter_c1(ii), x_additive + scatter_offset + jitter_c1(ii)], ...
                [selected_cond1(ci), additive_cond1(ci)], ...
                '-', 'Color', [clr_cond1, 0.15], 'LineWidth', 0.3);
        end
        
        % --- Boxplot for Condition 1, Selected side ---
        q_sel_c1 = quantile(selected_cond1(valid_c1), [0.25, 0.5, 0.75]);
        iqr_sel_c1 = q_sel_c1(3) - q_sel_c1(1);
        whi_lo_sel_c1 = max(min(selected_cond1(valid_c1)), q_sel_c1(1) - 1.5*iqr_sel_c1);
        whi_hi_sel_c1 = min(max(selected_cond1(valid_c1)), q_sel_c1(3) + 1.5*iqr_sel_c1);
        % Box
        rectangle('Position', [x_selected - box_w/2, q_sel_c1(1), box_w, max(q_sel_c1(3)-q_sel_c1(1), 0.01)], ...
            'EdgeColor', clr_cond1, 'FaceColor', [clr_cond1, 0.2], 'LineWidth', 1.2);
        % Median
        plot([x_selected - box_w/2, x_selected + box_w/2], [q_sel_c1(2), q_sel_c1(2)], '-', 'Color', clr_cond1, 'LineWidth', 2);
        % Whiskers
        plot([x_selected, x_selected], [whi_lo_sel_c1, q_sel_c1(1)], '-', 'Color', clr_cond1, 'LineWidth', 1);
        plot([x_selected, x_selected], [q_sel_c1(3), whi_hi_sel_c1], '-', 'Color', clr_cond1, 'LineWidth', 1);
        plot([x_selected - box_w/4, x_selected + box_w/4], [whi_lo_sel_c1, whi_lo_sel_c1], '-', 'Color', clr_cond1, 'LineWidth', 1);
        plot([x_selected - box_w/4, x_selected + box_w/4], [whi_hi_sel_c1, whi_hi_sel_c1], '-', 'Color', clr_cond1, 'LineWidth', 1);
        
        % --- Boxplot for Condition 1, Additive side ---
        q_add_c1 = quantile(additive_cond1(valid_c1), [0.25, 0.5, 0.75]);
        iqr_add_c1 = q_add_c1(3) - q_add_c1(1);
        whi_lo_add_c1 = max(min(additive_cond1(valid_c1)), q_add_c1(1) - 1.5*iqr_add_c1);
        whi_hi_add_c1 = min(max(additive_cond1(valid_c1)), q_add_c1(3) + 1.5*iqr_add_c1);
        % Box
        rectangle('Position', [x_additive - box_w/2, q_add_c1(1), box_w, max(q_add_c1(3)-q_add_c1(1), 0.01)], ...
            'EdgeColor', clr_cond1, 'FaceColor', [clr_cond1, 0.2], 'LineWidth', 1.2);
        % Median
        plot([x_additive - box_w/2, x_additive + box_w/2], [q_add_c1(2), q_add_c1(2)], '-', 'Color', clr_cond1, 'LineWidth', 2);
        % Whiskers
        plot([x_additive, x_additive], [whi_lo_add_c1, q_add_c1(1)], '-', 'Color', clr_cond1, 'LineWidth', 1);
        plot([x_additive, x_additive], [q_add_c1(3), whi_hi_add_c1], '-', 'Color', clr_cond1, 'LineWidth', 1);
        plot([x_additive - box_w/4, x_additive + box_w/4], [whi_lo_add_c1, whi_lo_add_c1], '-', 'Color', clr_cond1, 'LineWidth', 1);
        plot([x_additive - box_w/4, x_additive + box_w/4], [whi_hi_add_c1, whi_hi_add_c1], '-', 'Color', clr_cond1, 'LineWidth', 1);
        
        % --- Scatter points on the side of boxes ---
        for mci = 1:length(all_model_cats)
            mc = all_model_cats{mci};
            mc_mask = strcmp(tc_corr_results.model_category(valid_idx_c1), mc);
            if ~any(mc_mask), continue; end
            
            if model_cat_markers.isKey(mc)
                mkr = model_cat_markers(mc);
            else
                mkr = 'o';
            end
            
            mc_idx = valid_idx_c1(mc_mask);
            % Selected side (left of box)
            scatter(x_selected - scatter_offset + jitter_c1(mc_mask), selected_cond1(mc_idx), ...
                15, clr_cond1, mkr, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', clr_cond1, 'LineWidth', 0.6);
            % Additive side (right of box)
            scatter(x_additive + scatter_offset + jitter_c1(mc_mask), additive_cond1(mc_idx), ...
                15, clr_cond1, mkr, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', clr_cond1, 'LineWidth', 0.6);
        end
    end
    
    % ====== Condition 2 (VT - blue) ======
    valid_c2 = ~isnan(selected_cond2) & ~isnan(additive_cond2);
    valid_idx_c2 = find(valid_c2);
    
    if sum(valid_c2) > 0
        % Draw connecting lines
        jitter_c2 = 0.03 * randn(sum(valid_c2), 1);
        for ii = 1:length(valid_idx_c2)
            ci = valid_idx_c2(ii);
            plot([x_selected - scatter_offset + jitter_c2(ii), x_additive + scatter_offset + jitter_c2(ii)], ...
                [selected_cond2(ci), additive_cond2(ci)], ...
                '-', 'Color', [clr_cond2, 0.15], 'LineWidth', 0.3);
        end
        
        % --- Boxplot for Condition 2, Selected side ---
        q_sel_c2 = quantile(selected_cond2(valid_c2), [0.25, 0.5, 0.75]);
        iqr_sel_c2 = q_sel_c2(3) - q_sel_c2(1);
        whi_lo_sel_c2 = max(min(selected_cond2(valid_c2)), q_sel_c2(1) - 1.5*iqr_sel_c2);
        whi_hi_sel_c2 = min(max(selected_cond2(valid_c2)), q_sel_c2(3) + 1.5*iqr_sel_c2);
        % Box
        rectangle('Position', [x_selected - box_w/2, q_sel_c2(1), box_w, max(q_sel_c2(3)-q_sel_c2(1), 0.01)], ...
            'EdgeColor', clr_cond2, 'FaceColor', [clr_cond2, 0.2], 'LineWidth', 1.2);
        % Median
        plot([x_selected - box_w/2, x_selected + box_w/2], [q_sel_c2(2), q_sel_c2(2)], '-', 'Color', clr_cond2, 'LineWidth', 2);
        % Whiskers
        plot([x_selected, x_selected], [whi_lo_sel_c2, q_sel_c2(1)], '-', 'Color', clr_cond2, 'LineWidth', 1);
        plot([x_selected, x_selected], [q_sel_c2(3), whi_hi_sel_c2], '-', 'Color', clr_cond2, 'LineWidth', 1);
        plot([x_selected - box_w/4, x_selected + box_w/4], [whi_lo_sel_c2, whi_lo_sel_c2], '-', 'Color', clr_cond2, 'LineWidth', 1);
        plot([x_selected - box_w/4, x_selected + box_w/4], [whi_hi_sel_c2, whi_hi_sel_c2], '-', 'Color', clr_cond2, 'LineWidth', 1);
        
        % --- Boxplot for Condition 2, Additive side ---
        q_add_c2 = quantile(additive_cond2(valid_c2), [0.25, 0.5, 0.75]);
        iqr_add_c2 = q_add_c2(3) - q_add_c2(1);
        whi_lo_add_c2 = max(min(additive_cond2(valid_c2)), q_add_c2(1) - 1.5*iqr_add_c2);
        whi_hi_add_c2 = min(max(additive_cond2(valid_c2)), q_add_c2(3) + 1.5*iqr_add_c2);
        % Box
        rectangle('Position', [x_additive - box_w/2, q_add_c2(1), box_w, max(q_add_c2(3)-q_add_c2(1), 0.01)], ...
            'EdgeColor', clr_cond2, 'FaceColor', [clr_cond2, 0.2], 'LineWidth', 1.2);
        % Median
        plot([x_additive - box_w/2, x_additive + box_w/2], [q_add_c2(2), q_add_c2(2)], '-', 'Color', clr_cond2, 'LineWidth', 2);
        % Whiskers
        plot([x_additive, x_additive], [whi_lo_add_c2, q_add_c2(1)], '-', 'Color', clr_cond2, 'LineWidth', 1);
        plot([x_additive, x_additive], [q_add_c2(3), whi_hi_add_c2], '-', 'Color', clr_cond2, 'LineWidth', 1);
        plot([x_additive - box_w/4, x_additive + box_w/4], [whi_lo_add_c2, whi_lo_add_c2], '-', 'Color', clr_cond2, 'LineWidth', 1);
        plot([x_additive - box_w/4, x_additive + box_w/4], [whi_hi_add_c2, whi_hi_add_c2], '-', 'Color', clr_cond2, 'LineWidth', 1);
        
        % --- Scatter points on the side of boxes ---
        for mci = 1:length(all_model_cats)
            mc = all_model_cats{mci};
            mc_mask = strcmp(tc_corr_results.model_category(valid_idx_c2), mc);
            if ~any(mc_mask), continue; end
            
            if model_cat_markers.isKey(mc)
                mkr = model_cat_markers(mc);
            else
                mkr = 'o';
            end
            
            mc_idx = valid_idx_c2(mc_mask);
            % Selected side (left of box)
            scatter(x_selected - scatter_offset + jitter_c2(mc_mask), selected_cond2(mc_idx), ...
                15, clr_cond2, mkr, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', clr_cond2, 'LineWidth', 0.6);
            % Additive side (right of box)
            scatter(x_additive + scatter_offset + jitter_c2(mc_mask), additive_cond2(mc_idx), ...
                15, clr_cond2, mkr, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', clr_cond2, 'LineWidth', 0.6);
        end
    end
    
    % Reference line at 0 (for correlation) or formatting
    if is_r2
        plot([0.4, 2.6], [0, 0], 'k--', 'LineWidth', 0.5);
        ylim([0, 1]);
        ylabel('R² (variance explained)');
    else
        plot([0.4, 2.6], [0, 0], 'k--', 'LineWidth', 0.5);
        ylim([-1, 1]);
        ylabel('Pearson r');
    end
    
    hold off;
    xlim([0.4, 2.6]);
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'Selected', 'Additive'});
    title(sprintf('%s - %s', stim_type, metric_type), 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'box', 'off', 'FontSize', 9);
    
    % Print stats
    if sum(valid_c1) > 0
        delta_c1 = selected_cond1(valid_c1) - additive_cond1(valid_c1);
        fprintf('  %s %s %s: n=%d, mean delta=%.3f (Selected - Additive)\n', ...
            stim_type, metric_type, label_cond1, sum(valid_c1), mean(delta_c1));
    end
    if sum(valid_c2) > 0
        delta_c2 = selected_cond2(valid_c2) - additive_cond2(valid_c2);
        fprintf('  %s %s %s: n=%d, mean delta=%.3f (Selected - Additive)\n', ...
            stim_type, metric_type, label_cond2, sum(valid_c2), mean(delta_c2));
    end
end

% Add combined legend at bottom
% --- Condition legend (left) ---
ax_leg_paired_cond = axes(fig_paired, 'Position', [0.05 0.01 0.25 0.04], 'Visible', 'off');
hold(ax_leg_paired_cond, 'on');
h_cond_p = gobjects(3, 1);
h_cond_p(1) = plot(ax_leg_paired_cond, NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr_T, 'LineWidth', 1);
h_cond_p(2) = plot(ax_leg_paired_cond, NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr_V, 'LineWidth', 1);
h_cond_p(3) = plot(ax_leg_paired_cond, NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr_VT, 'LineWidth', 1);
legend(h_cond_p, {'T', 'V', 'VT'}, ...
    'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 8);
title(ax_leg_paired_cond, 'Condition:', 'FontSize', 8);
hold(ax_leg_paired_cond, 'off');

% --- Model type markers legend (right) ---
ax_leg_paired_model = axes(fig_paired, 'Position', [0.32 0.01 0.63 0.04], 'Visible', 'off');
hold(ax_leg_paired_model, 'on');

n_cats_p = length(all_model_cats);
h_leg_p = gobjects(n_cats_p, 1);
model_cat_labels_p = cell(n_cats_p, 1);
neutral_clr_p = [0.3 0.3 0.3];

for mci = 1:n_cats_p
    mc = all_model_cats{mci};
    if model_cat_markers.isKey(mc)
        mkr = model_cat_markers(mc);
    else
        mkr = 'o';
    end
    h_leg_p(mci) = plot(ax_leg_paired_model, NaN, NaN, mkr, 'MarkerSize', 6, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', neutral_clr_p, 'LineWidth', 1);
    model_cat_labels_p{mci} = mc;
end

legend(h_leg_p, model_cat_labels_p, ...
    'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 7, 'NumColumns', min(n_cats_p, 5));
title(ax_leg_paired_model, 'Winner Model:', 'FontSize', 8);
hold(ax_leg_paired_model, 'off');

% Save figure
if save_figs
    fname_paired = 'selected_vs_additive_paired.png';
    drawnow;
    exportgraphics(fig_paired, fullfile(ctl.figs.curr_dir, fname_paired), 'Resolution', 200);
    fprintf('  Saved: %s\n', fname_paired);
end
close(fig_paired);

%% ====================================================================
%  Section 8b: Trial-level GLM predictions (inspect variability)
%  ====================================================================
% This section visualizes GLM predictions at the single-trial level to
% inspect how well the model captures firing rate dynamics in time.
%
% Cluster selection: Top 4 clusters by Selected_cv_bps with >= 2 variables
% Trial selection per cluster: 2 T_Vstatic, 2 V, 4 VT (varying SF/OR combos)
%
% Figure layout (8 rows x 3 cols):
%   Row = trial, Col1 = Speed(t), Col2 = TF(t), Col3 = Observed & Predicted FR(t)

fprintf('\n--- Trial-level GLM predictions ---\n');

% Use the time-bin GLM results
gt_tag_trial = 'time';
ml_trial = 'Selected';  % Use the forward-selected model

% --- Select top 4 clusters with >= 2 selected variables ---
n_vars_all = results.(gt_tag_trial).n_selected_vars;
cv_bps_all = results.(gt_tag_trial).Selected_cv_bps;

% Filter: need at least 2 variables selected (e.g., Speed + TF)
eligible_mask = n_vars_all >= 2 & ~isnan(cv_bps_all);
eligible_idx = find(eligible_mask);

if isempty(eligible_idx)
    fprintf('  No clusters with >= 2 selected variables found. Skipping trial-level plots.\n');
else
    % Sort eligible clusters by cv_bps (descending)
    [~, sort_order] = sort(cv_bps_all(eligible_idx), 'descend');
    top_cluster_idx = eligible_idx(sort_order(1:min(4, length(sort_order))));
    n_top_clusters = length(top_cluster_idx);
    
    fprintf('  Selected %d top clusters for trial-level visualization:\n', n_top_clusters);
    for tci = 1:n_top_clusters
        ci = top_cluster_idx(tci);
        pid_tc = results.(gt_tag_trial).probe_id{ci};
        cid_tc = results.(gt_tag_trial).cluster_id(ci);
        vars_str = results.(gt_tag_trial).selected_vars_str{ci};
        fprintf('    %d. %s cluster %d: cv_bps=%.3f, vars=[%s]\n', ...
            tci, pid_tc, cid_tc, cv_bps_all(ci), vars_str);
    end
    
    % --- Create one figure per cluster ---
    for tci = 1:n_top_clusters
        ci = top_cluster_idx(tci);
        pid_tc = results.(gt_tag_trial).probe_id{ci};
        cid_tc = results.(gt_tag_trial).cluster_id(ci);
        vars_str = results.(gt_tag_trial).selected_vars_str{ci};
        
        % Get cluster data
        idx_cl = T_master_time.probe_id == pid_tc & T_master_time.cluster_id == cid_tc;
        T_cluster_trial = T_master_time(idx_cl, :);
        
        % Get predictions (aligned with T_cluster_trial rows)
        if isempty(cluster_predictions.(gt_tag_trial){ci}) || ...
                ~isfield(cluster_predictions.(gt_tag_trial){ci}, ml_trial)
            fprintf('    Skipping cluster %d: no predictions available\n', ci);
            continue;
        end
        cv_pred_fr_cl = cluster_predictions.(gt_tag_trial){ci}.(ml_trial).cv_predicted_fr;
        
        % --- Select 8 trials: 2 T_Vstatic, 2 V, 4 VT (varying SF/OR) ---
        trial_info_all = unique(T_cluster_trial(:, {'trial_id', 'condition', 'sf', 'orientation'}), 'rows');
        
        % Get trials by condition
        t_trials = trial_info_all(trial_info_all.condition == "T_Vstatic", :);
        v_trials = trial_info_all(trial_info_all.condition == "V", :);
        vt_trials = trial_info_all(trial_info_all.condition == "VT", :);
        
        selected_trials = [];
        selected_conds = {};
        selected_labels = {};
        
        % 1. Select 2 T_Vstatic trials (random)
        n_t_select = min(2, height(t_trials));
        if n_t_select >= 1
            t_perm = randperm(height(t_trials), n_t_select);
            for ti = 1:n_t_select
                selected_trials(end+1) = t_trials.trial_id(t_perm(ti));
                selected_conds{end+1} = 'T_Vstatic';
                selected_labels{end+1} = 'T (speed only)';
            end
        end
        
        % 2. Select 2 V trials (random)
        n_v_select = min(2, height(v_trials));
        if n_v_select >= 1
            v_perm = randperm(height(v_trials), n_v_select);
            for vi = 1:n_v_select
                sf_val_tr = v_trials.sf(v_perm(vi));
                or_val_tr = v_trials.orientation(v_perm(vi));
                selected_trials(end+1) = v_trials.trial_id(v_perm(vi));
                selected_conds{end+1} = 'V';
                selected_labels{end+1} = sprintf('V: SF=%.3f, OR=%.2f', sf_val_tr, or_val_tr);
            end
        end
        
        % 3. Select 4 VT trials (try to vary SF/OR combinations)
        n_vt_select = min(4, height(vt_trials));
        if n_vt_select >= 1
            % Get unique SF/OR combos
            sf_or_combos = unique(vt_trials(:, {'sf', 'orientation'}), 'rows');
            n_combos = height(sf_or_combos);
            
            if n_combos >= n_vt_select
                % Select from different combos
                combo_perm = randperm(n_combos, n_vt_select);
                for vti = 1:n_vt_select
                    combo = sf_or_combos(combo_perm(vti), :);
                    matching_trials = vt_trials(vt_trials.sf == combo.sf & ...
                        vt_trials.orientation == combo.orientation, :);
                    if height(matching_trials) >= 1
                        rng_idx = randi(height(matching_trials));
                        selected_trials(end+1) = matching_trials.trial_id(rng_idx);
                        selected_conds{end+1} = 'VT';
                        selected_labels{end+1} = sprintf('VT: SF=%.3f, OR=%.2f', combo.sf, combo.orientation);
                    end
                end
            else
                % Fewer combos than desired, pick random VT trials
                vt_perm = randperm(height(vt_trials), n_vt_select);
                for vti = 1:n_vt_select
                    sf_val_tr = vt_trials.sf(vt_perm(vti));
                    or_val_tr = vt_trials.orientation(vt_perm(vti));
                    selected_trials(end+1) = vt_trials.trial_id(vt_perm(vti));
                    selected_conds{end+1} = 'VT';
                    selected_labels{end+1} = sprintf('VT: SF=%.3f, OR=%.2f', sf_val_tr, or_val_tr);
                end
            end
        end
        
        n_selected_trials = length(selected_trials);
        if n_selected_trials < 1
            fprintf('    Skipping cluster %d: no valid trials found\n', ci);
            continue;
        end
        
        % --- Create figure: n_trials rows x 3 columns ---
        n_rows_trial = n_selected_trials;
        n_cols_trial = 3;
        fig_trial = figure('Position', [50 50 1000 150*n_rows_trial], ...
            'Visible', 'off', 'PaperOrientation', 'landscape');
        
        sgtitle(sprintf('Trial-Level Predictions | %s Cluster %d | cv\\_bps=%.3f | Vars: %s', ...
            pid_tc, cid_tc, cv_bps_all(ci), vars_str), 'FontSize', 11, 'FontWeight', 'bold');
        
        % Color map for conditions
        cond_colors = containers.Map();
        cond_colors('T_Vstatic') = [0.2 0.6 0.2];  % green
        cond_colors('V') = [0.8 0.2 0.2];           % red
        cond_colors('VT') = [0.2 0.2 0.8];          % blue
        
        % First pass: collect data and find global TF max for consistent y-axis
        trial_data = cell(n_selected_trials, 1);
        tf_max_global = 0;
        
        % Smoothing kernel for observed FR (boxcar, Park 2014 style)
        % Number of bins in the smoothing window
        n_smooth_bins = round(obs_fr_smooth_width / time_bin_width);
        if n_smooth_bins < 1, n_smooth_bins = 1; end
        boxcar_kernel = ones(n_smooth_bins, 1) / n_smooth_bins;
        
        for tri = 1:n_selected_trials
            trial_id_sel = selected_trials(tri);
            trial_mask = T_cluster_trial.trial_id == trial_id_sel;
            T_trial = T_cluster_trial(trial_mask, :);
            pred_fr_trial = cv_pred_fr_cl(trial_mask);
            
            % Sort by time to avoid line artifacts
            [time_sorted, sort_idx] = sort(T_trial.time_in_trial);
            trial_data{tri}.time = time_sorted;
            trial_data{tri}.speed = T_trial.speed(sort_idx);
            trial_data{tri}.tf = T_trial.tf(sort_idx);
            
            % Smooth observed spike counts with boxcar, then convert to rate
            % This matches Park 2014: "smoothing the spike train with a 100-ms boxcar"
            spike_counts_sorted = T_trial.spike_count(sort_idx);
            obs_fr_raw = spike_counts_sorted / time_bin_width;  % instantaneous rate
            obs_fr_smooth = conv(obs_fr_raw, boxcar_kernel, 'same');  % smoothed rate
            trial_data{tri}.obs_fr = obs_fr_smooth;
            trial_data{tri}.pred_fr = pred_fr_trial(sort_idx);
            
            tf_max_global = max(tf_max_global, max(trial_data{tri}.tf));
        end
        tf_ylim = [0, tf_max_global * 1.1];
        
        % Second pass: create plots
        for tri = 1:n_selected_trials
            trial_id_sel = selected_trials(tri);
            cond_sel = selected_conds{tri};
            label_sel = selected_labels{tri};
            
            % Get pre-sorted data
            time_vec = trial_data{tri}.time;
            speed_vec_tr = trial_data{tri}.speed;
            tf_vec_tr = trial_data{tri}.tf;
            obs_fr_tr = trial_data{tri}.obs_fr;
            pred_fr_tr = trial_data{tri}.pred_fr;
            
            % Get condition color
            if cond_colors.isKey(cond_sel)
                cond_clr = cond_colors(cond_sel);
            else
                cond_clr = [0.3 0.3 0.3];
            end
            
            % --- Col 1: Speed vs time ---
            ax1 = subplot(n_rows_trial, n_cols_trial, (tri-1)*n_cols_trial + 1);
            plot(time_vec, speed_vec_tr, '-', 'Color', cond_clr, 'LineWidth', 1.2);
            ylabel('Speed (cm/s)');
            if tri == 1, title('Speed', 'FontWeight', 'bold'); end
            if tri == n_rows_trial, xlabel('Time (s)'); end
            set(ax1, 'FontSize', 8, 'Box', 'off');
            xlim([min(time_vec)-0.05 max(time_vec)+0.05]);
            
            % Add trial label on left
            text(-0.15, 0.5, sprintf('Trial %d\n%s', trial_id_sel, label_sel), ...
                'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle', 'FontSize', 7, 'Rotation', 0);
            
            % --- Col 2: TF vs time (consistent y-axis) ---
            ax2 = subplot(n_rows_trial, n_cols_trial, (tri-1)*n_cols_trial + 2);
            plot(time_vec, tf_vec_tr, '-', 'Color', cond_clr, 'LineWidth', 1.2);
            ylabel('TF (Hz)');
            if tri == 1, title('Temporal Freq', 'FontWeight', 'bold'); end
            if tri == n_rows_trial, xlabel('Time (s)'); end
            set(ax2, 'FontSize', 8, 'Box', 'off');
            xlim([min(time_vec)-0.05 max(time_vec)+0.05]);
            ylim(tf_ylim);  % Use consistent TF y-axis across all trials
            
            % --- Col 3: Observed & Predicted firing rate vs time (lines) ---
            ax3 = subplot(n_rows_trial, n_cols_trial, (tri-1)*n_cols_trial + 3);
            hold on;
            plot(time_vec, obs_fr_tr, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
            plot(time_vec, pred_fr_tr, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
            ylabel('FR (Hz)');
            if tri == 1
                title('Firing Rate', 'FontWeight', 'bold');
                legend({'Observed', 'Predicted'}, 'Location', 'northeast', 'FontSize', 6);
            end
            if tri == n_rows_trial, xlabel('Time (s)'); end
            hold off;
            set(ax3, 'FontSize', 8, 'Box', 'off');
            xlim([min(time_vec)-0.05 max(time_vec)+0.05]);
            max_fr = max([max(obs_fr_tr), max(pred_fr_tr), 1]);
            ylim([0 max_fr * 1.1]);
        end
        
        % Save figure
        if save_figs
            fname_trial = sprintf('trial_level_%s_%d.png', pid_tc, cid_tc);
            drawnow;
            exportgraphics(fig_trial, fullfile(ctl.figs.curr_dir, fname_trial), 'Resolution', 200);
            fprintf('    Saved: %s\n', fname_trial);
        end
        close(fig_trial);
    end
end

fprintf('--- Trial-level visualization complete ---\n');

%% ====================================================================
%  Section 9: CSV exports
%  ====================================================================
if csv_output
    fprintf('\n--- Exporting CSV files ---\n');
    
    csv_dir = fullfile(ctl.figs.curr_dir, 'csv');
    if ~exist(csv_dir, 'dir'), mkdir(csv_dir); end
    
    % --- CSV 1: Model comparison (Forward Selection GLM) ---
    T_comp = table();
    T_comp.probe_id = results.time.probe_id;
    T_comp.cluster_id = results.time.cluster_id;
    T_comp.n_trials = results.time.n_trials;
    
    % Add pre-filter category from decision tree
    T_comp.prefilter_category = cell(n_unique_clusters, 1);
    for ci_pf = 1:n_unique_clusters
        key = sprintf('%s_%d', results.time.probe_id{ci_pf}, results.time.cluster_id(ci_pf));
        if clusters_for_glm.isKey(key)
            pf_info = clusters_for_glm(key);
            T_comp.prefilter_category{ci_pf} = pf_info.category;
        else
            T_comp.prefilter_category{ci_pf} = 'unknown';
        end
    end
    
    for gt = 1:length(glm_types)
        gt_tag = glm_types{gt};
        prefix = [gt_tag '_'];
        
        T_comp.([prefix 'n_bins']) = results.(gt_tag).n_bins_total;
        T_comp.([prefix 'n_spikes']) = results.(gt_tag).n_spikes_total;
        
        % Forward selection results
        T_comp.([prefix 'selected_vars']) = results.(gt_tag).selected_vars_str;
        T_comp.([prefix 'n_selected_vars']) = results.(gt_tag).n_selected_vars;
        T_comp.([prefix 'selection_rounds']) = results.(gt_tag).selection_rounds;
        
        % CV performance for key models
        T_comp.([prefix 'Null_cv_bps']) = results.(gt_tag).Null_cv_bps;
        T_comp.([prefix 'Selected_cv_bps']) = results.(gt_tag).Selected_cv_bps;
        T_comp.([prefix 'Additive_cv_bps']) = results.(gt_tag).Additive_cv_bps;
        T_comp.([prefix 'FullInteraction_cv_bps']) = results.(gt_tag).FullInteraction_cv_bps;
        
        % Delta comparisons
        T_comp.([prefix 'delta_selected_vs_null']) = results.(gt_tag).delta_bps_selected_vs_null;
        T_comp.([prefix 'delta_additive_vs_null']) = results.(gt_tag).delta_bps_additive_vs_null;
        T_comp.([prefix 'delta_interaction']) = results.(gt_tag).delta_bps_interaction;
        T_comp.([prefix 'delta_selected_vs_additive']) = results.(gt_tag).delta_bps_selected_vs_additive;
        
        % Classification flags (based on forward selection)
        T_comp.([prefix 'is_spprofileeed_tuned']) = results.(gt_tag).is_spprofileeed_tuned;
        T_comp.([prefix 'is_tf_tuned']) = results.(gt_tag).is_tf_tuned;
        T_comp.([prefix 'is_sf_tuned']) = results.(gt_tag).is_sf_tuned;
        T_comp.([prefix 'is_or_tuned']) = results.(gt_tag).is_or_tuned;
        T_comp.([prefix 'has_interaction']) = results.(gt_tag).has_interaction;
        
        % Specific interactions selected
        T_comp.([prefix 'has_spprofileeed_x_tf']) = results.(gt_tag).has_spprofileeed_x_tf;
        T_comp.([prefix 'has_spprofileeed_x_sf']) = results.(gt_tag).has_spprofileeed_x_sf;
        T_comp.([prefix 'has_spprofileeed_x_or']) = results.(gt_tag).has_spprofileeed_x_or;
        T_comp.([prefix 'has_tf_x_sf']) = results.(gt_tag).has_tf_x_sf;
        T_comp.([prefix 'has_tf_x_or']) = results.(gt_tag).has_tf_x_or;
        T_comp.([prefix 'has_sf_x_or']) = results.(gt_tag).has_sf_x_or;
    end
    
    writetable(T_comp, fullfile(csv_dir, 'glm_model_comparison.csv'));
    fprintf('  Saved glm_model_comparison.csv (%d rows)\n', height(T_comp));
    
    % --- CSV 2: Selection history ---
    if ~isempty(all_selection_history)
        T_sel_hist = cell2table(vertcat(all_selection_history{:}), ...
            'VariableNames', {'probe_id', 'cluster_id', 'round', 'best_candidate', ...
            'delta_bps', 'added', 'cv_bps_after'});
        writetable(T_sel_hist, fullfile(csv_dir, 'glm_selection_history.csv'));
        fprintf('  Saved glm_selection_history.csv (%d rows)\n', height(T_sel_hist));
    end
    
    % --- CSV 3: Coefficients ---
    if n_coef_stored > 0
        all_coefficients = all_coefficients(1:n_coef_stored);
        T_coef = cell2table(vertcat(all_coefficients{:}), ...
            'VariableNames', {'probe_id', 'cluster_id', 'glm_type', 'model', ...
            'coefficient', 'estimate', 'se'});
        writetable(T_coef, fullfile(csv_dir, 'glm_coefficients.csv'));
        fprintf('  Saved glm_coefficients.csv (%d rows)\n', height(T_coef));
    end
    
    % --- CSV 4: Trial-level data (time-bin) ---
    writetable(T_master_time, fullfile(csv_dir, 'glm_time_bin_data.csv'));
    fprintf('  Saved glm_time_bin_data.csv (%d rows)\n', height(T_master_time));
    
    % --- Update prefilter_decision_tree.csv with GLM results ---
    % Add forward selection results for each cluster
    fprintf('\n--- Updating prefilter_decision_tree.csv with GLM results ---\n');
    
    % Initialize new columns with NaN/empty for all prefilter clusters
    glm_selected_vars = cell(n_prefilter_clusters, 1);
    glm_n_selected = nan(n_prefilter_clusters, 1);
    glm_selected_cv_bps = nan(n_prefilter_clusters, 1);
    glm_is_spprofileeed_tuned = false(n_prefilter_clusters, 1);
    glm_is_tf_tuned = false(n_prefilter_clusters, 1);
    glm_is_sf_tuned = false(n_prefilter_clusters, 1);
    glm_is_or_tuned = false(n_prefilter_clusters, 1);
    glm_has_interaction = false(n_prefilter_clusters, 1);
    glm_delta_bps_interaction = nan(n_prefilter_clusters, 1);
    glm_additive_cv_bps = nan(n_prefilter_clusters, 1);
    glm_null_cv_bps = nan(n_prefilter_clusters, 1);
    
    % Fill in default values
    for pfi = 1:n_prefilter_clusters
        glm_selected_vars{pfi} = '';
    end
    
    % Match GLM results back to prefilter clusters (using Time-bin GLM)
    for ci = 1:n_unique_clusters
        glm_probe = results.time.probe_id{ci};
        glm_cid = results.time.cluster_id(ci);
        
        % Find matching prefilter row
        for pfi = 1:n_prefilter_clusters
            if strcmp(prefilter_results.probe_id{pfi}, glm_probe) && ...
               prefilter_results.cluster_id(pfi) == glm_cid
                % Copy forward selection results
                glm_selected_vars{pfi} = results.time.selected_vars_str{ci};
                glm_n_selected(pfi) = results.time.n_selected_vars(ci);
                glm_selected_cv_bps(pfi) = results.time.Selected_cv_bps(ci);
                glm_null_cv_bps(pfi) = results.time.Null_cv_bps(ci);
                
                % Copy classification flags (based on forward selection)
                glm_is_spprofileeed_tuned(pfi) = results.time.is_spprofileeed_tuned(ci);
                glm_is_tf_tuned(pfi) = results.time.is_tf_tuned(ci);
                glm_is_sf_tuned(pfi) = results.time.is_sf_tuned(ci);
                glm_is_or_tuned(pfi) = results.time.is_or_tuned(ci);
                glm_has_interaction(pfi) = results.time.has_interaction(ci);
                
                % Copy additional metrics
                glm_delta_bps_interaction(pfi) = results.time.delta_bps_interaction(ci);
                glm_additive_cv_bps(pfi) = results.time.Additive_cv_bps(ci);
                break;
            end
        end
    end
    
    % Update prefilter_table with new columns
    prefilter_table_updated = table(...
        string(prefilter_results.probe_id'), ...
        prefilter_results.cluster_id', ...
        prefilter_results.is_T_significant', ...
        prefilter_results.is_V_significant', ...
        prefilter_results.is_VT_significant', ...
        prefilter_results.n_significant', ...
        prefilter_results.is_spprofileeed_tuned', ...
        prefilter_results.should_run_glm', ...
        string(prefilter_results.category'), ...
        prefilter_results.p_T_svm', ...
        prefilter_results.p_V_svm', ...
        prefilter_results.p_VT_svm', ...
        prefilter_results.direction_T', ...
        prefilter_results.direction_V', ...
        prefilter_results.direction_VT', ...
        string(glm_selected_vars), ...
        glm_n_selected, ...
        glm_null_cv_bps, ...
        glm_selected_cv_bps, ...
        glm_additive_cv_bps, ...
        glm_is_spprofileeed_tuned, ...
        glm_is_tf_tuned, ...
        glm_is_sf_tuned, ...
        glm_is_or_tuned, ...
        glm_has_interaction, ...
        glm_delta_bps_interaction, ...
        'VariableNames', {'probe_id', 'cluster_id', 'is_T_significant', 'is_V_significant', ...
            'is_VT_significant', 'n_significant', 'is_spprofileeed_tuned_prefilter', 'should_run_glm', 'category', ...
            'p_T_svm', 'p_V_svm', 'p_VT_svm', ...
            'direction_T', 'direction_V', 'direction_VT', ...
            'glm_selected_vars', 'glm_n_selected', 'glm_null_cv_bps', 'glm_selected_cv_bps', 'glm_additive_cv_bps', ...
            'glm_is_spprofileeed_tuned', 'glm_is_tf_tuned', 'glm_is_sf_tuned', 'glm_is_or_tuned', ...
            'glm_has_interaction', 'glm_delta_bps_interaction'});
    
    % Save updated prefilter table
    writetable(prefilter_table_updated, prefilter_csv_path);
    fprintf('  Updated prefilter_decision_tree.csv with GLM results (%d rows)\n', height(prefilter_table_updated));
    fprintf('  Added columns: selected vars, CV bps, tuning classification\n');
end

%% ====================================================================
%  Section 10: Summary printout
%  ====================================================================
fprintf('\n====================================================================\n');
fprintf('  GLM SINGLE CLUSTER ANALYSIS (v10: Forward Selection) — COMPLETE\n');
fprintf('====================================================================\n');
fprintf('  Total VISp clusters scanned: %d\n', n_prefilter_clusters);
fprintf('  Clusters after pre-filtering: %d (%.1f%%)\n', n_should_run_glm, 100*n_should_run_glm/n_prefilter_clusters);
fprintf('  Clusters analysed with GLM: %d\n', n_unique_clusters);
fprintf('  Time-bin data points: %d\n', height(T_master_time));
fprintf('  Forward selection: %d candidates, threshold=%.3f bps\n', length(all_candidates), forward_selection_threshold);

fprintf('\n  --- Pre-filter Category Breakdown (GLM clusters) ---\n');
pf_categories = unique(T_comp.prefilter_category);
for pf_i = 1:length(pf_categories)
    pf_cat = pf_categories{pf_i};
    n_pf_cat = sum(strcmp(T_comp.prefilter_category, pf_cat));
    fprintf('    %s: %d (%.1f%%)\n', pf_cat, n_pf_cat, 100*n_pf_cat/n_unique_clusters);
end

for gt = 1:length(glm_types)
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    fprintf('\n  --- %s ---\n', gt_label);
    
    % Selection statistics
    n_selected = results.(gt_tag).n_selected;
    fprintf('  Variables selected per cluster:\n');
    fprintf('    Mean: %.2f, Median: %d, Range: [%d, %d]\n', ...
        mean(n_selected), median(n_selected), min(n_selected), max(n_selected));
    
    % Most selected variables
    fprintf('  Most commonly selected variables:\n');
    all_selected = results.(gt_tag).selected_vars;
    var_counts = containers.Map;
    for ci = 1:length(all_selected)
        sel = all_selected{ci};
        for vi = 1:length(sel)
            v = sel{vi};
            if isKey(var_counts, v)
                var_counts(v) = var_counts(v) + 1;
            else
                var_counts(v) = 1;
            end
        end
    end
    var_names = keys(var_counts);
    var_cnts = cell2mat(values(var_counts));
    [~, sort_idx] = sort(var_cnts, 'descend');
    for vi = 1:min(5, length(var_names))
        v = var_names{sort_idx(vi)};
        fprintf('    %s: %d (%.0f%%)\n', v, var_counts(v), 100*var_counts(v)/n_unique_clusters);
    end
    
    % CV performance
    fprintf('  CV performance (bits/spike):\n');
    fprintf('    Null: median=%.4f\n', median(results.(gt_tag).null_cv_bps));
    fprintf('    Selected: median=%.4f\n', median(results.(gt_tag).selected_cv_bps));
    fprintf('    Additive: median=%.4f\n', median(results.(gt_tag).additive_cv_bps));
    fprintf('  Interaction delta: median = %.4f bps\n', ...
        median(results.(gt_tag).delta_bps_interaction));
    fprintf('  Regularisation: YES (lambda=1.0)\n');
end
fprintf('====================================================================\n');

% --- Close diary and confirm log file ---
fprintf('\nAnalysis completed: %s\n', datestr(now));
fprintf('Full log saved to: %s\n', log_file);
diary off;

%% ====================================================================
%  Local functions
%  ====================================================================

function out = ternary(cond_val, val_true, val_false)
%TERNARY Inline conditional helper
    if cond_val
        out = val_true;
    else
        out = val_false;
    end
end


function [sf_val, vx_val, or_val] = parse_cloud_name(cname, sf_keys_arg, sf_values_map_arg, or_keys_arg, or_values_map_arg)
%PARSE_CLOUD_NAME Extract SF, VX, orientation from a cloud name string
    sf_val = NaN;
    for ki = 1:length(sf_keys_arg)
        if contains(cname, sf_keys_arg{ki})
            sf_val = sf_values_map_arg(sf_keys_arg{ki});
            break;
        end
    end
    
    vx_val = NaN;
    vx_tok = regexp(cname, 'VX(\d+p\d+)', 'tokens');
    if ~isempty(vx_tok)
        vx_str = strrep(vx_tok{1}{1}, 'p', '.');
        vx_val = str2double(vx_str);
    end
    
    % Use startsWith for orientation to avoid matching 'Btheta' bandwidth parameter
    % The main orientation is always at the START of the cloud name
    or_val = NaN;
    for ki = 1:length(or_keys_arg)
        if startsWith(cname, or_keys_arg{ki})
            or_val = or_values_map_arg(or_keys_arg{ki});
            break;
        end
    end
end


function B = make_raised_cosine_basis(x, n_bases, x_min, x_max)
%MAKE_RAISED_COSINE_BASIS Raised cosine bases on log-shifted axis (Weber-law scaling)
    epsilon = 0.5;
    log_x = log(x + epsilon);
    log_min = log(x_min + epsilon);
    log_max = log(x_max + epsilon);
    
    centers = linspace(log_min, log_max, n_bases);
    
    if n_bases > 1
        delta = (log_max - log_min) / (n_bases - 1);
    else
        delta = (log_max - log_min);
    end
    width = delta * 1.5;
    
    B = zeros(length(x), n_bases);
    for bi = 1:n_bases
        z = (log_x - centers(bi)) / width * pi;
        z = max(min(z, pi), -pi);
        B(:, bi) = 0.5 * (1 + cos(z));
    end
end


function B = make_onset_kernel_basis(t_since_onset, n_bases, t_max)
%MAKE_ONSET_KERNEL_BASIS Causal raised cosine bases for onset dynamics (Park et al. 2014 style)
%   Creates raised cosine bases that are NON-ZERO only for t >= 0 (causal).
%   For t < 0 (stationary periods before motion onset), all bases are zero.
%   This captures the transient response at motion onset and adaptation.
%
%   t_since_onset: time since motion onset (negative for stationary, positive for motion)
%   n_bases: number of raised cosine bases
%   t_max: maximum time covered by the kernel (e.g., 2.0 s)
%
%   Following Park et al. 2014: bases are raised cosine bumps on linear axis,
%   spanning from 0 to t_max, with overlapping coverage.

    n = length(t_since_onset);
    B = zeros(n, n_bases);
    
    % Centers span from 0 to t_max
    centers = linspace(0, t_max, n_bases);
    
    if n_bases > 1
        delta = t_max / (n_bases - 1);
    else
        delta = t_max;
    end
    width = delta * 1.5;  % overlap factor
    
    % Mask: only compute for t >= 0 (motion period)
    motion_mask = t_since_onset >= 0;
    t_motion = t_since_onset(motion_mask);
    
    if ~isempty(t_motion)
        for bi = 1:n_bases
            z = (t_motion - centers(bi)) / width * pi;
            z = max(min(z, pi), -pi);
            B(motion_mask, bi) = 0.5 * (1 + cos(z));
        end
    end
    % B remains zero for t < 0 (stationary periods)
end


function [X, col_names] = assemble_design_matrix(B_speed, B_tf, B_onset, sf_vals, or_vals, model_label, sf_ref_levels, or_ref_levels)
%ASSEMBLE_DESIGN_MATRIX Build design matrix from pre-computed basis matrices
%   This avoids recomputing raised cosine bases for each model of the same cluster.
%   B_onset: onset dynamics kernel bases (causal, zero for stationary periods)
%   sf_ref_levels, or_ref_levels: optional reference levels for dummy coding
%   (used for prediction to match training design matrix structure)
%   When ref levels are provided, zero-variance removal is skipped (for prediction).

    n = size(B_speed, 1);
    n_speed_b = size(B_speed, 2);
    n_tf_b = size(B_tf, 2);
    n_onset_b = size(B_onset, 2);
    
    % Determine if this is for prediction (ref levels provided)
    is_prediction = (nargin >= 7 && ~isempty(sf_ref_levels)) || (nargin >= 8 && ~isempty(or_ref_levels));
    
    % Column name lists
    spd_names = arrayfun(@(si) sprintf('Speed_%d', si), 1:n_speed_b, 'UniformOutput', false);
    tf_names_list = arrayfun(@(ti) sprintf('TF_%d', ti), 1:n_tf_b, 'UniformOutput', false);
    onset_names = arrayfun(@(oi) sprintf('Onset_%d', oi), 1:n_onset_b, 'UniformOutput', false);
    
    % SF dummy coding - use reference levels if provided
    if nargin >= 7 && ~isempty(sf_ref_levels)
        unique_sf = sf_ref_levels(:)';
    else
        unique_sf = sort(unique(sf_vals(~isnan(sf_vals))));
    end
    n_sf_levels = length(unique_sf);
    D_sf = zeros(n, max(n_sf_levels - 1, 0));
    sf_names = {};
    for si = 2:n_sf_levels
        D_sf(:, si-1) = double(sf_vals == unique_sf(si));
        sf_names{end+1} = sprintf('SF_%d', si);  %#ok<AGROW>
    end
    
    % Orientation dummy coding - use reference levels if provided
    if nargin >= 8 && ~isempty(or_ref_levels)
        unique_or = or_ref_levels(:)';
    else
        unique_or = sort(unique(or_vals(~isnan(or_vals))));
    end
    n_or_levels = length(unique_or);
    D_or = zeros(n, max(n_or_levels - 1, 0));
    or_names = {};
    for oi = 2:n_or_levels
        D_or(:, oi-1) = double(or_vals == unique_or(oi));
        or_names{end+1} = sprintf('OR_%d', oi);  %#ok<AGROW>
    end
    
    switch model_label
        case 'Null'
            % Null model: intercept + onset kernel only (baseline for forward selection)
            X = [ones(n, 1), B_onset];
            col_names = [{'Intercept'}, onset_names];
        case 'M0'
            % Null model: intercept only (no onset kernel)
            X = ones(n, 1);
            col_names = {'Intercept'};
        case 'M0_Speed'
            % Null + Speed + Onset kernel
            X = [ones(n,1), B_speed, B_onset];
            col_names = [{'Intercept'}, spd_names, onset_names];
        case 'M0_Speed_TF'
            % Null + Speed + TF + Onset kernel
            X = [ones(n,1), B_speed, B_tf, B_onset];
            col_names = [{'Intercept'}, spd_names, tf_names_list, onset_names];
        case 'M0_Speed_TF_SF'
            % Null + Speed + TF + SF + Onset kernel
            X = [ones(n,1), B_speed, B_tf, D_sf, B_onset];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, onset_names];
        case 'Additive'
            % Full main effects: Speed + TF + SF + OR + Onset kernel
            X = [ones(n,1), B_speed, B_tf, D_sf, D_or, B_onset];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, or_names, onset_names];
        case 'FullInteraction'
            % Speed x TF interaction (vectorised via kron-like expansion)
            n_inter_st = n_speed_b * n_tf_b;
            B_inter_st = zeros(n, n_inter_st);
            inter_st_names = cell(1, n_inter_st);
            col = 0;
            for si = 1:n_speed_b
                for ti = 1:n_tf_b
                    col = col + 1;
                    B_inter_st(:, col) = B_speed(:, si) .* B_tf(:, ti);
                    inter_st_names{col} = sprintf('Spd%d_x_TF%d', si, ti);
                end
            end
            
            % Speed x SF interaction
            n_sf_cols = size(D_sf, 2);
            B_inter_ssf = zeros(n, n_speed_b * n_sf_cols);
            inter_ssf_names = cell(1, n_speed_b * n_sf_cols);
            col = 0;
            for si = 1:n_speed_b
                for fi = 1:n_sf_cols
                    col = col + 1;
                    B_inter_ssf(:, col) = B_speed(:, si) .* D_sf(:, fi);
                    inter_ssf_names{col} = sprintf('Spd%d_x_%s', si, sf_names{fi});
                end
            end
            
            % Speed x Orientation interaction
            n_or_cols = size(D_or, 2);
            B_inter_sor = zeros(n, n_speed_b * n_or_cols);
            inter_sor_names = cell(1, n_speed_b * n_or_cols);
            col = 0;
            for si = 1:n_speed_b
                for oi = 1:n_or_cols
                    col = col + 1;
                    B_inter_sor(:, col) = B_speed(:, si) .* D_or(:, oi);
                    inter_sor_names{col} = sprintf('Spd%d_x_%s', si, or_names{oi});
                end
            end
            
            % TF x SF interaction
            B_inter_tf_sf = zeros(n, n_tf_b * n_sf_cols);
            inter_tf_sf_names = cell(1, n_tf_b * n_sf_cols);
            col = 0;
            for ti = 1:n_tf_b
                for fi = 1:n_sf_cols
                    col = col + 1;
                    B_inter_tf_sf(:, col) = B_tf(:, ti) .* D_sf(:, fi);
                    inter_tf_sf_names{col} = sprintf('TF%d_x_%s', ti, sf_names{fi});
                end
            end
            
            % TF x Orientation interaction
            B_inter_tf_or = zeros(n, n_tf_b * n_or_cols);
            inter_tf_or_names = cell(1, n_tf_b * n_or_cols);
            col = 0;
            for ti = 1:n_tf_b
                for oi = 1:n_or_cols
                    col = col + 1;
                    B_inter_tf_or(:, col) = B_tf(:, ti) .* D_or(:, oi);
                    inter_tf_or_names{col} = sprintf('TF%d_x_%s', ti, or_names{oi});
                end
            end
            
            % SF x Orientation interaction
            B_inter_sf_or = zeros(n, n_sf_cols * n_or_cols);
            inter_sf_or_names = cell(1, n_sf_cols * n_or_cols);
            col = 0;
            for fi = 1:n_sf_cols
                for oi = 1:n_or_cols
                    col = col + 1;
                    B_inter_sf_or(:, col) = D_sf(:, fi) .* D_or(:, oi);
                    inter_sf_or_names{col} = sprintf('%s_x_%s', sf_names{fi}, or_names{oi});
                end
            end
            
            % Full interaction model with onset kernel
            X = [ones(n,1), B_speed, B_tf, D_sf, D_or, B_onset, ...
                 B_inter_st, B_inter_ssf, B_inter_sor, ...
                 B_inter_tf_sf, B_inter_tf_or, B_inter_sf_or];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, or_names, onset_names, ...
                         inter_st_names, inter_ssf_names, inter_sor_names, ...
                         inter_tf_sf_names, inter_tf_or_names, inter_sf_or_names];
        case 'Additive_no_Speed'
            % Additive model without Speed: TF + SF + OR + Onset
            X = [ones(n,1), B_tf, D_sf, D_or, B_onset];
            col_names = [{'Intercept'}, tf_names_list, sf_names, or_names, onset_names];
        case 'Additive_no_TF'
            % Additive model without TF: Speed + SF + OR + Onset
            X = [ones(n,1), B_speed, D_sf, D_or, B_onset];
            col_names = [{'Intercept'}, spd_names, sf_names, or_names, onset_names];
        case 'Additive_no_SF'
            % Additive model without SF: Speed + TF + OR + Onset
            X = [ones(n,1), B_speed, B_tf, D_or, B_onset];
            col_names = [{'Intercept'}, spd_names, tf_names_list, or_names, onset_names];
        case 'Additive_no_OR'
            % Additive model without OR: Speed + TF + SF + Onset
            X = [ones(n,1), B_speed, B_tf, D_sf, B_onset];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, onset_names];
        otherwise
            error('Unknown model label: %s', model_label);
    end
    
    % Remove zero-variance columns (except intercept) - only during training, not prediction
    if ~is_prediction
        keep_cols = true(1, size(X, 2));
        col_vars = var(X(:, 2:end));
        keep_cols(2:end) = col_vars >= 1e-12;
        if ~all(keep_cols)
            X = X(:, keep_cols);
            col_names = col_names(keep_cols);
        end
    end
end


function result = fit_poisson_glm(X, y, offset, lambda_ridge)
%FIT_POISSON_GLM Fit Poisson GLM with log link via IRLS
    if nargin < 4, lambda_ridge = 0; end
    
    [n, p] = size(X);
    
    lambda_min = 1e-6;
    lambda_eff = max(lambda_ridge, lambda_min);
    
    beta = zeros(p, 1);
    beta(1) = log(max(mean(y), 0.1));
    
    max_iter = 100;
    tol = 1e-8;
    
    ws = warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    
    ridge_mat = lambda_eff * eye(p);
    ridge_mat(1,1) = 0;  % don't penalise intercept
    
    for iter = 1:max_iter
        eta = X * beta + offset;
        eta = max(min(eta, 20), -20);
        mu = exp(eta);
        mu = max(mu, 1e-10);
        
        W = mu;
        z = eta + (y - mu) ./ mu - offset;
        
        XtWX = X' * bsxfun(@times, X, W) + ridge_mat;
        XtWz = X' * (W .* z);
        
        beta_new = XtWX \ XtWz;
        
        if any(~isfinite(beta_new))
            break;
        end
        
        if max(abs(beta_new - beta)) < tol
            beta = beta_new;
            break;
        end
        beta = beta_new;
    end
    
    warning(ws);
    
    % Final predictions
    eta = X * beta + offset;
    eta = max(min(eta, 20), -20);
    mu = exp(eta);
    mu = max(mu, 1e-10);
    
    % Log-likelihood
    ll = sum(y .* log(mu) - mu - gammaln(y + 1));
    
    % Deviance
    y_pos = max(y, 1e-10);
    dev = 2 * sum(y .* log(y_pos ./ mu) - (y - mu));
    
    % Standard errors (Cholesky is ~3x faster than pinv for SPD matrices)
    W = mu;
    XtWX = X' * bsxfun(@times, X, W) + ridge_mat;
    se = zeros(p, 1);
    try
        R_chol = chol(XtWX);
        inv_R = R_chol \ eye(p);
        se = sqrt(sum(inv_R.^2, 2));
    catch
        try
            cov_beta = pinv(XtWX);
            se = sqrt(abs(diag(cov_beta)));
        catch
            % Leave SE as zeros
        end
    end
    
    result.beta = beta;
    result.se = se;
    result.log_likelihood = ll;
    result.aic = -2*ll + 2*p;
    result.bic = -2*ll + p*log(n);
    result.deviance = dev;
    result.n_params = p;
    result.n_obs = n;
    % predicted_count = mu = E[Y] = expected spike count per bin
    % To get firing rate lambda(t), divide by exposure time: lambda = mu / t
    % Since offset = log(t), we have: mu = exp(X*beta + offset) = t * exp(X*beta)
    % Therefore: lambda(t) = exp(X*beta) = mu / t
    result.predicted_count = mu;
    result.pearson_residuals = (y - mu) ./ sqrt(max(mu, 1e-10));
    result.dispersion = dev / max(n - p, 1);
    result.converged = (iter < max_iter) && all(isfinite(beta));
    result.lambda_ridge = lambda_ridge;
end


function [cv_ll, cv_bits_per_spike, fold_ll, cv_predicted_count] = cross_validate_glm(X, y, offset, fold_ids, lambda_ridge)
%CROSS_VALIDATE_GLM K-fold cross-validated log-likelihood for Poisson GLM
%  cv_predicted_count returns mu = E[Y] (expected spike count), not firing rate.
%  To get firing rate lambda(t), divide by exposure time: lambda = mu / t
    if nargin < 5, lambda_ridge = 0; end
    
    folds = unique(fold_ids);
    K = length(folds);
    
    total_ll = 0;
    total_spikes = 0;
    fold_ll = zeros(K, 1);
    cv_predicted_count = nan(size(y));
    
    for fi = 1:K
        test_idx = fold_ids == folds(fi);
        train_idx = ~test_idx;
        
        res = fit_poisson_glm(X(train_idx,:), y(train_idx), offset(train_idx), lambda_ridge);
        
        eta_test = X(test_idx,:) * res.beta + offset(test_idx);
        eta_test = max(min(eta_test, 20), -20);
        mu_test = exp(eta_test);
        mu_test = max(mu_test, 1e-10);
        
        cv_predicted_count(test_idx) = mu_test;
        
        ll_test = sum(y(test_idx) .* log(mu_test) - mu_test - gammaln(y(test_idx) + 1));
        fold_ll(fi) = ll_test;
        total_ll = total_ll + ll_test;
        total_spikes = total_spikes + sum(y(test_idx));
    end
    
    cv_ll = total_ll;
    if total_spikes > 0
        cv_bits_per_spike = (total_ll / total_spikes) / log(2);
    else
        cv_bits_per_spike = NaN;
    end
end


function [selected_vars, selection_history, final_cv_bps, null_cv_bps] = forward_select_model(...
    B_speed, B_tf, B_onset, sf_v, or_v, y, offset, fold_ids, threshold, ~, ...
    sf_ref_levels, or_ref_levels)
%FORWARD_SELECT_MODEL Hierarchical forward selection for GLM
%
%   Two-phase selection process:
%   Phase 1: Select main effects (Speed, TF, SF, OR) one at a time
%   Phase 2: Select interactions, but ONLY if both parent main effects
%            were selected in Phase 1
%
%   This ensures proper model hierarchy and interpretability.
%
%   Inputs:
%       B_speed, B_tf, B_onset - Pre-computed basis matrices
%       sf_v, or_v             - Categorical variable vectors
%       y                      - Spike counts
%       offset                 - Log of bin width
%       fold_ids               - CV fold assignments
%       threshold              - Minimum Δ bps for inclusion
%       ~                      - (unused, kept for API compatibility)
%       sf_ref_levels, or_ref_levels - Reference levels for dummy coding
%
%   Outputs:
%       selected_vars    - Cell array of selected variable names
%       selection_history - Struct array with round details
%       final_cv_bps     - CV bits/spike of final selected model
%       null_cv_bps      - CV bits/spike of null model (baseline)

    % Define main effects and interactions separately
    main_effects = {'Speed', 'TF', 'SF', 'OR'};
    interactions = {'Speed_x_TF', 'Speed_x_SF', 'Speed_x_OR', 'TF_x_SF', 'TF_x_OR', 'SF_x_OR'};
    
    % Map interactions to their parent main effects
    interaction_parents = containers.Map();
    interaction_parents('Speed_x_TF') = {'Speed', 'TF'};
    interaction_parents('Speed_x_SF') = {'Speed', 'SF'};
    interaction_parents('Speed_x_OR') = {'Speed', 'OR'};
    interaction_parents('TF_x_SF')    = {'TF', 'SF'};
    interaction_parents('TF_x_OR')    = {'TF', 'OR'};
    interaction_parents('SF_x_OR')    = {'SF', 'OR'};

    % --- Fit Null model (intercept + onset kernel) ---
    [X_null, ~] = assemble_design_matrix_selected(B_speed, B_tf, B_onset, ...
        sf_v, or_v, {}, sf_ref_levels, or_ref_levels);
    [~, null_cv_bps] = cross_validate_glm(X_null, y, offset, fold_ids, 0);
    
    selected_vars = {};
    current_cv_bps = null_cv_bps;
    selection_history = struct('round', {}, 'phase', {}, 'tested', {}, 'delta_bps', {}, ...
        'best_candidate', {}, 'added', {}, 'cv_bps_after', {});
    
    round_num = 0;
    
    % =====================================================================
    % PHASE 1: Main effects only
    % =====================================================================
    fprintf('    --- Phase 1: Main Effects ---\n');
    remaining_main = main_effects;
    
    while ~isempty(remaining_main)
        round_num = round_num + 1;
        
        best_delta = -Inf;
        best_candidate = '';
        tested_results = struct('candidate', {}, 'delta_bps', {}, 'cv_bps', {});
        
        for ci = 1:length(remaining_main)
            candidate = remaining_main{ci};
            
            % Build model with current selected + this candidate
            test_vars = [selected_vars, {candidate}];
            [X_test, ~] = assemble_design_matrix_selected(B_speed, B_tf, B_onset, ...
                sf_v, or_v, test_vars, sf_ref_levels, or_ref_levels);
            
            % Check if design matrix is valid
            if size(X_test, 2) >= length(y)
                tested_results(ci).candidate = candidate;
                tested_results(ci).delta_bps = -Inf;
                tested_results(ci).cv_bps = -Inf;
                continue;
            end
            
            % Cross-validate
            [~, test_cv_bps] = cross_validate_glm(X_test, y, offset, fold_ids, 0);
            delta = test_cv_bps - current_cv_bps;
            
            tested_results(ci).candidate = candidate;
            tested_results(ci).delta_bps = delta;
            tested_results(ci).cv_bps = test_cv_bps;
            
            if delta > best_delta
                best_delta = delta;
                best_candidate = candidate;
            end
        end
        
        % Record this round
        selection_history(round_num).round = round_num;
        selection_history(round_num).phase = 1;
        selection_history(round_num).tested = tested_results;
        selection_history(round_num).best_candidate = best_candidate;
        selection_history(round_num).delta_bps = best_delta;
        
        % Diagnostic logging
        fprintf('    Round %d (main): ', round_num);
        for ti = 1:length(tested_results)
            fprintf('%s=%.4f ', tested_results(ti).candidate, tested_results(ti).delta_bps);
        end
        fprintf('\n');
        fprintf('    Best: %s (Δbps=%.4f), threshold=%.4f\n', best_candidate, best_delta, threshold);
        
        % Add best candidate if it exceeds threshold
        if best_delta > threshold && ~isempty(best_candidate)
            fprintf('    -> ADDED %s\n', best_candidate);
            selected_vars{end+1} = best_candidate; %#ok<AGROW>
            remaining_main = setdiff(remaining_main, {best_candidate});
            % Update current_cv_bps
            for ti = 1:length(tested_results)
                if strcmp(tested_results(ti).candidate, best_candidate)
                    current_cv_bps = tested_results(ti).cv_bps;
                    break;
                end
            end
            selection_history(round_num).added = true;
            selection_history(round_num).cv_bps_after = current_cv_bps;
        else
            selection_history(round_num).added = false;
            selection_history(round_num).cv_bps_after = current_cv_bps;
            break;  % Stop Phase 1: no main effect improved enough
        end
    end
    
    % =====================================================================
    % PHASE 2: Interactions (only if both parents were selected)
    % =====================================================================
    fprintf('    --- Phase 2: Interactions (eligible based on selected main effects) ---\n');
    
    % Determine eligible interactions
    eligible_interactions = {};
    for ii = 1:length(interactions)
        int_name = interactions{ii};
        parents = interaction_parents(int_name);
        % Check if BOTH parents are in selected_vars
        if all(ismember(parents, selected_vars))
            eligible_interactions{end+1} = int_name; %#ok<AGROW>
        end
    end
    
    if isempty(eligible_interactions)
        fprintf('    No eligible interactions (need both parent main effects selected)\n');
    else
        fprintf('    Eligible interactions: %s\n', strjoin(eligible_interactions, ', '));
    end
    
    remaining_interactions = eligible_interactions;
    
    while ~isempty(remaining_interactions)
        round_num = round_num + 1;
        
        best_delta = -Inf;
        best_candidate = '';
        tested_results = struct('candidate', {}, 'delta_bps', {}, 'cv_bps', {});
        
        for ci = 1:length(remaining_interactions)
            candidate = remaining_interactions{ci};
            
            % Build model with current selected + this interaction
            test_vars = [selected_vars, {candidate}];
            [X_test, ~] = assemble_design_matrix_selected(B_speed, B_tf, B_onset, ...
                sf_v, or_v, test_vars, sf_ref_levels, or_ref_levels);
            
            % Check if design matrix is valid
            if size(X_test, 2) >= length(y)
                tested_results(ci).candidate = candidate;
                tested_results(ci).delta_bps = -Inf;
                tested_results(ci).cv_bps = -Inf;
                continue;
            end
            
            % Cross-validate
            [~, test_cv_bps] = cross_validate_glm(X_test, y, offset, fold_ids, 0);
            delta = test_cv_bps - current_cv_bps;
            
            tested_results(ci).candidate = candidate;
            tested_results(ci).delta_bps = delta;
            tested_results(ci).cv_bps = test_cv_bps;
            
            if delta > best_delta
                best_delta = delta;
                best_candidate = candidate;
            end
        end
        
        % Record this round
        selection_history(round_num).round = round_num;
        selection_history(round_num).phase = 2;
        selection_history(round_num).tested = tested_results;
        selection_history(round_num).best_candidate = best_candidate;
        selection_history(round_num).delta_bps = best_delta;
        
        % Diagnostic logging
        fprintf('    Round %d (interaction): ', round_num);
        for ti = 1:length(tested_results)
            fprintf('%s=%.4f ', tested_results(ti).candidate, tested_results(ti).delta_bps);
        end
        fprintf('\n');
        fprintf('    Best: %s (Δbps=%.4f), threshold=%.4f\n', best_candidate, best_delta, threshold);
        
        % Warn about extreme negative deltas (numerical instability)
        if best_delta < -1
            fprintf('    WARNING: Large negative delta (%.2f) suggests numerical instability\n', best_delta);
        end
        
        % Add best candidate if it exceeds threshold
        if best_delta > threshold && ~isempty(best_candidate)
            fprintf('    -> ADDED %s\n', best_candidate);
            selected_vars{end+1} = best_candidate; %#ok<AGROW>
            remaining_interactions = setdiff(remaining_interactions, {best_candidate});
            % Update current_cv_bps
            for ti = 1:length(tested_results)
                if strcmp(tested_results(ti).candidate, best_candidate)
                    current_cv_bps = tested_results(ti).cv_bps;
                    break;
                end
            end
            selection_history(round_num).added = true;
            selection_history(round_num).cv_bps_after = current_cv_bps;
        else
            selection_history(round_num).added = false;
            selection_history(round_num).cv_bps_after = current_cv_bps;
            break;  % Stop Phase 2: no interaction improved enough
        end
    end
    
    final_cv_bps = current_cv_bps;
    
    % Summary
    selected_main = intersect(selected_vars, main_effects, 'stable');
    selected_int = intersect(selected_vars, interactions, 'stable');
    fprintf('    FINAL: Main effects: [%s], Interactions: [%s]\n', ...
        strjoin(selected_main, ', '), strjoin(selected_int, ', '));
end


function [X, col_names] = assemble_design_matrix_selected(B_speed, B_tf, B_onset, ...
    sf_vals, or_vals, selected_vars, sf_ref_levels, or_ref_levels)
%ASSEMBLE_DESIGN_MATRIX_SELECTED Build design matrix for forward-selected model
%
%   Creates design matrix containing:
%   - Intercept (always)
%   - Onset kernel (always)
%   - Only the variables specified in selected_vars
%
%   For interactions, automatically includes parent main effect bases in the
%   interaction columns (but NOT as separate main effects unless explicitly selected).
%
%   Inputs:
%       B_speed, B_tf, B_onset - Pre-computed basis matrices
%       sf_vals, or_vals       - Categorical variable vectors
%       selected_vars          - Cell array of selected variable names
%       sf_ref_levels, or_ref_levels - Reference levels for dummy coding
%
%   Outputs:
%       X         - Design matrix [n_obs x n_params]
%       col_names - Cell array of column names

    n = size(B_speed, 1);
    n_speed_bases = size(B_speed, 2);
    n_tf_bases = size(B_tf, 2);
    n_onset_bases = size(B_onset, 2);
    
    % Determine SF levels for dummy coding
    % When sf_ref_levels is provided (prediction mode), use it as the complete set of levels
    % and exclude the first level (reference). This matches assemble_design_matrix behavior.
    if nargin >= 7 && ~isempty(sf_ref_levels)
        unique_sf = sf_ref_levels(:)';
        sf_levels = unique_sf(2:end);  % Exclude first (reference) level
    else
        sf_unique = unique(sf_vals(~isnan(sf_vals)));
        sf_levels = sf_unique(2:end);  % First level is reference
    end
    
    % Determine OR levels for dummy coding (same logic as SF)
    if nargin >= 8 && ~isempty(or_ref_levels)
        unique_or = or_ref_levels(:)';
        or_levels = unique_or(2:end);  % Exclude first (reference) level
    else
        or_unique = unique(or_vals(~isnan(or_vals)));
        or_levels = or_unique(2:end);  % First level is reference
    end
    
    % Start with intercept and onset kernel
    X = ones(n, 1);
    col_names = {'Intercept'};
    
    % Add onset kernel bases
    X = [X, B_onset];
    for bi = 1:n_onset_bases
        col_names{end+1} = sprintf('Onset_%d', bi); %#ok<AGROW>
    end
    
    % Check which variables are selected
    has_spprofileeed = ismember('Speed', selected_vars);
    has_tf = ismember('TF', selected_vars);
    has_sf = ismember('SF', selected_vars);
    has_or = ismember('OR', selected_vars);
    
    has_spprofileeed_x_tf = ismember('Speed_x_TF', selected_vars);
    has_spprofileeed_x_sf = ismember('Speed_x_SF', selected_vars);
    has_spprofileeed_x_or = ismember('Speed_x_OR', selected_vars);
    has_tf_x_sf = ismember('TF_x_SF', selected_vars);
    has_tf_x_or = ismember('TF_x_OR', selected_vars);
    has_sf_x_or = ismember('SF_x_OR', selected_vars);
    
    % --- Add main effects ---
    if has_spprofileeed
        X = [X, B_speed];
        for bi = 1:n_speed_bases
            col_names{end+1} = sprintf('Speed_%d', bi); %#ok<AGROW>
        end
    end
    
    if has_tf
        X = [X, B_tf];
        for bi = 1:n_tf_bases
            col_names{end+1} = sprintf('TF_%d', bi); %#ok<AGROW>
        end
    end
    
    if has_sf
        for li = 1:length(sf_levels)
            sf_dummy = double(sf_vals == sf_levels(li));
            sf_dummy(isnan(sf_vals)) = 0;
            X = [X, sf_dummy]; %#ok<AGROW>
            col_names{end+1} = sprintf('SF_%.4f', sf_levels(li)); %#ok<AGROW>
        end
    end
    
    if has_or
        for li = 1:length(or_levels)
            or_dummy = double(or_vals == or_levels(li));
            or_dummy(isnan(or_vals)) = 0;
            X = [X, or_dummy]; %#ok<AGROW>
            col_names{end+1} = sprintf('OR_%.3f', or_levels(li)); %#ok<AGROW>
        end
    end
    
    % --- Add interactions ---
    % For interactions, we create the interaction columns directly
    % (products of bases/dummies), without requiring main effects to be selected
    
    if has_spprofileeed_x_tf
        for si = 1:n_speed_bases
            for ti = 1:n_tf_bases
                X = [X, B_speed(:,si) .* B_tf(:,ti)]; %#ok<AGROW>
                col_names{end+1} = sprintf('Spd%d_x_TF%d', si, ti); %#ok<AGROW>
            end
        end
    end
    
    if has_spprofileeed_x_sf
        for si = 1:n_speed_bases
            for li = 1:length(sf_levels)
                sf_dummy = double(sf_vals == sf_levels(li));
                sf_dummy(isnan(sf_vals)) = 0;
                X = [X, B_speed(:,si) .* sf_dummy]; %#ok<AGROW>
                col_names{end+1} = sprintf('Spd%d_x_SF%.4f', si, sf_levels(li)); %#ok<AGROW>
            end
        end
    end
    
    if has_spprofileeed_x_or
        for si = 1:n_speed_bases
            for li = 1:length(or_levels)
                or_dummy = double(or_vals == or_levels(li));
                or_dummy(isnan(or_vals)) = 0;
                X = [X, B_speed(:,si) .* or_dummy]; %#ok<AGROW>
                col_names{end+1} = sprintf('Spd%d_x_OR%.3f', si, or_levels(li)); %#ok<AGROW>
            end
        end
    end
    
    if has_tf_x_sf
        for ti = 1:n_tf_bases
            for li = 1:length(sf_levels)
                sf_dummy = double(sf_vals == sf_levels(li));
                sf_dummy(isnan(sf_vals)) = 0;
                X = [X, B_tf(:,ti) .* sf_dummy]; %#ok<AGROW>
                col_names{end+1} = sprintf('TF%d_x_SF%.4f', ti, sf_levels(li)); %#ok<AGROW>
            end
        end
    end
    
    if has_tf_x_or
        for ti = 1:n_tf_bases
            for li = 1:length(or_levels)
                or_dummy = double(or_vals == or_levels(li));
                or_dummy(isnan(or_vals)) = 0;
                X = [X, B_tf(:,ti) .* or_dummy]; %#ok<AGROW>
                col_names{end+1} = sprintf('TF%d_x_OR%.3f', ti, or_levels(li)); %#ok<AGROW>
            end
        end
    end
    
    if has_sf_x_or
        for li_sf = 1:length(sf_levels)
            sf_dummy = double(sf_vals == sf_levels(li_sf));
            sf_dummy(isnan(sf_vals)) = 0;
            for li_or = 1:length(or_levels)
                or_dummy = double(or_vals == or_levels(li_or));
                or_dummy(isnan(or_vals)) = 0;
                X = [X, sf_dummy .* or_dummy]; %#ok<AGROW>
                col_names{end+1} = sprintf('SF%.4f_x_OR%.3f', sf_levels(li_sf), or_levels(li_or)); %#ok<AGROW>
            end
        end
    end
    
    % Remove zero-variance columns
    col_vars = var(X, 0, 1);
    keep_cols = col_vars > 1e-10 | (1:size(X,2)) == 1;  % Always keep intercept
    X = X(:, keep_cols);
    col_names = col_names(keep_cols);
end


function [grp_mean, grp_sem, grp_clr, grp_lbl, n_grps, grp_betas, grp_cn] = compute_grouped_params(b_all, se_all, cn_all)
%COMPUTE_GROUPED_PARAMS Group coefficients by feature type for swarm plot
%  Returns per-group individual betas and column names for basis labeling
    n_coef = length(b_all);
    
    grp_tags = cell(n_coef, 1);
    for ki = 1:n_coef
        cn = cn_all{ki};
        if strcmp(cn, 'Intercept')
            grp_tags{ki} = 'Intercept';
        elseif contains(cn, '_x_TF')
            grp_tags{ki} = 'Spd x TF';
        elseif contains(cn, '_x_SF')
            grp_tags{ki} = 'Spd x SF';
        elseif contains(cn, '_x_OR')
            grp_tags{ki} = 'Spd x OR';
        elseif startsWith(cn, 'Speed_')
            grp_tags{ki} = 'Speed';
        elseif startsWith(cn, 'TF_')
            grp_tags{ki} = 'TF';
        elseif startsWith(cn, 'SF_')
            grp_tags{ki} = 'SF';
        elseif startsWith(cn, 'OR_')
            grp_tags{ki} = 'OR';
        elseif startsWith(cn, 'Time_')
            grp_tags{ki} = 'Time';
        else
            grp_tags{ki} = 'Other';
        end
    end
    
    grp_order = {'Intercept','Speed','TF','SF','OR','Time', ...
                 'Spd x TF','Spd x SF','Spd x OR','Other'};
    % Colors matched to model barplot color scheme (line 2069-2072)
    grp_colors = [0.5 0.5 0.5;   ... % Intercept (gray)
                  0.17 0.63 0.17;... % Speed (green - matches color_speed)
                  1.0 0.50 0.05; ... % TF (orange - matches color_tf)
                  0.95 0.85 0.10;... % SF (yellow - matches color_sf)
                  0.84 0.15 0.16;... % OR (red - matches color_or)
                  0.3 0.8 0.8;   ... % Time (cyan)
                  0.6 0.35 0.05; ... % Spd x TF (darker orange-brown)
                  0.55 0.50 0.05;... % Spd x SF (olive)
                  0.50 0.10 0.10;... % Spd x OR (dark red)
                  0.7 0.7 0.7];      % Other (gray)
    
    grp_mean = []; grp_sem = []; grp_clr = []; grp_lbl = {};
    grp_betas = {}; grp_cn = {};
    n_grps = 0;
    for gi = 1:length(grp_order)
        mask = strcmp(grp_tags, grp_order{gi});
        if ~any(mask), continue; end
        n_grps = n_grps + 1;
        betas_g = b_all(mask);
        grp_mean(n_grps) = mean(betas_g);                      %#ok<AGROW>
        grp_sem(n_grps)  = std(betas_g) / sqrt(sum(mask));     %#ok<AGROW>
        grp_clr(n_grps,:) = grp_colors(gi,:);                  %#ok<AGROW>
        grp_lbl{n_grps}  = grp_order{gi};                      %#ok<AGROW>
        grp_betas{n_grps} = betas_g(:)';                       %#ok<AGROW>
        grp_cn{n_grps}   = cn_all(mask);                       %#ok<AGROW>
    end
end


function lbl = extract_basis_label(col_name)
%EXTRACT_BASIS_LABEL Extract basis index label from column name
%   'Speed_3' -> '3'
%   'TF_2' -> '2'
%   'Spd2_x_TF4' -> '2-4'
%   'Spd1_x_SF_2' -> '1-2'
%   'Intercept' -> ''
    lbl = '';
    
    if strcmp(col_name, 'Intercept')
        lbl = '';
        return;
    end
    
    % Check for interaction terms first (contain '_x_')
    if contains(col_name, '_x_')
        % Parse interaction: e.g., 'Spd2_x_TF4' or 'Spd1_x_SF_2'
        parts = strsplit(col_name, '_x_');
        if length(parts) == 2
            % First part: extract number from Spd#
            num1 = regexp(parts{1}, '\d+', 'match');
            % Second part: extract number (TF#, SF_#, OR_#)
            num2 = regexp(parts{2}, '\d+', 'match');
            if ~isempty(num1) && ~isempty(num2)
                lbl = sprintf('%s-%s', num1{end}, num2{end});
            end
        end
    else
        % Main effect: Speed_#, TF_#, SF_#, OR_#, Time_#
        nums = regexp(col_name, '\d+', 'match');
        if ~isempty(nums)
            lbl = nums{end};
        end
    end
end


function str = format_radians(val)
%FORMAT_RADIANS Format angle value as a fraction of pi for display
%   Converts common angles to nice π notation (e.g., -π/4, 0, π/4, π/2)
    
    tol = 1e-6;
    if abs(val) < tol
        str = '0';
    elseif abs(val - pi/2) < tol
        str = 'π/2';
    elseif abs(val + pi/2) < tol
        str = '-π/2';
    elseif abs(val - pi/4) < tol
        str = 'π/4';
    elseif abs(val + pi/4) < tol
        str = '-π/4';
    elseif abs(val - pi) < tol
        str = 'π';
    elseif abs(val + pi) < tol
        str = '-π';
    elseif abs(val - 3*pi/4) < tol
        str = '3π/4';
    elseif abs(val + 3*pi/4) < tol
        str = '-3π/4';
    else
        % General case: show decimal radians
        str = sprintf('%.2f', val);
    end
end


function X = build_design_matrix_from_colnames(col_names, B_speed, B_tf, B_onset, sf_val, or_val, sf_all_levels, or_all_levels)
%BUILD_DESIGN_MATRIX_FROM_COLNAMES Build design matrix matching stored column names
%   This function guarantees the prediction design matrix exactly matches the
%   training design matrix by building columns based on the stored column names.
%
%   Inputs:
%       col_names      - Cell array of column names from training
%       B_speed        - Speed basis matrix [n x n_speed_bases]
%       B_tf           - TF basis matrix [n x n_tf_bases]
%       B_onset        - Onset basis matrix [n x n_onset_bases]
%       sf_val         - SF value(s) for this prediction [n x 1] or scalar
%       or_val         - OR value(s) for this prediction [n x 1] or scalar
%       sf_all_levels  - All SF levels used during training
%       or_all_levels  - All OR levels used during training
%
%   Output:
%       X - Design matrix with columns matching col_names

    n = size(B_speed, 1);
    n_cols = length(col_names);
    X = zeros(n, n_cols);
    
    % Ensure sf_val and or_val are column vectors
    if isscalar(sf_val), sf_val = repmat(sf_val, n, 1); end
    if isscalar(or_val), or_val = repmat(or_val, n, 1); end
    
    % SF levels for dummy coding (exclude reference = first level)
    sf_dummy_levels = sf_all_levels(2:end);
    or_dummy_levels = or_all_levels(2:end);
    
    for ci = 1:n_cols
        cn = col_names{ci};
        
        if strcmp(cn, 'Intercept')
            X(:, ci) = 1;
            
        elseif startsWith(cn, 'Onset_')
            % Parse onset basis index: 'Onset_3' -> 3
            idx = sscanf(cn, 'Onset_%d');
            if ~isempty(idx) && idx <= size(B_onset, 2)
                X(:, ci) = B_onset(:, idx);
            end
            
        elseif startsWith(cn, 'Speed_')
            % Parse speed basis index: 'Speed_2' -> 2
            idx = sscanf(cn, 'Speed_%d');
            if ~isempty(idx) && idx <= size(B_speed, 2)
                X(:, ci) = B_speed(:, idx);
            end
            
        elseif startsWith(cn, 'TF_')
            % Parse TF basis index: 'TF_4' -> 4
            idx = sscanf(cn, 'TF_%d');
            if ~isempty(idx) && idx <= size(B_tf, 2)
                X(:, ci) = B_tf(:, idx);
            end
            
        elseif startsWith(cn, 'SF_')
            % SF dummy: 'SF_2' means second non-reference level, or 'SF_0.0060' format
            % Try numeric format first
            sf_level = sscanf(cn, 'SF_%f');
            if ~isempty(sf_level)
                X(:, ci) = double(abs(sf_val - sf_level) < 1e-6);
            else
                % Try index format: 'SF_2' -> second dummy level
                idx = sscanf(cn, 'SF_%d');
                if ~isempty(idx) && idx <= length(sf_dummy_levels)
                    X(:, ci) = double(abs(sf_val - sf_dummy_levels(idx)) < 1e-6);
                end
            end
            
        elseif startsWith(cn, 'OR_')
            % OR dummy: 'OR_2' or 'OR_0.785' format
            or_level = sscanf(cn, 'OR_%f');
            if ~isempty(or_level)
                X(:, ci) = double(abs(or_val - or_level) < 0.01);
            else
                idx = sscanf(cn, 'OR_%d');
                if ~isempty(idx) && idx <= length(or_dummy_levels)
                    X(:, ci) = double(abs(or_val - or_dummy_levels(idx)) < 0.01);
                end
            end
            
        elseif contains(cn, '_x_')
            % Interaction terms: 'Spd2_x_TF3', 'Spd1_x_SF0.0060', etc.
            parts = strsplit(cn, '_x_');
            if length(parts) == 2
                % Get first component
                comp1 = get_component_column(parts{1}, B_speed, B_tf, sf_val, or_val, sf_dummy_levels, or_dummy_levels);
                % Get second component
                comp2 = get_component_column(parts{2}, B_speed, B_tf, sf_val, or_val, sf_dummy_levels, or_dummy_levels);
                X(:, ci) = comp1 .* comp2;
            end
        end
    end
end


function col = get_component_column(name, B_speed, B_tf, sf_val, or_val, sf_levels, or_levels)
%GET_COMPONENT_COLUMN Get a single column for interaction term component
    n = size(B_speed, 1);
    col = zeros(n, 1);
    
    if startsWith(name, 'Spd')
        idx = sscanf(name, 'Spd%d');
        if ~isempty(idx) && idx <= size(B_speed, 2)
            col = B_speed(:, idx);
        end
    elseif startsWith(name, 'TF')
        idx = sscanf(name, 'TF%d');
        if ~isempty(idx) && idx <= size(B_tf, 2)
            col = B_tf(:, idx);
        end
    elseif startsWith(name, 'SF')
        sf_level = sscanf(name, 'SF%f');
        if ~isempty(sf_level)
            col = double(abs(sf_val - sf_level) < 1e-6);
        end
    elseif startsWith(name, 'OR')
        or_level = sscanf(name, 'OR%f');
        if ~isempty(or_level)
            col = double(abs(or_val - or_level) < 0.01);
        end
    end
end
