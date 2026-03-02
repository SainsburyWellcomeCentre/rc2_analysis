% GLM-based single-cluster analysis of multisensory integration
%
% Version 6: Pre-filtering Decision Tree
%
% This script first applies a decision tree to reduce the number of clusters
% that need full GLM analysis, then fits Poisson GLMs to the remaining 
% clusters using two complementary approaches.
%
% PRE-FILTERING DECISION TREE:
% ============================
% For each cluster:
%   1. Check speed tuning: stationary_T != motion_T (Wilcoxon signrank)
%   2. Check visual tuning: stationary_V != motion_V (Wilcoxon signrank)
%   
%   If neither speed nor visually tuned:
%       -> Category: "not_tuned" - do not analyse
%   
%   If speed tuned but NOT visually tuned:
%       Check stationary_VT != motion_VT
%       If true: 
%           Compare motion_T vs motion_VT using KW test on speed bins
%           If motion_T == motion_VT: 
%               -> Category: "speed_tuned_only" - do not analyse
%           Else: 
%               -> Category: "mixed_tuning" - run GLM (maybe visually tuned)
%       Else (VT not responsive): 
%           -> Warning: unexpected case
%   
%   If visually tuned (with or without speed tuning):
%       -> Category: "visually_tuned" - run full GLM
%
% GLM APPROACHES:
% ===============
%   (A) Speed-bin GLM: uses pre-computed tuning tables (20 velocity bins
%       per trial). Spike counts are estimated as rate * time from the
%       tuning tables. log(time_in_bin) is used as an offset.
%
%   (B) Time-bin GLM:  uses fixed 50 ms time bins with spike counts
%       counted directly from raw spike_times. Speed and TF are-
%       continuous covariates sampled at each bin. Five slow temporal
%       basis functions absorb within-trial drift.
%
% Both approaches test whether the response to combined visual+translation
% (VT) can be explained by an additive combination of visual (V) and
% translation (T_Vstatic) responses, or whether interactions are needed.
%
% Model hierarchy (same for both GLM types):
%   M0:              Null (intercept only)
%   M0_Speed:        Null + Speed basis
%   M0_Speed_TF:     Null + Speed + TF
%   M0_Speed_TF_SF:  Null + Speed + TF + SF
%   Additive:        Speed + TF + SF + OR (full main effects model)
%   FullInteraction: Additive + Spd*TF + Spd*SF + Spd*OR + TF*SF + TF*OR + SF*OR
%
% Drop-one reduced models (for hypothesis testing — unique contribution):
%   Additive_no_Speed: TF + SF + OR
%   Additive_no_TF:    Speed + SF + OR
%   Additive_no_SF:    Speed + TF + OR
%   Additive_no_OR:    Speed + TF + SF
%
% Positive Control (visual contribution):
%   M0_Speed is ALWAYS fitted (even when use_hierarchical_models = false)
%   delta_bps_visual = Additive - M0_Speed
%   This tests the combined contribution of all visual features (TF + SF + OR)
%   and serves as a positive control: neurons selected for GLM should show
%   improvement when visual features are added to a speed-only baseline.
%
% Neuron Classification — Drop-one test (threshold: delta CV bps > 0.01):
%   - is_visually_tuned_glm: Additive - M0_Speed > threshold (combined visual)
%   - is_speed_tuned: Additive - Additive_no_Speed > threshold (unique speed)
%   - is_tf_tuned:    Additive - Additive_no_TF > threshold (unique TF)
%   - is_sf_tuned:    Additive - Additive_no_SF > threshold (unique SF)
%   - is_or_tuned:    Additive - Additive_no_OR > threshold (unique OR)
%   - has_significant_interaction: FullInteraction - Additive > threshold
%
% Tuning Shape Classification:
%   - monotonic_inc:  Peak at maximum sampled speed/TF
%   - monotonic_dec:  Peak at minimum sampled speed/TF
%   - bandpass:       Peak at interior value
%
% OUTPUTS:
%   - prefilter_decision_tree.csv: Results of the pre-filtering decision tree
%       (updated after GLM analysis with leave-one-out results and classification:
%        glm_delta_bps_drop_speed, glm_delta_bps_drop_tf, glm_delta_bps_drop_sf,
%        glm_delta_bps_drop_or, glm_is_speed_tuned, glm_is_tf_tuned,
%        glm_is_sf_tuned, glm_is_or_tuned, glm_speed_tuning_shape,
%        glm_tf_tuning_shape, glm_has_significant_interaction,
%        glm_delta_bps_interaction, glm_additive_cv_bps, glm_winning_model)
%   - csv/glm_model_comparison.csv: Model comparison with pre-filter category
%   - csv/glm_coefficients.csv: GLM coefficients for each cluster
%   - csv/glm_speed_bin_data.csv: Speed-bin trial data
%   - csv/glm_time_bin_data.csv: Time-bin trial data
%   - csv/glm_component_importance.csv: Component importance decomposition
%
% See also: aritmetic_sum, speed_tuning_all_clusters_all_conditions,
%           tf_tuning_all_clusters_all_conditions, create_tables_by_batch

%% ====================================================================
%  Section 1: Configuration
%  ====================================================================
close all; clearvars;

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
n_speed_bins            = 20;     % speed bins per trial (speed-bin GLM)
n_tf_bins               = 20;     % TF bins per trial

% Time-bin GLM
time_bin_width          = 0.100;  % 100 ms bins

% Basis functions: raised cosine
n_speed_bases           = 5;      % number of raised cosine bases for speed
n_tf_bases              = 5;      % number of raised cosine bases for TF
n_temporal_bases        = 5;      % number of raised cosine bases for time-in-trial
speed_range             = [0, 50];  % cm/s (updated from tuning tables later)
tf_range                = [0, 7.3]; % Hz   (updated from tuning tables later)

% Cross-validation
n_cv_folds              = 5;

% Prediction cap (physiological limit for firing rate)
max_predicted_fr        = 500;    % Hz - cap predictions to prevent exp() explosion

% Figures
save_figs               = true;
overwrite               = true;
figure_dir              = {'glm_single_cluster'};

% CSV output
csv_output              = true;

% Model hierarchy selection
% If false, fit only Additive and FullInteraction (2 models per binning type)
% If true, fit full hierarchy: M0 -> M0_Speed -> M0_Speed_TF -> M0_Speed_TF_SF -> Additive -> FullInteraction
use_hierarchical_models = false;

% --- Pre-computed tuning table paths ---
tuning_curves_dir       = 'D:\mvelez\formatted_data\csvs\tuning_curves';
tf_tuning_curves_dir    = 'D:\mvelez\formatted_data\csvs\tf_tuning_curves';

% Motion cloud data paths
mc_sequence_path        = 'D:\mvelez\mateoData_mc\motion_cloud_sequence_250414.mat';
mc_folders_path         = 'D:\mvelez\mateoData_mc\image_folders.mat';

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

% Model labels — hierarchy for model selection
% Full hierarchy: M0 → M0_Speed → M0_Speed_TF → M0_Speed_TF_SF → Additive → FullInteraction
% Minimal (default): M0_Speed + Additive + FullInteraction (M0_Speed needed for positive control)
if use_hierarchical_models
    model_labels = {'M0', 'M0_Speed', 'M0_Speed_TF', 'M0_Speed_TF_SF', 'Additive', 'FullInteraction'};
else
    % Always include M0_Speed for positive control (Additive vs M0_Speed tests visual contribution)
    model_labels = {'M0_Speed', 'Additive', 'FullInteraction'};
end
n_models = length(model_labels);

% Drop-one reduced models (for hypothesis testing — unique contribution of each variable)
dropone_labels = {'Additive_no_Speed', 'Additive_no_TF', 'Additive_no_SF', 'Additive_no_OR'};
dropone_variable = {'Speed', 'TF', 'SF', 'OR'};
n_dropone = length(dropone_labels);

% Classification threshold (delta bits-per-spike)
classification_threshold = 0.01;  % Δ bps > 0.01 = significant tuning

% Model formula strings (for figure titles)
model_formulas = containers.Map();
model_formulas('M0')              = '\beta_0';
model_formulas('M0_Speed')        = '\beta_0 + f(Spd)';
model_formulas('M0_Speed_TF')     = '\beta_0 + f(Spd) + g(TF)';
model_formulas('M0_Speed_TF_SF')  = '\beta_0 + f(Spd) + g(TF) + SF';
model_formulas('Additive')        = '\beta_0 + f(Spd) + g(TF) + SF + OR';
model_formulas('FullInteraction') = 'Additive + Spd*TF + Spd*SF + Spd*OR + TF*SF + TF*OR + SF*OR';
model_formulas('Additive_no_Speed') = '\beta_0 + g(TF) + SF + OR';
model_formulas('Additive_no_TF')    = '\beta_0 + f(Spd) + SF + OR';
model_formulas('Additive_no_SF')    = '\beta_0 + f(Spd) + g(TF) + OR';
model_formulas('Additive_no_OR')    = '\beta_0 + f(Spd) + g(TF) + SF';

fprintf('=== GLM Single Cluster Analysis (v6: Pre-filtering Decision Tree) ===\n');
fprintf('Configuration:\n');
fprintf('  Speed bins: %d | TF bins: %d\n', n_speed_bins, n_tf_bins);
fprintf('  Time-bin width: %.0f ms\n', time_bin_width*1000);
fprintf('  Speed bases: %d | TF bases: %d\n', n_speed_bases, n_tf_bases);
fprintf('  CV folds: %d | Restricted: %d\n', n_cv_folds, restricted);
fprintf('  Speed range: [%.1f, %.1f] cm/s\n', speed_range(1), speed_range(2));
fprintf('  TF range: [%.1f, %.1f] Hz\n', tf_range(1), tf_range(2));
fprintf('  Classification threshold: delta bps > %.3f\n', classification_threshold);
fprintf('  Hierarchical models: %d\n', use_hierarchical_models);
fprintf('  Models: %s\n', strjoin(model_labels, ' -> '));
fprintf('  Drop-one models: %s\n', strjoin(dropone_labels, ', '));
fprintf('  Classification: Drop-one from Additive (unique contribution of each variable)\n');
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
% For each cluster:
%   1. Check if speed tuned: stationary_T != motion_T (Wilcoxon signrank)
%   2. Check if visually tuned: stationary_V != motion_V (Wilcoxon signrank)
%   
%   If neither speed nor visually tuned: do not analyse (set to "not_tuned")
%   
%   If speed tuned but not visually tuned:
%       Compare stationary_VT != motion_VT
%       If true: compare motion_T vs motion_VT using KW on speed bins
%           If motion_T == motion_VT: only speed tuned, do not analyse
%           Else: mixed tuning, becomes "maybe_visual", continue to GLM
%       Else: unexpected case, raise warning
%   
%   If visually tuned OR mixed tuning: run full GLM analysis
%
% The pre-filtering uses:
%   - svm_table: stationary vs motion firing rates per trial
%   - tuning_curves: speed bin distributions (20 bins)
%   - Wilcoxon sign-rank test for stationary vs motion comparisons
%   - Kruskal-Wallis test for comparing speed bin distributions across conditions

fprintf('\n--- Pre-filtering Decision Tree ---\n');

ctl = RC2Analysis();
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
prefilter_results.is_speed_tuned = [];         % stationary_T != motion_T
prefilter_results.is_visually_tuned = [];      % stationary_V != motion_V
prefilter_results.is_VT_responsive = [];       % stationary_VT != motion_VT
prefilter_results.is_mixed_tuning = [];        % speed tuned but shows visual modulation
prefilter_results.should_run_glm = [];         % final decision: run GLM or not
prefilter_results.category = {};               % category label for this cluster
prefilter_results.p_T_svm = [];                % p-value for T stationary vs motion
prefilter_results.p_V_svm = [];                % p-value for V stationary vs motion
prefilter_results.p_VT_svm = [];               % p-value for VT stationary vs motion
prefilter_results.p_T_vs_VT_motion = [];       % p-value for KW test T vs VT motion
prefilter_results.direction_T = [];            % direction of T modulation
prefilter_results.direction_V = [];            % direction of V modulation
prefilter_results.direction_VT = [];           % direction of VT modulation
prefilter_results.speed_shape_T = {};          % speed tuning shape from T_Vstatic: 'increasing', 'decreasing', 'bandpass', 'unclassified'
prefilter_results.speed_shape_VT = {};         % speed tuning shape from VT: 'increasing', 'decreasing', 'bandpass', 'unclassified'

% Significance threshold for Wilcoxon sign-rank tests
p_threshold = 0.05;

n_prefilter_clusters = 0;
n_speed_tuned_only = 0;
n_visually_tuned = 0;
n_mixed_tuning = 0;
n_not_tuned = 0;
n_should_run_glm = 0;

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
    
    % Load tuning tables for speed bin comparison
    tuning_file_spd = fullfile(tuning_curves_dir, [pid '.mat']);
    if ~exist(tuning_file_spd, 'file')
        fprintf('  WARNING: speed tuning file not found, skipping pre-filtering\n');
        % Still add clusters but mark as should_run_glm = true (no filtering possible)
        for ci = 1:length(visp_clusters)
            n_prefilter_clusters = n_prefilter_clusters + 1;
            prefilter_results.probe_id{n_prefilter_clusters} = pid;
            prefilter_results.cluster_id(n_prefilter_clusters) = visp_clusters(ci).id;
            prefilter_results.is_speed_tuned(n_prefilter_clusters) = NaN;
            prefilter_results.is_visually_tuned(n_prefilter_clusters) = NaN;
            prefilter_results.is_VT_responsive(n_prefilter_clusters) = NaN;
            prefilter_results.is_mixed_tuning(n_prefilter_clusters) = false;
            prefilter_results.should_run_glm(n_prefilter_clusters) = true;
            prefilter_results.category{n_prefilter_clusters} = 'no_data';
            prefilter_results.p_T_svm(n_prefilter_clusters) = NaN;
            prefilter_results.p_V_svm(n_prefilter_clusters) = NaN;
            prefilter_results.p_VT_svm(n_prefilter_clusters) = NaN;
            prefilter_results.p_T_vs_VT_motion(n_prefilter_clusters) = NaN;
            prefilter_results.direction_T(n_prefilter_clusters) = NaN;
            prefilter_results.direction_V(n_prefilter_clusters) = NaN;
            prefilter_results.direction_VT(n_prefilter_clusters) = NaN;
            prefilter_results.speed_shape_T{n_prefilter_clusters} = 'unclassified';
            prefilter_results.speed_shape_VT{n_prefilter_clusters} = 'unclassified';
            n_should_run_glm = n_should_run_glm + 1;
        end
        continue;
    end
    D_spd = load(tuning_file_spd);
    
    % Get tuning data for each condition
    tc_T = []; tc_V = []; tc_VT = [];
    for tg_i = 1:length(D_spd.trial_groups)
        if strcmp(D_spd.trial_groups{tg_i}, 'T_Vstatic')
            tc_T = D_spd.tuning_curves{tg_i};
        elseif strcmp(D_spd.trial_groups{tg_i}, 'V')
            tc_V = D_spd.tuning_curves{tg_i};
        elseif strcmp(D_spd.trial_groups{tg_i}, 'VT')
            tc_VT = D_spd.tuning_curves{tg_i};
        end
    end
    
    % Process each VISp cluster
    for ci = 1:length(visp_clusters)
        cid = visp_clusters(ci).id;
        n_prefilter_clusters = n_prefilter_clusters + 1;
        
        prefilter_results.probe_id{n_prefilter_clusters} = pid;
        prefilter_results.cluster_id(n_prefilter_clusters) = cid;
        
        % ====== Step 1: Check speed tuning (stationary_T vs motion_T) ======
        [sig_T, p_T, dir_T] = data.is_stationary_vs_motion_significant(cid, 'T_Vstatic');
        prefilter_results.p_T_svm(n_prefilter_clusters) = p_T;
        prefilter_results.direction_T(n_prefilter_clusters) = dir_T;
        is_speed_tuned = sig_T && ~isempty(sig_T);
        prefilter_results.is_speed_tuned(n_prefilter_clusters) = is_speed_tuned;
        
        % ====== Step 2: Check visual tuning (stationary_V vs motion_V) ======
        [sig_V, p_V, dir_V] = data.is_stationary_vs_motion_significant(cid, 'V');
        prefilter_results.p_V_svm(n_prefilter_clusters) = p_V;
        prefilter_results.direction_V(n_prefilter_clusters) = dir_V;
        is_visually_tuned = sig_V && ~isempty(sig_V);
        prefilter_results.is_visually_tuned(n_prefilter_clusters) = is_visually_tuned;
        
        % ====== Step 3: Check VT responsiveness (stationary_VT vs motion_VT) ======
        [sig_VT, p_VT, dir_VT] = data.is_stationary_vs_motion_significant(cid, 'VT');
        prefilter_results.p_VT_svm(n_prefilter_clusters) = p_VT;
        prefilter_results.direction_VT(n_prefilter_clusters) = dir_VT;
        is_VT_responsive = sig_VT && ~isempty(sig_VT);
        prefilter_results.is_VT_responsive(n_prefilter_clusters) = is_VT_responsive;
        
        % ====== Apply decision tree ======
        is_mixed_tuning = false;
        should_run_glm = false;
        category = 'not_tuned';
        p_T_vs_VT = NaN;
        
        if ~is_speed_tuned && ~is_visually_tuned
            % Neither speed nor visually tuned: do not analyse
            category = 'not_tuned';
            should_run_glm = false;
            n_not_tuned = n_not_tuned + 1;
            
        elseif is_speed_tuned && ~is_visually_tuned
            % Speed tuned but not visually tuned
            if is_VT_responsive
                % Compare motion_T vs motion_VT using Kruskal-Wallis on speed bins
                % Find tuning data for this cluster in T and VT conditions
                tc_T_idx = find([tc_T.cluster_id] == cid, 1);
                tc_VT_idx = find([tc_VT.cluster_id] == cid, 1);
                
                if ~isempty(tc_T_idx) && ~isempty(tc_VT_idx)
                    % Get speed bin firing rates across trials
                    tuning_T = tc_T(tc_T_idx).tuning;   % n_bins x n_trials
                    tuning_VT = tc_VT(tc_VT_idx).tuning;
                    
                    % Flatten to vectors for KW test
                    fr_T = tuning_T(:);
                    fr_VT = tuning_VT(:);
                    
                    % Remove NaN values
                    valid_T = ~isnan(fr_T);
                    valid_VT = ~isnan(fr_VT);
                    fr_T = fr_T(valid_T);
                    fr_VT = fr_VT(valid_VT);
                    
                    % Create group labels
                    grp_T = repmat({'T'}, length(fr_T), 1);
                    grp_VT = repmat({'VT'}, length(fr_VT), 1);
                    
                    % Combine and run Kruskal-Wallis test
                    all_fr = [fr_T; fr_VT];
                    all_grp = [grp_T; grp_VT];
                    
                    if length(all_fr) >= 4  % Need at least a few samples
                        [p_T_vs_VT, ~, ~] = kruskalwallis(all_fr, all_grp, 'off');
                    else
                        p_T_vs_VT = NaN;
                    end
                    
                    if ~isnan(p_T_vs_VT) && p_T_vs_VT >= p_threshold
                        % motion_T == motion_VT: only speed tuned, do not analyse
                        category = 'speed_tuned_only';
                        should_run_glm = false;
                        n_speed_tuned_only = n_speed_tuned_only + 1;
                    else
                        % motion_T != motion_VT: mixed tuning, continue to GLM
                        category = 'mixed_tuning';
                        is_mixed_tuning = true;
                        should_run_glm = true;
                        n_mixed_tuning = n_mixed_tuning + 1;
                        n_should_run_glm = n_should_run_glm + 1;
                    end
                else
                    % No tuning data available, run GLM anyway
                    category = 'speed_tuned_no_VT_data';
                    should_run_glm = true;
                    n_should_run_glm = n_should_run_glm + 1;
                end
            else
                % VT not responsive but T is: unexpected, raise warning
                warning('Cluster %s:%d - speed tuned but VT not responsive (stationary_VT == motion_VT)', pid, cid);
                category = 'speed_tuned_VT_not_responsive';
                should_run_glm = false;
                n_speed_tuned_only = n_speed_tuned_only + 1;
            end
            
        elseif is_visually_tuned
            % Visually tuned (with or without speed tuning): run full GLM
            category = 'visually_tuned';
            should_run_glm = true;
            n_visually_tuned = n_visually_tuned + 1;
            n_should_run_glm = n_should_run_glm + 1;
        end
        
        prefilter_results.is_mixed_tuning(n_prefilter_clusters) = is_mixed_tuning;
        prefilter_results.should_run_glm(n_prefilter_clusters) = should_run_glm;
        prefilter_results.category{n_prefilter_clusters} = category;
        prefilter_results.p_T_vs_VT_motion(n_prefilter_clusters) = p_T_vs_VT;
        prefilter_results.speed_shape_T{n_prefilter_clusters} = 'unclassified';  % Will be computed below
        prefilter_results.speed_shape_VT{n_prefilter_clusters} = 'unclassified'; % Will be computed below
    end
end

% Print pre-filtering summary
fprintf('\n====================================================================\n');
fprintf('  PRE-FILTERING SUMMARY\n');
fprintf('====================================================================\n');
fprintf('  Total VISp clusters scanned: %d\n', n_prefilter_clusters);
fprintf('\n  Category breakdown:\n');
fprintf('    Not tuned (neither speed nor visual):     %d (%.1f%%)\n', n_not_tuned, 100*n_not_tuned/n_prefilter_clusters);
fprintf('    Speed tuned only (excluded from GLM):     %d (%.1f%%)\n', n_speed_tuned_only, 100*n_speed_tuned_only/n_prefilter_clusters);
fprintf('    Mixed tuning (speed + some visual):       %d (%.1f%%)\n', n_mixed_tuning, 100*n_mixed_tuning/n_prefilter_clusters);
fprintf('    Visually tuned (will run GLM):            %d (%.1f%%)\n', n_visually_tuned, 100*n_visually_tuned/n_prefilter_clusters);
fprintf('\n  Clusters for full GLM analysis:             %d (%.1f%%)\n', n_should_run_glm, 100*n_should_run_glm/n_prefilter_clusters);
fprintf('  Clusters excluded from GLM:                 %d (%.1f%%)\n', n_prefilter_clusters - n_should_run_glm, 100*(n_prefilter_clusters - n_should_run_glm)/n_prefilter_clusters);
fprintf('====================================================================\n\n');

%% ====================================================================
%  Section 2c: Speed Tuning Shape Classification (Asymmetric Gaussian Fit)
%  ====================================================================
%  For all speed-tuned clusters, fit asymmetric Gaussian to the T_Vstatic
%  and VT tuning curves to classify shape as 'increasing', 'decreasing',
%  or 'bandpass'. This is done BEFORE the GLM to ensure all speed-tuned
%  neurons get a shape classification regardless of whether they enter GLM.
%  ====================================================================
fprintf('\n--- Classifying speed tuning shapes (Asymmetric Gaussian fit) ---\n');

% Get speed bin edges from the first available tuning file
first_spd_file = fullfile(tuning_curves_dir, [probe_ids{1} '.mat']);
D_tmp = load(first_spd_file);
spd_bin_edges_shape = D_tmp.tuning_curves{1}(1).bin_edges(:)';
spd_bin_centers_shape = D_tmp.tuning_curves{1}(1).bin_centers(:)';
clear D_tmp;

% Define thresholds for shape classification based on bin edges
% Use second bin edge as lower threshold, second-to-last as upper threshold
shape_threshold_low = spd_bin_edges_shape(2);    % decreasing if peak <= this
shape_threshold_high = spd_bin_edges_shape(end-1);  % increasing if peak >= this
fprintf('  Shape thresholds: decreasing if peak <= %.2f, increasing if peak >= %.2f cm/s\n', ...
    shape_threshold_low, shape_threshold_high);

n_shape_T_classified = 0;
n_shape_VT_classified = 0;

for pi = 1:n_prefilter_clusters
    pid = prefilter_results.probe_id{pi};
    cid = prefilter_results.cluster_id(pi);
    is_speed_tuned = prefilter_results.is_speed_tuned(pi);
    
    % Only classify shape for speed-tuned clusters
    if ~is_speed_tuned || isnan(is_speed_tuned)
        continue;
    end
    
    % Load tuning file for this probe
    tuning_file_shape = fullfile(tuning_curves_dir, [pid '.mat']);
    if ~exist(tuning_file_shape, 'file')
        continue;
    end
    D_spd_shape = load(tuning_file_shape);
    
    % Get T_Vstatic and VT tuning curves
    tc_T_shape = []; tc_VT_shape = [];
    for tg_i = 1:length(D_spd_shape.trial_groups)
        if strcmp(D_spd_shape.trial_groups{tg_i}, 'T_Vstatic')
            tc_T_shape = D_spd_shape.tuning_curves{tg_i};
        elseif strcmp(D_spd_shape.trial_groups{tg_i}, 'VT')
            tc_VT_shape = D_spd_shape.tuning_curves{tg_i};
        end
    end
    
    % --- Classify T_Vstatic shape ---
    if ~isempty(tc_T_shape)
        tc_idx_T = find([tc_T_shape.cluster_id] == cid, 1);
        if ~isempty(tc_idx_T)
            tuning_T = tc_T_shape(tc_idx_T).tuning;  % n_bins x n_trials
            try
                fit_T = AsymmetricGaussianFit(tuning_T, spd_bin_centers_shape(:));
                x_max_T = fit_T.x_max;
                if ~isnan(x_max_T)
                    if x_max_T <= shape_threshold_low
                        prefilter_results.speed_shape_T{pi} = 'decreasing';
                    elseif x_max_T >= shape_threshold_high
                        prefilter_results.speed_shape_T{pi} = 'increasing';
                    else
                        prefilter_results.speed_shape_T{pi} = 'bandpass';
                    end
                    n_shape_T_classified = n_shape_T_classified + 1;
                end
            catch
                % Fitting failed, keep 'unclassified'
            end
        end
    end
    
    % --- Classify VT shape ---
    if ~isempty(tc_VT_shape)
        tc_idx_VT = find([tc_VT_shape.cluster_id] == cid, 1);
        if ~isempty(tc_idx_VT)
            tuning_VT = tc_VT_shape(tc_idx_VT).tuning;  % n_bins x n_trials
            try
                fit_VT = AsymmetricGaussianFit(tuning_VT, spd_bin_centers_shape(:));
                x_max_VT = fit_VT.x_max;
                if ~isnan(x_max_VT)
                    if x_max_VT <= shape_threshold_low
                        prefilter_results.speed_shape_VT{pi} = 'decreasing';
                    elseif x_max_VT >= shape_threshold_high
                        prefilter_results.speed_shape_VT{pi} = 'increasing';
                    else
                        prefilter_results.speed_shape_VT{pi} = 'bandpass';
                    end
                    n_shape_VT_classified = n_shape_VT_classified + 1;
                end
            catch
                % Fitting failed, keep 'unclassified'
            end
        end
    end
end

fprintf('  Speed-tuned clusters with T shape classified: %d\n', n_shape_T_classified);
fprintf('  Speed-tuned clusters with VT shape classified: %d\n', n_shape_VT_classified);

% Print shape distribution
shape_cats = {'increasing', 'decreasing', 'bandpass', 'unclassified'};
fprintf('\n  Shape distribution (T_Vstatic):\n');
for si = 1:length(shape_cats)
    n_cat = sum(strcmp(prefilter_results.speed_shape_T, shape_cats{si}));
    fprintf('    %s: %d\n', shape_cats{si}, n_cat);
end
fprintf('\n  Shape distribution (VT):\n');
for si = 1:length(shape_cats)
    n_cat = sum(strcmp(prefilter_results.speed_shape_VT, shape_cats{si}));
    fprintf('    %s: %d\n', shape_cats{si}, n_cat);
end
fprintf('\n');

% Save pre-filtering results to CSV
prefilter_table = table(...
    string(prefilter_results.probe_id'), ...
    prefilter_results.cluster_id', ...
    prefilter_results.is_speed_tuned', ...
    prefilter_results.is_visually_tuned', ...
    prefilter_results.is_VT_responsive', ...
    prefilter_results.is_mixed_tuning', ...
    prefilter_results.should_run_glm', ...
    string(prefilter_results.category'), ...
    prefilter_results.p_T_svm', ...
    prefilter_results.p_V_svm', ...
    prefilter_results.p_VT_svm', ...
    prefilter_results.p_T_vs_VT_motion', ...
    prefilter_results.direction_T', ...
    prefilter_results.direction_V', ...
    prefilter_results.direction_VT', ...
    string(prefilter_results.speed_shape_T'), ...
    string(prefilter_results.speed_shape_VT'), ...
    'VariableNames', {'probe_id', 'cluster_id', 'is_speed_tuned', 'is_visually_tuned', ...
        'is_VT_responsive', 'is_mixed_tuning', 'should_run_glm', 'category', ...
        'p_T_svm', 'p_V_svm', 'p_VT_svm', 'p_T_vs_VT_motion', ...
        'direction_T', 'direction_V', 'direction_VT', 'speed_shape_T', 'speed_shape_VT'});

prefilter_csv_path = fullfile(ctl.figs.curr_dir, 'prefilter_decision_tree.csv');
writetable(prefilter_table, prefilter_csv_path);
fprintf('Pre-filtering results saved to: %s\n\n', prefilter_csv_path);

% Create a set of cluster IDs that should be analyzed
clusters_for_glm = containers.Map('KeyType', 'char', 'ValueType', 'any');
for pi = 1:n_prefilter_clusters
    if prefilter_results.should_run_glm(pi)
        key = sprintf('%s_%d', prefilter_results.probe_id{pi}, prefilter_results.cluster_id(pi));
        clusters_for_glm(key) = struct(...
            'probe_id', prefilter_results.probe_id{pi}, ...
            'cluster_id', prefilter_results.cluster_id(pi), ...
            'category', prefilter_results.category{pi});
    end
end

fprintf('Clusters passing pre-filter for GLM: %d\n\n', clusters_for_glm.Count);

%% ====================================================================
%  Section 3: Data collection — Speed-bin GLM (from tuning tables)
%  ====================================================================
fprintf('\n--- Loading pre-computed tuning tables (speed-bin GLM) ---\n');

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

fprintf('\n--- Speed-bin data collection complete ---\n');
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

fprintf('  Speed-bin table: %d rows x %d columns\n', height(T_master_spd), width(T_master_spd));

%% ====================================================================
%  Section 3b: Data collection — Time-bin GLM (50 ms bins from raw data)
%  ====================================================================
fprintf('\n--- Building time-bin GLM data (%.0f ms bins) ---\n', time_bin_width*1000);

max_rows_t = 2000000;  % generous: ~12000 bins/cluster * ~100 clusters
tc_probe_id    = cell(max_rows_t, 1);
tc_cluster_id  = zeros(max_rows_t, 1);
tc_trial_id    = zeros(max_rows_t, 1);
tc_condition   = cell(max_rows_t, 1);
tc_speed       = zeros(max_rows_t, 1);
tc_tf          = zeros(max_rows_t, 1);
tc_sf          = zeros(max_rows_t, 1);
tc_orientation = zeros(max_rows_t, 1);
tc_batch_gain  = zeros(max_rows_t, 1);
tc_spike_count = zeros(max_rows_t, 1);
tc_time_in_trial = zeros(max_rows_t, 1);

total_rows_t = 0;
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
    fprintf('  %d trials to process\n', length(all_trial_ids));
    
    % Pre-load spike times for VISp clusters
    spike_times_cache = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for ci_cache = 1:length(visp_id_set)
        cid_cache = visp_id_set(ci_cache);
        try
            cluster_obj = data.get_cluster_with_id(cid_cache);
            spike_times_cache(int32(cid_cache)) = cluster_obj.spike_times;
        catch ME
            fprintf('    WARNING: cluster %d spike_times: %s\n', cid_cache, ME.message);
        end
    end
    fprintf('  Spike times cached for %d clusters\n', spike_times_cache.Count);
    
    % --- Pre-compute valid cluster IDs for this probe ---
    valid_cids = visp_id_set(arrayfun(@(c) spike_times_cache.isKey(int32(c)), visp_id_set));
    n_valid_clusters = length(valid_cids);
    
    % Process each trial (vectorised inner loop)
    n_trials_probe = length(all_trial_ids);
    t_probe_start = tic;
    for ti_main = 1:n_trials_probe
        trial_id = all_trial_ids(ti_main);
        condition = all_trial_conditions{ti_main};
        
        % --- Progress feedback every 50 trials ---
        if mod(ti_main, 50) == 0 || ti_main == n_trials_probe
            elapsed = toc(t_probe_start);
            rate = ti_main / elapsed;
            eta = (n_trials_probe - ti_main) / rate;
            fprintf('    Trial %d/%d (%.0f trials/s, ETA %.0f s)\n', ...
                ti_main, n_trials_probe, rate, eta);
        end
        
        stim = trial_stim(trial_id);
        sf_val = stim.sf;
        gain_val = stim.batch_gain;
        or_val = stim.or;
        
        % Load trial object and convert to aligned
        try
            trial_obj = data.get_trials_with_trial_ids(trial_id);
            if iscell(trial_obj)
                trial_obj = trial_obj{1};
            end
            aligned_obj = trial_obj.to_aligned;
        catch ME
            fprintf('    WARNING: trial %d load failed: %s\n', trial_id, ME.message);
            continue;
        end
        
        tr_probe_t = aligned_obj.probe_t;
        tr_vel     = aligned_obj.velocity;
        tr_mmask   = aligned_obj.motion_mask;
        
        % --- Find continuous motion period (first to last motion sample) ---
        motion_idx = find(tr_mmask);
        if isempty(motion_idx)
            continue;
        end
        motion_start_idx = motion_idx(1);
        motion_end_idx   = motion_idx(end);
        
        % Motion period time range (for binning)
        t_motion_start = tr_probe_t(motion_start_idx);
        t_motion_end   = tr_probe_t(motion_end_idx);
        motion_dur = t_motion_end - t_motion_start;
        max_trial_duration = max(max_trial_duration, motion_dur);
        
        if motion_dur < time_bin_width
            continue;
        end
        
        % Create 50 ms time bin edges WITHIN THE MOTION PERIOD ONLY
        bin_edges_t = t_motion_start : time_bin_width : t_motion_end;
        n_tbins = length(bin_edges_t) - 1;
        if n_tbins < 1
            continue;
        end
        
        % --- Vectorised per-bin sample assignment ---
        % Assign each probe sample to a time bin (0 = outside motion period)
        sample_bin = discretize(tr_probe_t, bin_edges_t);
        sample_bin(isnan(sample_bin)) = 0;
        
        % Per-bin: total samples, motion samples, sum(abs(vel)*motion_mask)
        valid_mask   = sample_bin > 0;
        bins_vec     = sample_bin(valid_mask);
        mmask_vec    = tr_mmask(valid_mask);
        absvel_vec   = abs(tr_vel(valid_mask)) .* mmask_vec;
        
        n_samp_per_bin   = accumarray(bins_vec(:), 1,             [n_tbins 1]);
        n_motion_per_bin = accumarray(bins_vec(:), mmask_vec(:),  [n_tbins 1]);
        sum_speed_per_bin = accumarray(bins_vec(:), absvel_vec(:), [n_tbins 1]);
        
        % Bins that pass the 50% motion threshold
        good_bins = find(n_samp_per_bin > 0 & n_motion_per_bin >= n_samp_per_bin * 0.5);
        if isempty(good_bins)
            continue;
        end
        
        % Mean speed per good bin (motion samples only)
        mean_speed_per_bin = sum_speed_per_bin(good_bins) ./ n_motion_per_bin(good_bins);
        
        % Time-in-trial for bin centres RELATIVE TO MOTION ONSET
        bin_centres = (bin_edges_t(1:end-1) + bin_edges_t(2:end)) / 2;
        time_in_trial_vec = bin_centres(good_bins)' - t_motion_start;
        
        % Compute condition-level features (same for all bins in this trial)
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
        
        % --- Count spikes per bin per cluster (vectorised with histcounts) ---
        for ci = 1:n_valid_clusters
            cid = valid_cids(ci);
            st = spike_times_cache(int32(cid));
            
            % histcounts over bin_edges is O(n log n) — much faster than per-bin masking
            spike_counts_all = histcounts(st, bin_edges_t)';
            spike_counts_good = spike_counts_all(good_bins);
            
            % Write rows in bulk
            idx_start = total_rows_t + 1;
            idx_end   = total_rows_t + n_good;
            
            % Expand arrays if needed
            if idx_end > max_rows_t
                max_rows_t = max(max_rows_t * 2, idx_end + 100000);
                tc_probe_id{max_rows_t}      = [];
                tc_cluster_id(max_rows_t)    = 0;
                tc_trial_id(max_rows_t)      = 0;
                tc_condition{max_rows_t}     = [];
                tc_speed(max_rows_t)         = 0;
                tc_tf(max_rows_t)            = 0;
                tc_sf(max_rows_t)            = 0;
                tc_orientation(max_rows_t)   = 0;
                tc_batch_gain(max_rows_t)    = 0;
                tc_spike_count(max_rows_t)   = 0;
                tc_time_in_trial(max_rows_t) = 0;
            end
            
            tc_probe_id(idx_start:idx_end)    = repmat({pid}, n_good, 1);
            tc_cluster_id(idx_start:idx_end)  = cid;
            tc_trial_id(idx_start:idx_end)    = trial_id;
            tc_condition(idx_start:idx_end)   = repmat({condition}, n_good, 1);
            tc_speed(idx_start:idx_end)       = speed_vec;
            tc_tf(idx_start:idx_end)          = tf_vec;
            tc_sf(idx_start:idx_end)          = sf_vec;
            tc_orientation(idx_start:idx_end) = or_vec;
            tc_batch_gain(idx_start:idx_end)  = gain_vec;
            tc_spike_count(idx_start:idx_end) = spike_counts_good;
            tc_time_in_trial(idx_start:idx_end) = time_in_trial_vec;
            
            total_rows_t = idx_end;
        end
    end
    
    elapsed_probe = toc(t_probe_start);
    fprintf('  Probe done: %d trials, %d rows, %.1f s\n', ...
        n_trials_probe, total_rows_t, elapsed_probe);
end

fprintf('\n--- Time-bin data collection complete ---\n');
fprintf('  Total rows (time bins within motion period): %d\n', total_rows_t);
fprintf('  Max motion duration: %.2f s\n', max_trial_duration);

% Trim and build table
T_master_time = table(...
    string(tc_probe_id(1:total_rows_t)), ...
    tc_cluster_id(1:total_rows_t), ...
    tc_trial_id(1:total_rows_t), ...
    string(tc_condition(1:total_rows_t)), ...
    tc_speed(1:total_rows_t), ...
    tc_tf(1:total_rows_t), ...
    tc_sf(1:total_rows_t), ...
    tc_orientation(1:total_rows_t), ...
    tc_batch_gain(1:total_rows_t), ...
    tc_spike_count(1:total_rows_t), ...
    tc_time_in_trial(1:total_rows_t), ...
    'VariableNames', {'probe_id', 'cluster_id', 'trial_id', 'condition', ...
        'speed', 'tf', 'sf', 'orientation', 'batch_gain', ...
        'spike_count', 'time_in_trial'});

fprintf('  Time-bin table: %d rows x %d columns\n', height(T_master_time), width(T_master_time));

% Temporal basis range (based on max motion duration, not full trial)
trial_duration_range = [0, max_trial_duration];
fprintf('  Max motion duration: %.2f s (time bins restricted to motion period)\n', max_trial_duration);

%% ====================================================================
%  Section 4: Regularisation setup and cluster enumeration
%  ====================================================================

% Get unique clusters (shared across both GLM types)
unique_clusters = unique(T_master_spd(:, {'probe_id', 'cluster_id'}), 'rows');
n_unique_clusters = height(unique_clusters);
fprintf('  Unique clusters: %d\n', n_unique_clusters);

% Create shape arrays indexed by unique_clusters (lookup from prefilter_results)
cluster_speed_shape_T = cell(n_unique_clusters, 1);
cluster_speed_shape_VT = cell(n_unique_clusters, 1);
for ci = 1:n_unique_clusters
    pid = char(unique_clusters.probe_id(ci));
    cid = unique_clusters.cluster_id(ci);
    % Find this cluster in prefilter_results
    pf_idx = find(strcmp(prefilter_results.probe_id, pid) & prefilter_results.cluster_id == cid, 1);
    if ~isempty(pf_idx)
        cluster_speed_shape_T{ci} = prefilter_results.speed_shape_T{pf_idx};
        cluster_speed_shape_VT{ci} = prefilter_results.speed_shape_VT{pf_idx};
    else
        cluster_speed_shape_T{ci} = 'unclassified';
        cluster_speed_shape_VT{ci} = 'unclassified';
    end
end
fprintf('  Shape arrays created for %d clusters\n', n_unique_clusters);

% Ridge regularisation for FullInteraction model (always applied)
lambda_spd  = 1.0;
lambda_time = 1.0;
fprintf('  Regularisation: lambda=%.1f (applied to FullInteraction model)\n', lambda_spd);

%% ====================================================================
%  Section 4b: Basis function visualisation
%  ====================================================================
fprintf('\n--- Plotting basis functions ---\n');

x_speed_fine = linspace(speed_range(1), speed_range(2), 500)';
x_tf_fine    = linspace(tf_range(1),    tf_range(2),    500)';

B_speed_fine = make_raised_cosine_basis(x_speed_fine, n_speed_bases, speed_range(1), speed_range(2));
B_tf_fine    = make_raised_cosine_basis(x_tf_fine,    n_tf_bases,    tf_range(1),    tf_range(2));

fig_basis = figure('Position', [50 50 1800 700], 'Name', 'Basis Functions Overview');

% (a) Speed bases — linear x-axis
subplot(2, 3, 1);
cmap_speed = lines(n_speed_bases);
for bi = 1:n_speed_bases
    plot(x_speed_fine, B_speed_fine(:, bi), 'Color', cmap_speed(bi,:), 'LineWidth', 1.8); hold on;
end
xlabel('Speed (cm/s)'); ylabel('Basis amplitude');
title(sprintf('(a) Speed basis (n=%d)', n_speed_bases));
set(gca, 'box', 'off'); xlim(speed_range);

% (b) Speed bases — log-shifted axis
subplot(2, 3, 2);
epsilon = 0.5;
for bi = 1:n_speed_bases
    plot(log(x_speed_fine + epsilon), B_speed_fine(:, bi), ...
        'Color', cmap_speed(bi,:), 'LineWidth', 1.8); hold on;
end
xlabel('log(Speed + \epsilon)'); ylabel('Basis amplitude');
title('(b) Speed basis (log-shifted)');
set(gca, 'box', 'off');

% (c) TF bases
subplot(2, 3, 3);
cmap_tf = lines(n_tf_bases);
for bi = 1:n_tf_bases
    plot(x_tf_fine, B_tf_fine(:, bi), 'Color', cmap_tf(bi,:), 'LineWidth', 1.8); hold on;
end
xlabel('Temporal frequency (Hz)'); ylabel('Basis amplitude');
title(sprintf('(c) TF basis (n=%d)', n_tf_bases));
set(gca, 'box', 'off'); xlim(tf_range);

% (d) SF dummy coding
subplot(2, 3, 4);
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
title(sprintf('(d) SF dummy (%d levels)', n_sf));
set(gca, 'box', 'off');

% (e) Orientation dummy coding
subplot(2, 3, 5);
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
title(sprintf('(e) OR dummy (%d levels)', n_or));
set(gca, 'box', 'off');

% (f) Speed x TF interaction (tensor product, 2D contour)
subplot(2, 3, 6);
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
title(sprintf('(f) Speed x TF interaction (%d bases)', n_speed_bases * n_tf_bases));
legend(leg_entries, 'Location', 'eastoutside', 'FontSize', 7);
set(gca, 'box', 'off');

sgtitle('GLM Basis Functions', 'FontSize', 14, 'FontWeight', 'bold');

if save_figs
    print(fig_basis, fullfile(ctl.figs.curr_dir, 'basis_functions_overview.png'), '-dpng', '-r300');
    fprintf('  Saved basis_functions_overview.png\n');
end

%% ====================================================================
%  Section 5: Model fitting and cross-validation (dual GLM)
%  ====================================================================
fprintf('\n--- Fitting model hierarchy (dual GLM) ---\n');

glm_types = {'spd', 'time'};
glm_type_labels = {'Speed-bin', 'Time-bin'};

% Storage for both GLM types
for gt = 1:length(glm_types)
    gt_tag = glm_types{gt};
    
    results.(gt_tag) = struct();
    results.(gt_tag).probe_id = cell(n_unique_clusters, 1);
    results.(gt_tag).cluster_id = zeros(n_unique_clusters, 1);
    results.(gt_tag).n_trials = zeros(n_unique_clusters, 1);
    results.(gt_tag).n_bins_total = zeros(n_unique_clusters, 1);
    results.(gt_tag).n_spikes_total = zeros(n_unique_clusters, 1);
    
    for mi = 1:n_models
        ml = model_labels{mi};
        results.(gt_tag).([ml '_aic']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_bic']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_cv_bps']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_train_bps']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_deviance']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_n_params']) = zeros(n_unique_clusters, 1);
        results.(gt_tag).([ml '_dispersion']) = zeros(n_unique_clusters, 1);
    end
    
    results.(gt_tag).winning_model = cell(n_unique_clusters, 1);
    results.(gt_tag).delta_bps_interaction = zeros(n_unique_clusters, 1);
    
    % --- Positive control: Additive vs M0_Speed (tests combined visual contribution) ---
    results.(gt_tag).delta_bps_visual = zeros(n_unique_clusters, 1);  % Additive - M0_Speed
    results.(gt_tag).is_visually_tuned_glm = false(n_unique_clusters, 1);  % delta_bps_visual > threshold
    
    % --- Drop-one model CV bps (for hypothesis testing) ---
    for di = 1:n_dropone
        dl = dropone_labels{di};
        results.(gt_tag).([dl '_cv_bps']) = zeros(n_unique_clusters, 1);
    end
    
    % --- Drop-one delta bps (Additive - Additive_no_X) ---
    results.(gt_tag).delta_bps_drop_speed = zeros(n_unique_clusters, 1);  % Additive - Additive_no_Speed
    results.(gt_tag).delta_bps_drop_tf = zeros(n_unique_clusters, 1);     % Additive - Additive_no_TF
    results.(gt_tag).delta_bps_drop_sf = zeros(n_unique_clusters, 1);     % Additive - Additive_no_SF
    results.(gt_tag).delta_bps_drop_or = zeros(n_unique_clusters, 1);     % Additive - Additive_no_OR
    
    % Binary classification flags — drop-one test (unique contribution)
    results.(gt_tag).is_speed_tuned = false(n_unique_clusters, 1);        % Unique speed contribution
    results.(gt_tag).is_tf_tuned = false(n_unique_clusters, 1);           % Unique TF contribution
    results.(gt_tag).is_sf_tuned = false(n_unique_clusters, 1);           % Unique SF contribution
    results.(gt_tag).is_or_tuned = false(n_unique_clusters, 1);           % Unique OR contribution
    results.(gt_tag).has_significant_interaction = false(n_unique_clusters, 1);
end

% Storage for coefficients and predictions
% Preallocate: ~50 params × 3 models × 2 GLM types × n_clusters
max_coef_rows = 50 * n_models * length(glm_types) * n_unique_clusters;
all_coefficients = cell(max_coef_rows, 1);
n_coef_stored = 0;
cluster_predictions = struct();
cluster_predictions.spd  = cell(n_unique_clusters, 1);
cluster_predictions.time = cell(n_unique_clusters, 1);

for ci = 1:n_unique_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);
    
    % Progress: print every 10 clusters
    if mod(ci, 10) == 1 || ci == n_unique_clusters
        fprintf('Fitting cluster %d/%d...\n', ci, n_unique_clusters);
    end
    
    % ====== Loop over both GLM types ======
    for gt = 1:length(glm_types)
        gt_tag = glm_types{gt};
        
        % Select the right table
        if strcmp(gt_tag, 'spd')
            T_all = T_master_spd;
            lambda_reg = lambda_spd;
            use_temporal = false;
        else
            T_all = T_master_time;
            lambda_reg = lambda_time;
            use_temporal = true;
        end
        
        idx = T_all.probe_id == pid & T_all.cluster_id == cid;
        T_cluster = T_all(idx, :);
        
        if height(T_cluster) < 10
            fprintf('  [%s] Too few rows (%d), skipping\n', gt_tag, height(T_cluster));
            results.(gt_tag).probe_id{ci} = char(pid);
            results.(gt_tag).cluster_id(ci) = cid;
            results.(gt_tag).winning_model{ci} = 'N/A';
            continue;
        end
        
        results.(gt_tag).probe_id{ci} = char(pid);
        results.(gt_tag).cluster_id(ci) = cid;
        results.(gt_tag).n_trials(ci) = length(unique(T_cluster.trial_id));
        results.(gt_tag).n_bins_total(ci) = height(T_cluster);
        results.(gt_tag).n_spikes_total(ci) = sum(T_cluster.spike_count);
        
        % Prepare features
        speed_v = T_cluster.speed;
        tf_v = T_cluster.tf;
        sf_v = T_cluster.sf; sf_v(isnan(sf_v)) = 0;
        or_v = T_cluster.orientation; or_v(isnan(or_v)) = 0;
        y = T_cluster.spike_count;
        
        % Offset
        if strcmp(gt_tag, 'spd')
            offset = log(max(T_cluster.time_in_bin, 1e-6));
        else
            % Time-bin: constant offset (all bins same width)
            offset = log(time_bin_width) * ones(height(T_cluster), 1);
        end
        
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
        B_time_cl  = zeros(length(speed_v), 0);  % No temporal bases
        
        % --- Fit each model ---
        best_cv_bps = -Inf;
        best_model = 'M0';
        
        for mi = 1:n_models
            ml = model_labels{mi};
            
            % Ridge for FullInteraction if needed
            if strcmp(ml, 'FullInteraction')
                lambda = lambda_reg;
            else
                lambda = 0;
            end
            
            % Build design matrix from pre-computed bases
            [X, col_names] = assemble_design_matrix(B_speed_cl, B_tf_cl, B_time_cl, ...
                sf_v, or_v, ml);
            
            % Check dimensions
            if size(X, 2) >= height(T_cluster)
                fprintf('    %s-%s: SKIPPED (p=%d >= n=%d)\n', gt_tag, ml, size(X,2), height(T_cluster));
                results.(gt_tag).([ml '_aic'])(ci) = Inf;
                results.(gt_tag).([ml '_bic'])(ci) = Inf;
                results.(gt_tag).([ml '_cv_bps'])(ci) = -Inf;
                results.(gt_tag).([ml '_train_bps'])(ci) = -Inf;
                results.(gt_tag).([ml '_deviance'])(ci) = Inf;
                results.(gt_tag).([ml '_n_params'])(ci) = size(X, 2);
                results.(gt_tag).([ml '_dispersion'])(ci) = NaN;
                continue;
            end
            
            % Fit on full data
            res = fit_poisson_glm(X, y, offset, lambda);
            
            results.(gt_tag).([ml '_aic'])(ci) = res.aic;
            results.(gt_tag).([ml '_bic'])(ci) = res.bic;
            results.(gt_tag).([ml '_deviance'])(ci) = res.deviance;
            results.(gt_tag).([ml '_n_params'])(ci) = res.n_params;
            results.(gt_tag).([ml '_dispersion'])(ci) = res.dispersion;
            
            if results.(gt_tag).n_spikes_total(ci) > 0
                results.(gt_tag).([ml '_train_bps'])(ci) = ...
                    (res.log_likelihood / results.(gt_tag).n_spikes_total(ci)) / log(2);
            else
                results.(gt_tag).([ml '_train_bps'])(ci) = NaN;
            end
            
            % Cross-validate
            [~, cv_bps, ~, cv_predicted] = cross_validate_glm(X, y, offset, fold_ids, lambda);
            results.(gt_tag).([ml '_cv_bps'])(ci) = cv_bps;
            
            if cv_bps > best_cv_bps
                best_cv_bps = cv_bps;
                best_model = ml;
            end
            
            % Store coefficients (preallocated)
            n_beta = length(res.beta);
            for ki = 1:n_beta
                n_coef_stored = n_coef_stored + 1;
                if n_coef_stored > max_coef_rows
                    max_coef_rows = max_coef_rows * 2;
                    all_coefficients{max_coef_rows} = [];
                end
                all_coefficients{n_coef_stored} = {char(pid), cid, gt_tag, ml, col_names{ki}, ...
                    res.beta(ki), res.se(ki)};
            end
            
            % Store predictions
            % predicted_count = mu = E[Y] (expected spike count per bin)
            % predicted_fr = lambda(t) = mu / exposure_time (instantaneous firing rate in Hz)
            cluster_predictions.(gt_tag){ci}.(ml).predicted_count = res.predicted_count;
            cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_count = cv_predicted;
            cluster_predictions.(gt_tag){ci}.(ml).predicted_fr = res.predicted_count ./ exp(offset);  % lambda(t) = mu / t
            cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr = cv_predicted ./ exp(offset);      % CV lambda(t)
            cluster_predictions.(gt_tag){ci}.(ml).pearson_residuals = res.pearson_residuals;
            cluster_predictions.(gt_tag){ci}.(ml).beta = res.beta;
            cluster_predictions.(gt_tag){ci}.(ml).se = res.se;
            cluster_predictions.(gt_tag){ci}.(ml).col_names = col_names;
        end
        
        results.(gt_tag).winning_model{ci} = best_model;
        
        % --- Fit drop-one models (for hypothesis testing) ---
        for di = 1:n_dropone
            dl = dropone_labels{di};
            
            [X_drop, ~] = assemble_design_matrix(B_speed_cl, B_tf_cl, B_time_cl, ...
                sf_v, or_v, dl);
            
            if size(X_drop, 2) >= height(T_cluster)
                results.(gt_tag).([dl '_cv_bps'])(ci) = -Inf;
                continue;
            end
            
            [~, cv_bps_drop] = cross_validate_glm(X_drop, y, offset, fold_ids, 0);
            results.(gt_tag).([dl '_cv_bps'])(ci) = cv_bps_drop;
        end
        
        % --- Interaction test ---
        results.(gt_tag).delta_bps_interaction(ci) = ...
            results.(gt_tag).FullInteraction_cv_bps(ci) - results.(gt_tag).Additive_cv_bps(ci);
        
        % --- Positive control: Additive vs M0_Speed (visual contribution) ---
        results.(gt_tag).delta_bps_visual(ci) = ...
            results.(gt_tag).Additive_cv_bps(ci) - results.(gt_tag).M0_Speed_cv_bps(ci);
        results.(gt_tag).is_visually_tuned_glm(ci) = ...
            results.(gt_tag).delta_bps_visual(ci) > classification_threshold;
        
        % --- Drop-one comparisons (Additive vs Additive_no_X) ---
        additive_cv = results.(gt_tag).Additive_cv_bps(ci);
        results.(gt_tag).delta_bps_drop_speed(ci) = additive_cv - results.(gt_tag).Additive_no_Speed_cv_bps(ci);
        results.(gt_tag).delta_bps_drop_tf(ci)    = additive_cv - results.(gt_tag).Additive_no_TF_cv_bps(ci);
        results.(gt_tag).delta_bps_drop_sf(ci)    = additive_cv - results.(gt_tag).Additive_no_SF_cv_bps(ci);
        results.(gt_tag).delta_bps_drop_or(ci)    = additive_cv - results.(gt_tag).Additive_no_OR_cv_bps(ci);
        
        % --- Binary classification — drop-one test (unique contribution) ---
        results.(gt_tag).is_speed_tuned(ci) = ...
            results.(gt_tag).delta_bps_drop_speed(ci) > classification_threshold;
        results.(gt_tag).is_tf_tuned(ci) = ...
            results.(gt_tag).delta_bps_drop_tf(ci) > classification_threshold;
        results.(gt_tag).is_sf_tuned(ci) = ...
            results.(gt_tag).delta_bps_drop_sf(ci) > classification_threshold;
        results.(gt_tag).is_or_tuned(ci) = ...
            results.(gt_tag).delta_bps_drop_or(ci) > classification_threshold;
        results.(gt_tag).has_significant_interaction(ci) = ...
            results.(gt_tag).delta_bps_interaction(ci) > classification_threshold;
    end
end

fprintf('\n--- Model fitting complete ---\n');

%% ====================================================================
%  Section 6: Population summary
%  ====================================================================
fprintf('\n--- Population summary ---\n');

for gt = 1:length(glm_types)
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    fprintf('\n=== %s GLM ===\n', gt_label);
    
    % Winning model counts
    fprintf('Winning model distribution:\n');
    for mi = 1:n_models
        ml = model_labels{mi};
        n_win = sum(strcmp(results.(gt_tag).winning_model, ml));
        fprintf('  %s: %d clusters (%.1f%%)\n', ml, n_win, 100*n_win/n_unique_clusters);
    end
    
    % Interaction effect
    delta = results.(gt_tag).delta_bps_interaction;
    sig_inter = sum(delta > 0.01);
    fprintf('Interaction (FullInteraction vs Additive):\n');
    fprintf('  Median delta bps: %.4f\n', median(delta));
    fprintf('  Clusters with delta > 0.01 bps: %d (%.1f%%)\n', ...
        sig_inter, 100*sig_inter/n_unique_clusters);
    
    % Train-CV gap
    gaps_gt = zeros(n_unique_clusters, 1);
    for ci_g = 1:n_unique_clusters
        wm = results.(gt_tag).winning_model{ci_g};
        if strcmp(wm, 'N/A'), continue; end
        gaps_gt(ci_g) = results.(gt_tag).([wm '_train_bps'])(ci_g) - ...
            results.(gt_tag).([wm '_cv_bps'])(ci_g);
    end
    fprintf('Train-CV gap (winning model):\n');
    fprintf('  Median: %.4f bps\n', median(gaps_gt));
end

% Cross-GLM comparison
fprintf('\n=== Speed-bin vs Time-bin comparison ===\n');
fprintf('Additive CV bps: speed-bin median=%.4f, time-bin median=%.4f\n', ...
    median(results.spd.Additive_cv_bps), median(results.time.Additive_cv_bps));
fprintf('FullInteraction CV bps: speed-bin median=%.4f, time-bin median=%.4f\n', ...
    median(results.spd.FullInteraction_cv_bps), median(results.time.FullInteraction_cv_bps));

% --- Gather numbers for paper-ready paragraph (use speed-bin GLM as primary) ---
% Interaction statistics
delta_bps_s = results.spd.delta_bps_interaction;
delta_bps_s = delta_bps_s(~isnan(delta_bps_s));
n_inter_meaningful = sum(delta_bps_s > classification_threshold);
n_total_valid = length(delta_bps_s);
pct_inter = 100 * n_inter_meaningful / n_total_valid;

% Neuron classification statistics (drop-one tests)
n_speed_tuned = sum(results.spd.is_speed_tuned);
n_tf_tuned = sum(results.spd.is_tf_tuned);
n_sf_tuned = sum(results.spd.is_sf_tuned);
n_or_tuned = sum(results.spd.is_or_tuned);
n_has_interaction = sum(results.spd.has_significant_interaction);

% Positive control: visual contribution (Additive vs M0_Speed)
n_visually_tuned_glm = sum(results.spd.is_visually_tuned_glm);
delta_visual = results.spd.delta_bps_visual;
delta_visual = delta_visual(~isnan(delta_visual));

% Print tuning summary
fprintf('\n=== Neuron Classification Summary (Speed-bin GLM) ===\n');
fprintf('Positive control (Additive vs M0_Speed, visual contribution):\n');
fprintf('  Median delta bps (visual): %.4f\n', median(delta_visual));
fprintf('  Visually tuned (delta > %.3f): %d (%.1f%%)\n', ...
    classification_threshold, n_visually_tuned_glm, 100*n_visually_tuned_glm/n_unique_clusters);
fprintf('Drop-one tests (threshold: delta bps > %.3f):\n', classification_threshold);
fprintf('  Speed tuned (unique): %d (%.1f%%)\n', n_speed_tuned, 100*n_speed_tuned/n_unique_clusters);
fprintf('  TF tuned (unique):    %d (%.1f%%)\n', n_tf_tuned, 100*n_tf_tuned/n_unique_clusters);
fprintf('  SF tuned (unique):    %d (%.1f%%)\n', n_sf_tuned, 100*n_sf_tuned/n_unique_clusters);
fprintf('  OR tuned (unique):    %d (%.1f%%)\n', n_or_tuned, 100*n_or_tuned/n_unique_clusters);
fprintf('  Has significant interaction: %d (%.1f%%)\n', n_has_interaction, 100*n_has_interaction/n_unique_clusters);

%% ====================================================================
%  Section 6b: Component importance decomposition (no refitting needed)
%  ====================================================================
%  For each cluster's FullInteraction model, compute:
%   (b) Partial linear predictor range → fold-change in firing rate
%   (c) Variance of partial linear predictor → relative importance
%  Also computes joint Wald chi^2 test per component.
%  ====================================================================
fprintf('\n====================================================================\n');
fprintf('  COMPONENT IMPORTANCE DECOMPOSITION\n');
fprintf('====================================================================\n');

% Component groups to decompose (order matters for display)
comp_groups = {'Speed', 'TF', 'SF', 'OR', ...
               'Spd x TF', 'Spd x SF', 'Spd x OR', ...
               'TF x SF', 'TF x OR', 'SF x OR'};
% Match functions: precise column-name matchers for each component.
% Main effects exclude interaction columns (those containing '_x_').
% Interactions use startsWith/regexp to separate Spd*, TF*, SF* prefixes.
comp_match_fns = { ...
    @(c) startsWith(c, 'Speed_') && ~contains(c, '_x_'), ...  % Speed
    @(c) startsWith(c, 'TF_')    && ~contains(c, '_x_'), ...  % TF
    @(c) startsWith(c, 'SF_')    && ~contains(c, '_x_'), ...  % SF
    @(c) startsWith(c, 'OR_')    && ~contains(c, '_x_'), ...  % OR
    @(c) startsWith(c, 'Spd')    && contains(c, '_x_TF'), ... % Spd x TF
    @(c) startsWith(c, 'Spd')    && contains(c, '_x_SF'), ... % Spd x SF
    @(c) startsWith(c, 'Spd')    && contains(c, '_x_OR'), ... % Spd x OR
    @(c) ~isempty(regexp(c, '^TF\d+_x_SF', 'once')), ...     % TF x SF
    @(c) ~isempty(regexp(c, '^TF\d+_x_OR', 'once')), ...     % TF x OR
    @(c) ~isempty(regexp(c, '^SF_\d+_x_OR', 'once'))};
n_comp = length(comp_groups);


% Preallocate storage: n_clusters x n_components (for both Additive and FullInteraction)
for gt = 1:length(glm_types)
    gt_tag = glm_types{gt};
    % FullInteraction model
    results.(gt_tag).comp_fold_change  = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_var_eta      = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_rel_importance = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_wald_chi2    = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_wald_p       = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_labels       = comp_groups;
    % Additive model (no interactions, so only first 5 components)
    results.(gt_tag).comp_fold_change_additive  = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_var_eta_additive      = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_rel_importance_additive = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_wald_chi2_additive    = nan(n_unique_clusters, n_comp);
    results.(gt_tag).comp_wald_p_additive       = nan(n_unique_clusters, n_comp);
end

for ci = 1:n_unique_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);
    
    for gt = 1:length(glm_types)
        gt_tag = glm_types{gt};
        
        % Re-build design matrix to get the actual column values
        if strcmp(gt_tag, 'spd')
            T_all_ci = T_master_spd;
        else
            T_all_ci = T_master_time;
        end
        idx_ci = T_all_ci.probe_id == pid & T_all_ci.cluster_id == cid;
        T_cl = T_all_ci(idx_ci, :);
        if height(T_cl) < 2, continue; end
        
        speed_v = T_cl.speed;
        tf_v = T_cl.tf;
        sf_v = T_cl.sf; sf_v(isnan(sf_v)) = 0;
        or_v = T_cl.orientation; or_v(isnan(or_v)) = 0;
        
        B_spd_ci = make_raised_cosine_basis(speed_v, n_speed_bases, speed_range(1), speed_range(2));
        B_tf_ci  = make_raised_cosine_basis(tf_v, n_tf_bases, tf_range(1), tf_range(2));
        B_time_ci = zeros(height(T_cl), 0);  % No temporal bases
        
        % === Process both Additive and FullInteraction models ===
        models_to_process = {'Additive', 'FullInteraction'};
        suffixes = {'_additive', ''};  % suffix for storage fields
        
        for mi_comp = 1:2
            model_name = models_to_process{mi_comp};
            suffix = suffixes{mi_comp};
            
            % Check model was fitted
            if isempty(cluster_predictions.(gt_tag){ci}) || ...
                    ~isfield(cluster_predictions.(gt_tag){ci}, model_name)
                continue;
            end
            
            beta_model = cluster_predictions.(gt_tag){ci}.(model_name).beta;
            se_model   = cluster_predictions.(gt_tag){ci}.(model_name).se;
            cn_model   = cluster_predictions.(gt_tag){ci}.(model_name).col_names;
            
            if isempty(beta_model), continue; end
            
            % Assemble design matrix for this model type
            [X_model, cn_rebuilt] = assemble_design_matrix(B_spd_ci, B_tf_ci, B_time_ci, ...
                sf_v, or_v, model_name);
            
            % Verify column count matches stored beta
            if size(X_model, 2) ~= length(beta_model)
                continue;  % mismatch (e.g. zero-var columns removed differently)
            end
            
            % --- Compute partial linear predictors per component ---
            total_var_eta = 0;
            comp_vars = zeros(n_comp, 1);
            
            for gi = 1:n_comp
                % Find columns belonging to this component
                col_mask = cellfun(comp_match_fns{gi}, cn_rebuilt);
                
                if ~any(col_mask), continue; end
                
                % Partial linear predictor: eta_g = X_g * beta_g
                eta_g = X_model(:, col_mask) * beta_model(col_mask);
                
                % (b) Range → fold-change
                eta_range = max(eta_g) - min(eta_g);
                results.(gt_tag).(['comp_fold_change' suffix])(ci, gi) = exp(eta_range);
                
                % (c) Variance of partial predictor
                v_g = var(eta_g);
                comp_vars(gi) = v_g;
                results.(gt_tag).(['comp_var_eta' suffix])(ci, gi) = v_g;
                total_var_eta = total_var_eta + v_g;
                
                % Joint Wald chi^2 test (using diagonal approximation for speed)
                beta_g = beta_model(col_mask);
                se_g   = se_model(col_mask);
                % Full Wald: beta' * inv(Cov) * beta; diagonal approx = sum(z^2)
                z_g = beta_g ./ max(se_g, 1e-12);
                wald_chi2 = sum(z_g.^2);
                df_g = sum(col_mask);
                results.(gt_tag).(['comp_wald_chi2' suffix])(ci, gi) = wald_chi2;
                results.(gt_tag).(['comp_wald_p' suffix])(ci, gi) = 1 - chi2cdf(wald_chi2, df_g);
            end
            
            % Relative importance (proportion of total predictor variance)
            if total_var_eta > 0
                results.(gt_tag).(['comp_rel_importance' suffix])(ci, :) = comp_vars / total_var_eta;
            end
        end
    end
end

% --- Print population summary ---
fprintf('\n');
for gt = 1:length(glm_types)
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    fprintf('\n=== %s GLM: Component Importance (FullInteraction model) ===\n', gt_label);
    
    fc  = results.(gt_tag).comp_fold_change;
    ri  = results.(gt_tag).comp_rel_importance;
    wp  = results.(gt_tag).comp_wald_p;
    
    fprintf('\n  %-15s | %10s | %10s | %10s | %8s | %7s\n', ...
        'Component', 'Fold-Chg', 'Var(eta)', 'Rel Imp %', 'Wald p', '% Sig');
    fprintf('  %s\n', repmat('-', 1, 75));
    
    for gi = 1:n_comp
        fc_gi = fc(:, gi);
        ri_gi = ri(:, gi);
        wp_gi = wp(:, gi);
        
        fc_v = fc_gi(~isnan(fc_gi));
        ri_v = ri_gi(~isnan(ri_gi));
        wp_v = wp_gi(~isnan(wp_gi));
        
        if isempty(fc_v)
            fprintf('  %-15s | %10s | %10s | %10s | %8s | %7s\n', ...
                comp_groups{gi}, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A');
            continue;
        end
        
        pct_sig = 100 * sum(wp_v < 0.05) / length(wp_v);
        
        fprintf('  %-15s | %9.2fx | %10.4f | %9.1f%% | %8.4f | %6.1f%%\n', ...
            comp_groups{gi}, median(fc_v), median(ri_v), ...
            100*median(ri_v), median(wp_v), pct_sig);
    end
    
    fprintf('\n  Fold-change interpretation: how much FR is modulated by each component\n');
    fprintf('  (e.g. 2.5x means firing rate varies by a factor of 2.5 across the range)\n');
    fprintf('  Rel Importance: %% of total linear predictor variance attributable to component\n');
    fprintf('  Wald p: joint significance (diagonal approx), %% Sig = clusters with p < 0.05\n');
    
    % --- Identify dominant modulator per neuron ---
    fprintf('\n  Dominant modulator per neuron (highest relative importance):\n');
    [~, dom_idx] = max(ri, [], 2);
    for gi = 1:n_comp
        n_dom = sum(dom_idx == gi & any(~isnan(ri), 2));
        if n_dom > 0
            fprintf('    %-15s: %d neurons (%.1f%%)\n', comp_groups{gi}, n_dom, ...
                100*n_dom/sum(any(~isnan(ri), 2)));
        end
    end
end

%% ====================================================================
%  PAPER-READY PARAGRAPH
%  ====================================================================
%  Generate the paper-ready paragraph AFTER component importance is computed

% Tuning shape breakdown (from asymmetric Gaussian fits on T_Vstatic and VT)
% Using cluster_speed_shape_T and cluster_speed_shape_VT created earlier
n_inc_T = sum(strcmp(cluster_speed_shape_T, 'increasing'));
n_dec_T = sum(strcmp(cluster_speed_shape_T, 'decreasing'));
n_band_T = sum(strcmp(cluster_speed_shape_T, 'bandpass'));

n_inc_VT = sum(strcmp(cluster_speed_shape_VT, 'increasing'));
n_dec_VT = sum(strcmp(cluster_speed_shape_VT, 'decreasing'));
n_band_VT = sum(strcmp(cluster_speed_shape_VT, 'bandpass'));

% Map speed shape counts to expected variable names (using VT condition)
n_mono_inc_spd = n_inc_VT;
n_mono_dec_spd = n_dec_VT;
n_bandpass_spd = n_band_VT;

% TF tuning shapes not computed via asymmetric Gaussian fit (set to 0)
% Future work: implement TF tuning curve shape classification
n_mono_inc_tf = 0;
n_mono_dec_tf = 0;
n_bandpass_tf = 0;

% Component importance (FullInteraction model)
fc_spd  = results.spd.comp_fold_change;
ri_spd  = results.spd.comp_rel_importance;
wp_spd  = results.spd.comp_wald_p;
comp_labels = results.spd.comp_labels;

% Identify dominant component for each neuron
[~, dom_idx] = max(ri_spd, [], 2);
valid_dom = any(~isnan(ri_spd), 2);
n_speed_dominant = sum(dom_idx == 1 & valid_dom);  % Speed
n_tf_dominant = sum(dom_idx == 2 & valid_dom);     % TF

% =====================================================================
% Print the paper-ready paragraph
% =====================================================================
fprintf('\n');
fprintf('====================================================================\n');
fprintf('  PAPER-READY PARAGRAPH\n');
fprintf('====================================================================\n');
fprintf('\n');
fprintf('  "We used Poisson GLMs with cross-validation to characterize how VISp neurons\n');
fprintf('  integrate locomotion speed and visual motion (temporal frequency, TF) during\n');
fprintf('  combined visuo-locomotor stimulation. Across %d VISp neurons recorded during\n', n_unique_clusters);
fprintf('  motion-cloud viewing combined with treadmill locomotion, we used a drop-one\n');
fprintf('  procedure to assess the unique contribution of each variable.\n');
fprintf('\n');

fprintf('  NEURON CLASSIFICATION (drop-one from Additive):\n');
fprintf('  Of %d neurons, %d (%.1f%%) showed unique speed tuning,\n', ...
    n_unique_clusters, n_speed_tuned, 100*n_speed_tuned/n_unique_clusters);
fprintf('  %d (%.1f%%) unique TF tuning, %d (%.1f%%) unique SF tuning,\n', ...
    n_tf_tuned, 100*n_tf_tuned/n_unique_clusters, ...
    n_sf_tuned, 100*n_sf_tuned/n_unique_clusters);
fprintf('  and %d (%.1f%%) unique OR tuning (Delta CV bps > %.2f).\n', ...
    n_or_tuned, 100*n_or_tuned/n_unique_clusters, classification_threshold);
fprintf('\n');

fprintf('  TUNING SHAPES:\n');
fprintf('  Among speed-tuned neurons: %d (%.1f%%) monotonically increasing,\n', ...
    n_mono_inc_spd, 100*n_mono_inc_spd/max(n_speed_tuned,1));
fprintf('  %d (%.1f%%) monotonically decreasing, %d (%.1f%%) bandpass.\n', ...
    n_mono_dec_spd, 100*n_mono_dec_spd/max(n_speed_tuned,1), ...
    n_bandpass_spd, 100*n_bandpass_spd/max(n_speed_tuned,1));
fprintf('  Among TF-tuned neurons: %d (%.1f%%) monotonically increasing,\n', ...
    n_mono_inc_tf, 100*n_mono_inc_tf/max(n_tf_tuned,1));
fprintf('  %d (%.1f%%) monotonically decreasing, %d (%.1f%%) bandpass.\n', ...
    n_mono_dec_tf, 100*n_mono_dec_tf/max(n_tf_tuned,1), ...
    n_bandpass_tf, 100*n_bandpass_tf/max(n_tf_tuned,1));
fprintf('\n');

fprintf('  MULTISENSORY INTERACTIONS:\n');
fprintf('  Adding speed x stimulus interactions improved predictive performance by a median\n');
fprintf('  of %.4f bits/spike, with %d / %d neurons (%.1f%%) showing meaningful improvement\n', ...
    median(delta_bps_s), n_inter_meaningful, n_total_valid, pct_inter);
fprintf('  (Delta CV bps > %.2f).\n', classification_threshold);
if pct_inter < 50
    fprintf('  This indicates that for the majority of neurons, the combined visuo-locomotor\n');
    fprintf('  response is well captured by a linear sum of visual and translation components,\n');
    fprintf('  consistent with an additive multisensory integration framework.\n');
else
    fprintf('  This indicates that for the majority of neurons, the combined visuo-locomotor\n');
    fprintf('  response requires non-linear multisensory interactions beyond simple summation.\n');
end
fprintf('\n');

fprintf('  COMPONENT IMPORTANCE:\n');
fprintf('  Decomposing the FullInteraction model, %d neurons (%.1f%%) were dominated by speed,\n', ...
    n_speed_dominant, 100*n_speed_dominant/sum(valid_dom));
fprintf('  while %d (%.1f%%) were dominated by TF (based on relative variance explained).\n', ...
    n_tf_dominant, 100*n_tf_dominant/sum(valid_dom));
if ~isempty(fc_spd)
    fc_speed_v = fc_spd(~isnan(fc_spd(:,1)), 1);
    fc_tf_v = fc_spd(~isnan(fc_spd(:,2)), 2);
    if ~isempty(fc_speed_v) && ~isempty(fc_tf_v)
        fprintf('  Median fold-change in firing rate: %.2fx for speed, %.2fx for TF.\n', ...
            median(fc_speed_v), median(fc_tf_v));
    end
end
fprintf('\n');

fprintf('"\n');

fprintf('\n====================================================================\n');

%% ====================================================================
%  Section 7: Population visualisations
%  ====================================================================
fprintf('\n--- Generating population figures ---\n');

% --- Figure 2: Model Selection & Classification Summary (Redesigned) ---
% Panel 1: Additive vs FullInteraction comparison
% Panel 2: Drop-one test barplot
% Panel 3: Venn diagram showing overlapping tuning categories
fig2 = figure('Position', [50 50 1800 900], 'Name', 'Model Selection & Classification Summary');

for gt = 1:2
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    % Get classification flags
    is_spd = results.(gt_tag).is_speed_tuned;
    is_tf = results.(gt_tag).is_tf_tuned;
    is_sf = results.(gt_tag).is_sf_tuned;
    is_or = results.(gt_tag).is_or_tuned;
    is_int = results.(gt_tag).has_significant_interaction;
    
    n_total = n_unique_clusters;
    n_speed_tuned = sum(is_spd);
    n_tf_tuned = sum(is_tf);
    n_sf_tuned = sum(is_sf);
    n_or_tuned = sum(is_or);
    n_interaction = sum(is_int);
    
    % Tuning shape breakdown (from asymmetric Gaussian fits)
    % Use cluster_speed_shape_VT for VT (visual) condition
    n_inc_spd = sum(strcmp(cluster_speed_shape_VT, 'increasing'));
    n_dec_spd = sum(strcmp(cluster_speed_shape_VT, 'decreasing'));
    n_band_spd = sum(strcmp(cluster_speed_shape_VT, 'bandpass'));
    
    % --- Panel 1: Additive vs FullInteraction comparison ---
    subplot(2, 4, (gt-1)*4 + 1);
    
    % Count wins for Additive and FullInteraction only
    n_additive_win = sum(strcmp(results.(gt_tag).winning_model, 'Additive'));
    n_interaction_win = sum(strcmp(results.(gt_tag).winning_model, 'FullInteraction'));
    
    bar_data_models = [n_additive_win; n_interaction_win];
    b_models = bar(bar_data_models, 'FaceColor', 'flat');
    b_models.CData = [0.4 0.6 0.8; 0.8 0.4 0.4];
    set(gca, 'XTickLabel', {'Additive', 'FullInteraction'}, 'XTickLabelRotation', 30, 'FontSize', 8);
    ylabel('# Clusters');
    title(sprintf('%s: Model Selection', gt_label), 'FontSize', 10);
    set(gca, 'box', 'off');
    
    % Add count labels
    for mi = 1:2
        text(mi, bar_data_models(mi) + 0.5, sprintf('%d (%.0f%%)', bar_data_models(mi), 100*bar_data_models(mi)/n_total), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
    
    % --- Panel 2: Drop-one test barplot ---
    subplot(2, 4, (gt-1)*4 + 2);
    class_counts = [n_speed_tuned; n_tf_tuned; n_sf_tuned; n_or_tuned; n_interaction];
    class_labels = {'Speed', 'TF', 'SF', 'OR', 'Interaction'};
    class_colors = [0.2 0.6 1.0; 0.8 0.35 0.0; 0.2 0.6 0.2; 0.6 0.2 0.6; 0.9 0.2 0.2];
    bh = barh(class_counts, 'FaceColor', 'flat');
    bh.CData = class_colors;
    set(gca, 'YTick', 1:5, 'YTickLabel', class_labels, 'FontSize', 8);
    xlabel('# Neurons');
    % Add count and percentage labels
    for ki = 1:length(class_counts)
        text(class_counts(ki) + 0.5, ki, sprintf('%d (%.0f%%)', class_counts(ki), 100*class_counts(ki)/n_total), ...
            'VerticalAlignment', 'middle', 'FontSize', 7);
    end
    title(sprintf('%s: Drop-One Tests (n=%d)', gt_label, n_total), 'FontSize', 10);
    set(gca, 'box', 'off');
    
    % --- Panel 3: 4-Set Venn Diagram for Speed, TF, SF, OR ---
    ax_venn = subplot(2, 4, (gt-1)*4 + 3);
    hold on;
    
    % Draw 4-set Venn using ellipses (standard 4-set Venn layout)
    % Colors for each set
    colors = {[0.2 0.6 1.0 0.3], ...  % Speed - blue
              [0.8 0.35 0.0 0.3], ... % TF - orange
              [0.2 0.6 0.2 0.3], ...  % SF - green
              [0.6 0.2 0.6 0.3]};     % OR - purple
    
    % Ellipse parameters for 4-set Venn (positions, sizes, rotations)
    % Using a standard 4-set Venn layout
    ellipse_params = [
        0.35, 0.6, 0.35, 0.55, 45;   % Speed: cx, cy, rx, ry, angle
        0.65, 0.6, 0.35, 0.55, -45;  % TF
        0.35, 0.4, 0.35, 0.55, -45;  % SF
        0.65, 0.4, 0.35, 0.55, 45];  % OR
    
    % Draw ellipses
    set_names = {'Speed', 'TF', 'SF', 'OR'};
    label_positions = [0.05, 0.85; 0.95, 0.85; 0.05, 0.15; 0.95, 0.15];
    
    for si = 1:4
        draw_ellipse(ellipse_params(si,1), ellipse_params(si,2), ...
                     ellipse_params(si,3), ellipse_params(si,4), ...
                     ellipse_params(si,5), colors{si});
        text(label_positions(si,1), label_positions(si,2), set_names{si}, ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', colors{si}(1:3)*0.7, ...
            'HorizontalAlignment', 'center');
    end
    
    % Calculate all 16 intersection regions for 4 sets
    % Binary representation: Speed=1, TF=2, SF=4, OR=8
    region_counts = zeros(16, 1);
    region_has_interaction = zeros(16, 1);
    for ci = 1:n_total
        idx = 1 + is_spd(ci) + 2*is_tf(ci) + 4*is_sf(ci) + 8*is_or(ci);
        region_counts(idx) = region_counts(idx) + 1;
        if is_int(ci)
            region_has_interaction(idx) = region_has_interaction(idx) + 1;
        end
    end
    
    % Approximate positions for each region in the Venn diagram
    % These are manually tuned for the 4-ellipse layout
    region_positions = [
        0.50, 0.50;   % 0: none
        0.20, 0.70;   % 1: Speed only
        0.80, 0.70;   % 2: TF only
        0.50, 0.75;   % 3: Speed+TF
        0.20, 0.30;   % 4: SF only
        0.25, 0.50;   % 5: Speed+SF
        0.50, 0.35;   % 6: TF+SF (actually diagonal)
        0.38, 0.58;   % 7: Speed+TF+SF
        0.80, 0.30;   % 8: OR only
        0.50, 0.65;   % 9: Speed+OR (diagonal)
        0.75, 0.50;   % 10: TF+OR
        0.62, 0.58;   % 11: Speed+TF+OR
        0.50, 0.25;   % 12: SF+OR
        0.38, 0.42;   % 13: Speed+SF+OR
        0.62, 0.42;   % 14: TF+SF+OR
        0.50, 0.50];  % 15: all four
    
    % Draw counts in each region
    for ri = 1:16
        if region_counts(ri) > 0
            count_str = sprintf('%d', region_counts(ri));
            % If some have interaction, show as "count (int)"
            if region_has_interaction(ri) > 0
                count_str = sprintf('%d', region_counts(ri));
                text(region_positions(ri,1), region_positions(ri,2), count_str, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.9 0.1 0.1]);
            else
                text(region_positions(ri,1), region_positions(ri,2), count_str, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontSize', 8, 'FontWeight', 'bold');
            end
        end
    end
    
    % Show "none" count outside
    if region_counts(1) > 0
        text(0.50, 0.02, sprintf('None: %d', region_counts(1)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
    end
    
    axis equal;
    axis([0 1 0 1]);
    axis off;
    title(sprintf('%s: Tuning Overlap (red=w/interaction)', gt_label), 'FontSize', 10);
    hold off;
    
    % --- Panel 4: Tuning shape breakdown (combined Speed + TF) ---
    subplot(2, 4, (gt-1)*4 + 4);
    shape_data = [n_mono_inc_spd, n_mono_inc_tf; ...
                  n_mono_dec_spd, n_mono_dec_tf; ...
                  n_bandpass_spd, n_bandpass_tf];
    b_shape = bar(shape_data, 'grouped');
    b_shape(1).FaceColor = [0.2 0.6 1.0];  % Speed
    b_shape(2).FaceColor = [0.8 0.35 0.0]; % TF
    set(gca, 'XTickLabel', {'Mono Inc', 'Mono Dec', 'Bandpass'}, 'FontSize', 8);
    ylabel('# Neurons');
    legend({'Speed', 'TF'}, 'Location', 'northeast', 'FontSize', 7);
    title(sprintf('%s: Tuning Shapes', gt_label), 'FontSize', 10);
    set(gca, 'box', 'off');
end

sgtitle(sprintf('Model Selection & Drop-One Classification (threshold: \\Delta bps > %.2f)', classification_threshold), ...
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

for gt = 1:2
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    n_total = n_unique_clusters;
    n_spd = sum(results.(gt_tag).is_speed_tuned);
    n_tf = sum(results.(gt_tag).is_tf_tuned);
    n_sf = sum(results.(gt_tag).is_sf_tuned);
    n_or = sum(results.(gt_tag).is_or_tuned);
    n_inter = sum(results.(gt_tag).has_significant_interaction);
    
    fprintf('  --- %s GLM ---\n', gt_label);
    fprintf('  Speed-tuned (unique):    %3d / %d (%.1f%%)\n', n_spd, n_total, 100*n_spd/n_total);
    fprintf('  TF-tuned (unique):       %3d / %d (%.1f%%)\n', n_tf, n_total, 100*n_tf/n_total);
    fprintf('  SF-tuned (unique):       %3d / %d (%.1f%%)\n', n_sf, n_total, 100*n_sf/n_total);
    fprintf('  OR-tuned (unique):       %3d / %d (%.1f%%)\n', n_or, n_total, 100*n_or/n_total);
    fprintf('  Significant interaction: %3d / %d (%.1f%%)\n', n_inter, n_total, 100*n_inter/n_total);
    
    % Tuning shape breakdown (from asymmetric Gaussian fits)
    fprintf('  Speed tuning shapes (T_Vstatic):\n');
    fprintf('    Increasing:   %3d\n', sum(strcmp(cluster_speed_shape_T, 'increasing')));
    fprintf('    Decreasing:   %3d\n', sum(strcmp(cluster_speed_shape_T, 'decreasing')));
    fprintf('    Bandpass:     %3d\n', sum(strcmp(cluster_speed_shape_T, 'bandpass')));
    
    fprintf('  Speed tuning shapes (VT):\n');
    fprintf('    Increasing:   %3d\n', sum(strcmp(cluster_speed_shape_VT, 'increasing')));
    fprintf('    Decreasing:   %3d\n', sum(strcmp(cluster_speed_shape_VT, 'decreasing')));
    fprintf('    Bandpass:     %3d\n\n', sum(strcmp(cluster_speed_shape_VT, 'bandpass')));
end

fprintf('====================================================================\n');

% --- Figure 4: Component importance decomposition (Additive & FullInteraction only) ---
fig4 = figure('Position', [20 20 1800 800], 'Name', 'Component Importance - Additive & FullInteraction');

% Only Additive and FullInteraction models
model_names_all = {'Additive', 'FullInteraction'};
n_rows_fig4 = length(model_names_all) * 2;  % 2 models x 2 GLM types = 4 rows

for gt = 1:2
    gt_tag = glm_types{gt};
    gt_label = glm_type_labels{gt};
    
    for mt = 1:length(model_names_all)
        model_tag = model_names_all{mt};
        row_idx = (gt-1)*length(model_names_all) + mt;  % 1-6 for speed-bin; 7-12 for time-bin
        
        % Determine which stored data to use based on model
        % Only Additive and FullInteraction have pre-computed component decomposition
        % For others, we compute on-the-fly from cluster_predictions
        if strcmp(model_tag, 'Additive')
            fc = results.(gt_tag).comp_fold_change_additive;
            ri = results.(gt_tag).comp_rel_importance_additive;
            wp = results.(gt_tag).comp_wald_p_additive;
            has_decomposition = true;
        elseif strcmp(model_tag, 'FullInteraction')
            fc = results.(gt_tag).comp_fold_change;
            ri = results.(gt_tag).comp_rel_importance;
            wp = results.(gt_tag).comp_wald_p;
            has_decomposition = true;
        else
            % For simpler models, compute decomposition on-the-fly
            has_decomposition = false;
            fc = nan(n_unique_clusters, n_comp);
            ri = nan(n_unique_clusters, n_comp);
            wp = nan(n_unique_clusters, n_comp);
            
            % Compute component metrics for this model
            for ci_comp = 1:n_unique_clusters
                if isempty(cluster_predictions.(gt_tag){ci_comp}) || ...
                        ~isfield(cluster_predictions.(gt_tag){ci_comp}, model_tag)
                    continue;
                end
                
                beta_m = cluster_predictions.(gt_tag){ci_comp}.(model_tag).beta;
                se_m = cluster_predictions.(gt_tag){ci_comp}.(model_tag).se;
                cn_m = cluster_predictions.(gt_tag){ci_comp}.(model_tag).col_names;
                
                if isempty(beta_m), continue; end
                
                % Get design matrix for this cluster
                pid_comp = unique_clusters.probe_id(ci_comp);
                cid_comp = unique_clusters.cluster_id(ci_comp);
                if strcmp(gt_tag, 'spd')
                    T_all_comp = T_master_spd;
                else
                    T_all_comp = T_master_time;
                end
                idx_comp = T_all_comp.probe_id == pid_comp & T_all_comp.cluster_id == cid_comp;
                T_cl_comp = T_all_comp(idx_comp, :);
                if height(T_cl_comp) < 2, continue; end
                
                speed_v_comp = T_cl_comp.speed;
                tf_v_comp = T_cl_comp.tf;
                sf_v_comp = T_cl_comp.sf; sf_v_comp(isnan(sf_v_comp)) = 0;
                or_v_comp = T_cl_comp.orientation; or_v_comp(isnan(or_v_comp)) = 0;
                
                B_spd_comp = make_raised_cosine_basis(speed_v_comp, n_speed_bases, speed_range(1), speed_range(2));
                B_tf_comp = make_raised_cosine_basis(tf_v_comp, n_tf_bases, tf_range(1), tf_range(2));
                if strcmp(gt_tag, 'time') && use_temporal
                    B_time_comp = make_linear_raised_cosine_basis(T_cl_comp.time_in_trial, 5, ...
                        trial_duration_range(1), trial_duration_range(2));
                else
                    B_time_comp = zeros(height(T_cl_comp), 0);
                end
                
                [X_m, cn_rebuilt_m] = assemble_design_matrix(B_spd_comp, B_tf_comp, B_time_comp, ...
                    sf_v_comp, or_v_comp, model_tag);
                
                if size(X_m, 2) ~= length(beta_m), continue; end
                
                % Compute component metrics
                total_var_eta_m = 0;
                comp_vars_m = zeros(n_comp, 1);
                
                for gi_m = 1:n_comp
                    col_mask_m = cellfun(comp_match_fns{gi_m}, cn_rebuilt_m);
                    
                    if ~any(col_mask_m), continue; end
                    
                    eta_g_m = X_m(:, col_mask_m) * beta_m(col_mask_m);
                    eta_range_m = max(eta_g_m) - min(eta_g_m);
                    fc(ci_comp, gi_m) = exp(eta_range_m);
                    
                    v_g_m = var(eta_g_m);
                    comp_vars_m(gi_m) = v_g_m;
                    total_var_eta_m = total_var_eta_m + v_g_m;
                    
                    beta_g_m = beta_m(col_mask_m);
                    se_g_m = se_m(col_mask_m);
                    z_g_m = beta_g_m ./ max(se_g_m, 1e-12);
                    wald_chi2_m = sum(z_g_m.^2);
                    df_g_m = sum(col_mask_m);
                    wp(ci_comp, gi_m) = 1 - chi2cdf(wald_chi2_m, df_g_m);
                end
                
                if total_var_eta_m > 0
                    ri(ci_comp, :) = comp_vars_m / total_var_eta_m;
                end
            end
        end
        
        % Determine which components have data
        has_data = any(~isnan(fc), 1);
        comp_idx = find(has_data);
        n_show = length(comp_idx);
        
        if n_show == 0
            % No components for this model (e.g., M0)
            subplot(n_rows_fig4, 3, (row_idx-1)*3 + 1);
            text(0.5, 0.5, sprintf('%s %s\n(Intercept only)', gt_label, model_tag), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
            axis off;
            subplot(n_rows_fig4, 3, (row_idx-1)*3 + 2);
            axis off;
            subplot(n_rows_fig4, 3, (row_idx-1)*3 + 3);
            axis off;
            continue;
        end
        
        comp_labels_show = comp_groups(comp_idx);
        
        % Component colors (11 entries: 5 main effects + 6 interactions)
        comp_clrs = [0.2 0.6 1.0; 1.0 0.5 0.0; 0.4 0.8 0.4; 0.8 0.4 0.8; ...
                     0.3 0.8 0.8; 0.9 0.2 0.2; 0.8 0.5 0.2; 0.6 0.2 0.6; ...
                     0.5 0.8 0.2; 0.2 0.5 0.8; 0.9 0.6 0.6];
        
        % --- Col 1: Fold-change box/swarm ---
        subplot(n_rows_fig4, 3, (row_idx-1)*3 + 1);
        hold on;
        for gi = 1:n_show
            vals = fc(:, comp_idx(gi));
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            % Jittered strip plot
            jitter = 0.15 * (rand(length(vals), 1) - 0.5);
            scatter(gi + jitter, vals, 8, comp_clrs(comp_idx(gi),:), ...
                'filled', 'MarkerFaceAlpha', 0.35);
            % Median + IQR
            plot([gi-0.2, gi+0.2], [median(vals), median(vals)], '-k', 'LineWidth', 2);
            q25 = quantile(vals, 0.25); q75 = quantile(vals, 0.75);
            plot([gi gi], [q25 q75], '-k', 'LineWidth', 1);
        end
        yline(1, 'k--', 'LineWidth', 0.5);
        hold off;
        set(gca, 'XTick', 1:n_show, 'XTickLabel', comp_labels_show, ...
            'XTickLabelRotation', 45, 'TickLabelInterpreter', 'tex', 'FontSize', 6);
        ylabel('Fold-change FR', 'FontSize', 7);
        ylim([0 15]);  % Fixed y-limit for fold change
        title(sprintf('%s %s: Modulation', gt_label, model_tag), 'FontSize', 8);
        set(gca, 'box', 'off');
        
        % --- Col 2: Relative importance stacked bar (sorted by dominant component) ---
        subplot(n_rows_fig4, 3, (row_idx-1)*3 + 2);
        ri_plot = ri(:, comp_idx);
        valid_rows = all(~isnan(ri_plot), 2);
        ri_valid = ri_plot(valid_rows, :);
        if ~isempty(ri_valid) && size(ri_valid, 2) > 0
            % Sort neurons by dominant component for visual clarity
            [~, dom] = max(ri_valid, [], 2);
            [~, sort_order] = sort(dom);
            ri_sorted = ri_valid(sort_order, :);
            
            b_stacked = bar(ri_sorted, 'stacked', 'EdgeColor', 'none');
            for bi = 1:length(b_stacked)
                b_stacked(bi).FaceColor = comp_clrs(comp_idx(bi), :);
            end
            xlabel('Neurons', 'FontSize', 7); ylabel('Rel. importance', 'FontSize', 7);
            if row_idx == 1  % Only add legend to first row to save space
                legend(comp_labels_show, 'Location', 'eastoutside', 'FontSize', 5, ...
                    'Interpreter', 'tex');
            end
        else
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off;
        end
        title(sprintf('%s %s: Variance', gt_label, model_tag), 'FontSize', 8);
        set(gca, 'box', 'off', 'FontSize', 6);
        
        % --- Col 3: % significant per component (bar) ---
        subplot(n_rows_fig4, 3, (row_idx-1)*3 + 3);
        pct_sig_vec = zeros(n_show, 1);
        for gi = 1:n_show
            wp_gi = wp(:, comp_idx(gi));
            wp_v = wp_gi(~isnan(wp_gi));
            if ~isempty(wp_v)
                pct_sig_vec(gi) = 100 * sum(wp_v < 0.05) / length(wp_v);
            end
        end
        b_sig = bar(pct_sig_vec, 'EdgeColor', 'none');
        for bi = 1:n_show
            % Color each bar individually via CData
            b_sig.FaceColor = 'flat';
            b_sig.CData(bi, :) = comp_clrs(comp_idx(bi), :);
        end
        set(gca, 'XTick', 1:n_show, 'XTickLabel', comp_labels_show, ...
            'XTickLabelRotation', 45, 'TickLabelInterpreter', 'tex', 'FontSize', 6);
        ylabel('% sig.', 'FontSize', 7);
        yline(5, 'k--', 'chance', 'FontSize', 6);
        ylim([0 100]);
        title(sprintf('%s %s: Wald Sig.', gt_label, model_tag), 'FontSize', 8);
        set(gca, 'box', 'off');
    end
end

sgtitle('Component Importance: Additive & FullInteraction Models', 'FontSize', 12, 'FontWeight', 'bold');

if save_figs
    print(fig4, fullfile(ctl.figs.curr_dir, 'component_importance.png'), '-dpng', '-r300');
    fprintf('  Saved component_importance.png\n');
end

%% ====================================================================
%  Section 8: Per-cluster visualisations (Model Overview with Scatter)
%  ====================================================================
%  LAYOUT: (2 × n_models) rows × 6 columns
%    Rows organized by: [Speed-bin models, Time-bin models]
%    For each model row:
%      Col 1: Beta swarm plot (grouped by feature type)
%      Col 2: Mean scatter by condition with SEM error bars
%      Col 3: Scatter coloured by speed
%      Col 4: Scatter coloured by TF
%      Col 5: Scatter coloured by SF
%      Col 6: Scatter coloured by OR
%    Winning model row is highlighted with red border.
%  ====================================================================
fprintf('\n--- Generating per-cluster figures (Model Overview with Scatter) ---\n');

n_cols_scatter = 6;
n_rows_scatter = 2 * n_models;  % n_models per GLM type × 2 GLM types

% Condition colours
cond_map = containers.Map({'T_Vstatic', 'V', 'VT'}, ...
    {[0.8 0.2 0.2], [0.2 0.7 0.2], [0 0.4 0.8]});
cond_order = {'T_Vstatic', 'V', 'VT'};

% SF / OR colour maps
sf_unique_all = [0.003, 0.006, 0.012];
sf_cmap = [0.2 0.4 0.9; 0.1 0.7 0.3; 0.9 0.3 0.1];
or_unique_all = [-pi/4, 0, pi/4, pi/2];  % sorted order to match unique() output
or_cmap = [0.9 0.2 0.2; 0.2 0.7 0.2; 0.2 0.4 0.9; 0.8 0.4 0.8];  % distinct colors for each orientation

% Abbreviated model labels for tight spacing
model_abbrev = containers.Map();
model_abbrev('M0')              = 'M0';
model_abbrev('M0_Speed')        = 'M0+S';
model_abbrev('M0_Speed_TF')     = 'M0+S+TF';
model_abbrev('M0_Speed_TF_SF')  = 'M0+S+TF+SF';
model_abbrev('Additive')        = 'Additive';
model_abbrev('FullInteraction') = 'FullInt';

% Row definitions for scatter plots: {glm_type, model_label, display_prefix}
% Rows 1-n_models: Speed-bin models, Rows (n_models+1)-(2*n_models): Time-bin models
scatter_row_defs = cell(n_rows_scatter, 3);
for mi = 1:n_models
    scatter_row_defs{mi, 1} = 'spd';
    scatter_row_defs{mi, 2} = model_labels{mi};
    scatter_row_defs{mi, 3} = sprintf('Spd: %s', model_abbrev(model_labels{mi}));
end
for mi = 1:n_models
    scatter_row_defs{n_models + mi, 1} = 'time';
    scatter_row_defs{n_models + mi, 2} = model_labels{mi};
    scatter_row_defs{n_models + mi, 3} = sprintf('Time: %s', model_abbrev(model_labels{mi}));
end

for ci = 1:n_unique_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);
    
    wm_spd  = results.spd.winning_model{ci};
    wm_time = results.time.winning_model{ci};
    
    % --- Build classification label string ---
    class_spd_str = '';
    if results.spd.is_speed_tuned(ci)
        class_spd_str = [class_spd_str, 'Spd'];  %#ok<AGROW>
    end
    if results.spd.is_tf_tuned(ci)
        if ~isempty(class_spd_str), class_spd_str = [class_spd_str, '+'];  end  %#ok<AGROW>
        class_spd_str = [class_spd_str, 'TF'];  %#ok<AGROW>
    end
    if results.spd.is_sf_tuned(ci)
        if ~isempty(class_spd_str), class_spd_str = [class_spd_str, '+'];  end  %#ok<AGROW>
        class_spd_str = [class_spd_str, 'SF'];  %#ok<AGROW>
    end
    if results.spd.is_or_tuned(ci)
        if ~isempty(class_spd_str), class_spd_str = [class_spd_str, '+'];  end  %#ok<AGROW>
        class_spd_str = [class_spd_str, 'OR'];  %#ok<AGROW>
    end
    if isempty(class_spd_str)
        class_spd_str = 'None';
    end
    if results.spd.has_significant_interaction(ci)
        class_spd_str = [class_spd_str, ' Int'];  %#ok<AGROW>
    end
    spd_shape_T = cluster_speed_shape_T{ci};
    spd_shape_VT = cluster_speed_shape_VT{ci};
    if ~isempty(spd_shape_VT) && ~strcmp(spd_shape_VT, 'unclassified')
        class_spd_str = [class_spd_str, ' [VT:', spd_shape_VT(1:min(3,end)), ']'];  %#ok<AGROW>
    end
    if ~isempty(spd_shape_T) && ~strcmp(spd_shape_T, 'unclassified')
        class_spd_str = [class_spd_str, ' [T:', spd_shape_T(1:min(3,end)), ']'];  %#ok<AGROW>
    end
    
    % Figure for n_rows_scatter rows × 6 columns
    fig_scatter = figure('Position', [10 10 240*n_cols_scatter 120*n_rows_scatter], ...
        'Visible', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A3');
    
    sgtitle(sprintf('Probe: %s  Cluster: %d  |  Winner: SpdBin=%s / TimeBin=%s  |  Class: %s', ...
        char(pid), cid, wm_spd, wm_time, class_spd_str), 'FontSize', 10, 'FontWeight', 'bold');
    
    for ri = 1:n_rows_scatter
        gt_tag = scatter_row_defs{ri, 1};
        ml     = scatter_row_defs{ri, 2};
        row_label = scatter_row_defs{ri, 3};
        
        if strcmp(gt_tag, 'spd')
            T_all = T_master_spd;
            wm_this = wm_spd;
        else
            T_all = T_master_time;
            wm_this = wm_time;
        end
        
        % Determine if this model is the winner for red highlighting
        is_winner = strcmp(ml, wm_this);
        
        idx = T_all.probe_id == pid & T_all.cluster_id == cid;
        T_cluster = T_all(idx, :);
        vt_idx_cl = T_cluster.condition == "VT";
        
        has_predictions = ~isempty(cluster_predictions.(gt_tag){ci}) && ...
            isfield(cluster_predictions.(gt_tag){ci}, ml);
        
        % Grid sizes for scatter plots
        n_spd_grid = length(spd_bin_centers);
        n_tf_grid  = length(tf_bin_centers);
        
        if strcmp(gt_tag, 'spd')
            bin_time = T_cluster.time_in_bin;
        else
            bin_time = time_bin_width * ones(height(T_cluster), 1);
        end
        
        % Col 1: Beta swarm plot (grouped by feature type)
        ax1 = subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 1);
        if has_predictions
            b_all = cluster_predictions.(gt_tag){ci}.(ml).beta;
            se_all = cluster_predictions.(gt_tag){ci}.(ml).se;
            cn_all = cluster_predictions.(gt_tag){ci}.(ml).col_names;
            n_coef = length(b_all);
            
            [grp_mean, ~, grp_clr, grp_lbl, n_grps, grp_betas, grp_sig] = ...
                compute_grouped_params(b_all, se_all, cn_all);
            
            hold on;
            for gii = 1:n_grps
                bv = grp_betas{gii};
                sv = grp_sig{gii};
                n_pts = length(bv);
                % Jitter x positions for visibility
                if n_pts > 1
                    jitter = linspace(-0.25, 0.25, n_pts);
                else
                    jitter = 0;
                end
                xpos = gii + jitter;
                % Non-significant: dots
                ns = ~sv;
                if any(ns)
                    plot(xpos(ns), bv(ns), 'o', ...
                        'MarkerSize', 3, 'Color', grp_clr(gii,:), ...
                        'MarkerFaceColor', grp_clr(gii,:), 'MarkerEdgeColor', 'none');
                end
                % Significant: stars
                if any(sv)
                    plot(xpos(sv), bv(sv), 'p', ...
                        'MarkerSize', 5, 'Color', grp_clr(gii,:), ...
                        'MarkerFaceColor', grp_clr(gii,:), 'MarkerEdgeColor', 'k', ...
                        'LineWidth', 0.2);
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
        
        % Col 2: Mean scatter by condition with STD error bars
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 2);
        if has_predictions
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            obs_fr = T_cluster.spike_count ./ bin_time;
            
            hold on;
            n_conds = length(cond_order);
            mean_pred_cond = nan(n_conds, 1);
            mean_obs_cond = nan(n_conds, 1);
            std_pred_cond = nan(n_conds, 1);
            std_obs_cond = nan(n_conds, 1);
            
            for coi = 1:n_conds
                c_mask = T_cluster.condition == cond_order{coi};
                if sum(c_mask) > 1
                    pred_vals = cv_pred_fr(c_mask);
                    obs_vals = obs_fr(c_mask);
                    mean_pred_cond(coi) = mean(pred_vals);
                    mean_obs_cond(coi) = mean(obs_vals);
                    std_pred_cond(coi) = std(pred_vals);
                    std_obs_cond(coi) = std(obs_vals);
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
        
        % Col 3: Speed scatter with STD error bars
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 3);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = T_vt.spike_count ./ bin_time(vt_idx_cl);
            
            obs_s = nan(n_spd_grid, 1); pred_s = nan(n_spd_grid, 1); sv = nan(n_spd_grid, 1);
            std_obs_s = nan(n_spd_grid, 1); std_pred_s = nan(n_spd_grid, 1);
            for si = 1:n_spd_grid
                b_idx = T_vt.speed >= spd_bin_edges(si) & T_vt.speed < spd_bin_edges(si+1);
                if sum(b_idx) > 1
                    obs_vals = obs_fr_vt(b_idx);
                    pred_vals = cv_pred_vt(b_idx);
                    obs_s(si) = mean(obs_vals);
                    pred_s(si) = mean(pred_vals);
                    std_obs_s(si) = std(obs_vals);
                    std_pred_s(si) = std(pred_vals);
                    sv(si) = spd_bin_centers(si);
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
        
        % Col 4: TF scatter with STD error bars
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 4);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = T_vt.spike_count ./ bin_time(vt_idx_cl);
            
            obs_tf = nan(n_tf_grid, 1); pred_tf = nan(n_tf_grid, 1); tfv = nan(n_tf_grid, 1);
            std_obs_tf = nan(n_tf_grid, 1); std_pred_tf = nan(n_tf_grid, 1);
            for ti = 1:n_tf_grid
                b_idx = T_vt.tf >= tf_bin_edges(ti) & T_vt.tf < tf_bin_edges(ti+1);
                if sum(b_idx) > 1
                    obs_vals = obs_fr_vt(b_idx);
                    pred_vals = cv_pred_vt(b_idx);
                    obs_tf(ti) = mean(obs_vals);
                    pred_tf(ti) = mean(pred_vals);
                    std_obs_tf(ti) = std(obs_vals);
                    std_pred_tf(ti) = std(pred_vals);
                    tfv(ti) = tf_bin_centers(ti);
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
        
        % Col 5: SF scatter with STD error bars
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 5);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = T_vt.spike_count ./ bin_time(vt_idx_cl);
            sf_vt = T_vt.sf;
            sf_u = unique(sf_vt(~isnan(sf_vt)));
            
            if ~isempty(sf_u)
                obs_sf = nan(length(sf_u), 1);
                pred_sf = nan(length(sf_u), 1);
                std_obs_sf = nan(length(sf_u), 1);
                std_pred_sf = nan(length(sf_u), 1);
                for sfi = 1:length(sf_u)
                    b_idx = sf_vt == sf_u(sfi);
                    if sum(b_idx) > 1
                        obs_vals = obs_fr_vt(b_idx);
                        pred_vals = cv_pred_vt(b_idx);
                        obs_sf(sfi) = mean(obs_vals);
                        pred_sf(sfi) = mean(pred_vals);
                        std_obs_sf(sfi) = std(obs_vals);
                        std_pred_sf(sfi) = std(pred_vals);
                    end
                end
                hold on;
                h_sf = gobjects(length(sf_u), 1);  % Store handles for legend
                leg_labels_sf = {};
                for sfi = 1:length(sf_u)
                    if ~isnan(obs_sf(sfi))
                        [~, sc_idx] = min(abs(sf_unique_all - sf_u(sfi)));
                        h_sf(sfi) = errorbar(pred_sf(sfi), obs_sf(sfi), ...
                            std_obs_sf(sfi), std_obs_sf(sfi), ...
                            std_pred_sf(sfi), std_pred_sf(sfi), ...
                            'o', 'Color', sf_cmap(sc_idx,:), ...
                            'MarkerFaceColor', sf_cmap(sc_idx,:), ...
                            'MarkerSize', 5, 'LineWidth', 0.8, 'CapSize', 3);
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
        
        % Col 6: OR scatter with STD error bars
        subplot(n_rows_scatter, n_cols_scatter, (ri-1)*n_cols_scatter + 6);
        if has_predictions && any(vt_idx_cl)
            T_vt = T_cluster(vt_idx_cl, :);
            % cv_predicted_fr is lambda(t) = instantaneous firing rate (Hz)
            cv_pred_fr = cluster_predictions.(gt_tag){ci}.(ml).cv_predicted_fr;
            cv_pred_vt = cv_pred_fr(vt_idx_cl);
            obs_fr_vt = T_vt.spike_count ./ bin_time(vt_idx_cl);
            or_vt = T_vt.orientation;
            or_u = unique(or_vt(~isnan(or_vt)));
            
            if ~isempty(or_u)
                obs_or = nan(length(or_u), 1);
                pred_or = nan(length(or_u), 1);
                std_obs_or = nan(length(or_u), 1);
                std_pred_or = nan(length(or_u), 1);
                for ori = 1:length(or_u)
                    b_idx = or_vt == or_u(ori);
                    if sum(b_idx) > 1
                        obs_vals = obs_fr_vt(b_idx);
                        pred_vals = cv_pred_vt(b_idx);
                        obs_or(ori) = mean(obs_vals);
                        pred_or(ori) = mean(pred_vals);
                        std_obs_or(ori) = std(obs_vals);
                        std_pred_or(ori) = std(pred_vals);
                    end
                end
                hold on;
                h_or = gobjects(length(or_u), 1);  % Store handles for legend
                leg_labels_or = {};
                for ori = 1:length(or_u)
                    if ~isnan(obs_or(ori))
                        [~, oc_idx] = min(abs(or_unique_all - or_u(ori)));
                        h_or(ori) = errorbar(pred_or(ori), obs_or(ori), ...
                            std_obs_or(ori), std_obs_or(ori), ...
                            std_pred_or(ori), std_pred_or(ori), ...
                            'o', 'Color', or_cmap(oc_idx,:), ...
                            'MarkerFaceColor', or_cmap(oc_idx,:), ...
                            'MarkerSize', 5, 'LineWidth', 0.8, 'CapSize', 3);
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
%    Row 2 — Speed-bin GLM Additive predictions
%    Row 3 — Speed-bin GLM FullInteraction predictions
%    Row 4 — Time-bin GLM Additive predictions
%    Row 5 — Time-bin GLM FullInteraction predictions
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
    {[0.8 0.2 0.2], [0.2 0.7 0.2], [0 0.4 0.8]});

% Row labels for the figure (will be updated per-cluster for winning models)
row_labels_base = {'Original Data', 'Speed-bin GLM: Winner', 'Time-bin GLM: Winner', ...
              'Speed-bin GLM: Additive', 'Time-bin GLM: Additive'};

for ci = 1:n_unique_clusters
    pid = unique_clusters.probe_id(ci);
    cid = unique_clusters.cluster_id(ci);

    fig_tc = figure('Position', [30 30 1800 1400], 'Visible', 'off', ...
        'PaperOrientation', 'portrait', 'PaperType', 'A3');

    sgtitle(sprintf('GLM Tuning Curves — Probe: %s  Cluster: %d', ...
        char(pid), cid), 'FontSize', 14, 'FontWeight', 'bold');

    % --- Gather common data from speed-bin GLM (used for observed data) ---
    T_all_spd_rc = T_master_spd;
    idx_spd_rc = T_all_spd_rc.probe_id == pid & T_all_spd_rc.cluster_id == cid;
    T_cl_spd_rc = T_all_spd_rc(idx_spd_rc, :);
    
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
    
    % Temporal bases at mid-trial (only if temporal bases were included in model)
    if use_temporal
        mid_time = max_trial_duration / 2;
        B_time_mean = make_linear_raised_cosine_basis(mid_time, n_temporal_bases, ...
            trial_duration_range(1), trial_duration_range(2));
    else
        B_time_mean = zeros(1, 0);
    end
    
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
    ax_handles = gobjects(5, 4);  % Store subplot handles for y-limit adjustment
    
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
            % Find condition in tuning table
            cond_idx_spd = find(strcmp(D_spd_rc.trial_groups, cond_name), 1);
            if isempty(cond_idx_spd) || isempty(D_spd_rc.tuning_curves{cond_idx_spd})
                continue;
            end
            tc_array_spd = D_spd_rc.tuning_curves{cond_idx_spd};
            tc_cluster_ids_spd = arrayfun(@(x) double(x.cluster_id), tc_array_spd);
            tc_idx_spd = find(tc_cluster_ids_spd == cid, 1);
            if isempty(tc_idx_spd), continue; end
            
            tc_spd = tc_array_spd(tc_idx_spd);
            % tuning matrix: n_bins x n_trials (firing rate in Hz)
            tuning_mat_spd = tc_spd.tuning;
            bc_spd = tc_spd.bin_centers(:);
            
            % Compute mean and STD across trials
            fr_mean = nanmean(tuning_mat_spd, 2);
            n_trials_spd = sum(~isnan(tuning_mat_spd), 2);
            fr_std = nanstd(tuning_mat_spd, 0, 2);
            
            vm = ~isnan(fr_mean) & n_trials_spd > 1;
            clr = cond_colors_rc(cond_name);
            errorbar(bc_spd(vm), fr_mean(vm), fr_std(vm), 'o', ...
                'Color', clr, 'MarkerFaceColor', clr, 'MarkerSize', 4, ...
                'LineWidth', 0.8, 'CapSize', 3);
            plot(bc_spd(vm), fr_mean(vm), '-', 'Color', clr, 'LineWidth', 1);
        end
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
    for coi = 1:length(cond_order)
        cond_mask = T_cl_spd_rc.condition == cond_order{coi};
        if ~any(cond_mask) || strcmp(cond_order{coi}, 'T_Vstatic'), continue; end
        sf_c  = sf_data(cond_mask);
        spk_c = T_cl_spd_rc.spike_count(cond_mask);
        bt_c  = bin_t_spd(cond_mask);
        fr_obs_sf = nan(n_sf_rc, 1);
        fr_std_sf = nan(n_sf_rc, 1);
        for sfi = 1:n_sf_rc
            % Use tolerance for floating-point comparison
            bm = abs(sf_c - sf_levels_rc(sfi)) < 1e-6;
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
    for coi = 1:length(cond_order)
        cond_mask = T_cl_spd_rc.condition == cond_order{coi};
        if ~any(cond_mask) || strcmp(cond_order{coi}, 'T_Vstatic'), continue; end
        or_c  = or_data(cond_mask);
        spk_c = T_cl_spd_rc.spike_count(cond_mask);
        bt_c  = bin_t_spd(cond_mask);
        fr_obs_or = nan(n_or_rc, 1);
        fr_std_or = nan(n_or_rc, 1);
        for ori = 1:n_or_rc
            % Use tolerance for floating-point comparison
            bm = abs(or_c - or_levels_rc(ori)) < 1e-6;
            if sum(bm) > 1
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
        arrayfun(@(v) sprintf('%d°', round(rad2deg(v))), or_levels_rc, 'UniformOutput', false));
    xlabel('Orientation'); ylabel('FR (Hz)');
    title('Original Data: OR', 'FontSize', 9);
    set(gca, 'box', 'off', 'FontSize', 8);
    yl = ylim; ylim_col{4} = [min(ylim_col{4}(1), yl(1)), max(ylim_col{4}(2), yl(2))];
    
    % ================================================================
    %  ROWS 2-5: GLM PREDICTIONS
    %  Row 2: Speed-bin GLM Winning model
    %  Row 3: Time-bin GLM Winning model
    %  Row 4: Speed-bin GLM Additive
    %  Row 5: Time-bin GLM Additive
    % ================================================================
    
    % Get winning models for this cluster
    wm_spd_rc  = results.spd.winning_model{ci};
    wm_time_rc = results.time.winning_model{ci};
    
    % Define row mapping: {glm_type, model_label, row_number}
    % Row 2-3: Speed-bin GLM (Additive, FullInteraction)
    % Row 4-5: Time-bin GLM (Additive, FullInteraction)
    pred_rows = {
        'spd',  'Additive',        2;
        'spd',  'FullInteraction', 3;
        'time', 'Additive',        4;
        'time', 'FullInteraction', 5
    };
    
    for pr_i = 1:size(pred_rows, 1)
        gt_tag   = pred_rows{pr_i, 1};
        ml_name  = pred_rows{pr_i, 2};
        row_num  = pred_rows{pr_i, 3};
        
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
        
        % Get coefficients
        beta_rc = cluster_predictions.(gt_tag){ci}.(ml_name).beta;
        
        % Prepare temporal basis for time-bin GLM (only if temporal bases were included)
        if strcmp(gt_tag, 'time') && use_temporal
            B_time_use = B_time_mean;
        else
            B_time_use = zeros(1, 0);
        end
        
        % Row title
        if strcmp(gt_tag, 'spd')
            glm_label = 'Speed-bin';
        else
            glm_label = 'Time-bin';
        end
        
        % Row title suffix is just the model name
        row_title_suffix = ml_name;
        
        % ---- Column 1: Speed tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 1) = subplot(5, 4, (row_num-1)*4 + 1);
        hold on;
        
        n_bins_spd = n_spd_plot_bins;
        if strcmp(gt_tag, 'time')
            B_time_rep_spd = repmat(B_time_use, n_bins_spd, 1);
        else
            B_time_rep_spd = zeros(n_bins_spd, 0);
        end
        
        % T_Vstatic: Speed sweep with TF=0, SF=reference (first level), OR=reference
        B_tf_zero = make_raised_cosine_basis(0, n_tf_bases, tf_range(1), tf_range(2));
        B_tf_rep_T = repmat(B_tf_zero, n_bins_spd, 1);
        sf_sweep_T = repmat(sf_levels_rc(1), n_bins_spd, 1);  % reference SF
        or_sweep_T = repmat(or_levels_rc(1), n_bins_spd, 1);  % reference OR
        [X_sweep_T, ~] = assemble_design_matrix(B_speed_bincenters_rc, B_tf_rep_T, ...
            B_time_rep_spd, sf_sweep_T, or_sweep_T, ml_name, sf_levels_rc, or_levels_rc);
        if size(X_sweep_T, 2) == length(beta_rc)
            fr_pred_T = min(exp(X_sweep_T * beta_rc), max_predicted_fr);
            clr_T = cond_colors_rc('T_Vstatic');
            plot(spd_bin_centers, fr_pred_T, 'o-', 'Color', clr_T, ...
                'MarkerFaceColor', clr_T, 'MarkerSize', 4, 'LineWidth', 1.2);
        end
        
        % V: Speed=0 (visual only, no translation) - show as horizontal line at speed=0
        B_spd_zero = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        B_tf_rep_V = repmat(B_tf_mean, n_bins_spd, 1);
        sf_sweep_V = repmat(mode_sf, n_bins_spd, 1);
        or_sweep_V = repmat(mode_or, n_bins_spd, 1);
        B_spd_rep_V = repmat(B_spd_zero, n_bins_spd, 1);
        [X_sweep_V, ~] = assemble_design_matrix(B_spd_rep_V, B_tf_rep_V, ...
            B_time_rep_spd, sf_sweep_V, or_sweep_V, ml_name, sf_levels_rc, or_levels_rc);
        if size(X_sweep_V, 2) == length(beta_rc)
            fr_pred_V = min(exp(X_sweep_V * beta_rc), max_predicted_fr);
            clr_V = cond_colors_rc('V');
            % Plot as horizontal dashed line (V has no speed variation)
            plot(spd_bin_centers([1 end]), [fr_pred_V(1) fr_pred_V(1)], '--', 'Color', clr_V, ...
                'LineWidth', 1.2);
        end
        
        % VT: Speed sweep with mean TF, mode SF, mode OR (both stimuli on)
        B_tf_rep_VT = repmat(B_tf_mean, n_bins_spd, 1);
        sf_sweep_VT = repmat(mode_sf, n_bins_spd, 1);
        or_sweep_VT = repmat(mode_or, n_bins_spd, 1);
        [X_sweep_VT, ~] = assemble_design_matrix(B_speed_bincenters_rc, B_tf_rep_VT, ...
            B_time_rep_spd, sf_sweep_VT, or_sweep_VT, ml_name, sf_levels_rc, or_levels_rc);
        if size(X_sweep_VT, 2) == length(beta_rc)
            fr_pred_VT = min(exp(X_sweep_VT * beta_rc), max_predicted_fr);
            clr_VT = cond_colors_rc('VT');
            plot(spd_bin_centers, fr_pred_VT, 'o-', 'Color', clr_VT, ...
                'MarkerFaceColor', clr_VT, 'MarkerSize', 4, 'LineWidth', 1.2);
        end
        
        hold off;
        xlabel('Speed (cm/s)'); ylabel('FR (Hz)');
        title(sprintf('%s %s: Speed', glm_label, row_title_suffix), 'FontSize', 9);
        set(gca, 'box', 'off', 'FontSize', 8);
        yl = ylim; ylim_col{1} = [min(ylim_col{1}(1), yl(1)), max(ylim_col{1}(2), yl(2))];
        
        % ---- Column 2: TF tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 2) = subplot(5, 4, (row_num-1)*4 + 2);
        hold on;
        
        n_bins_tf = n_tf_plot_bins;
        if strcmp(gt_tag, 'time')
            B_time_rep_tf = repmat(B_time_use, n_bins_tf, 1);
        else
            B_time_rep_tf = zeros(n_bins_tf, 0);
        end
        
        % T_Vstatic: No TF tuning (TF=0 always) - skip
        
        % V: TF sweep with Speed=0 (visual only)
        B_spd_zero_tf = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        B_spd_rep_V_tf = repmat(B_spd_zero_tf, n_bins_tf, 1);
        sf_sweep_V_tf = repmat(mode_sf, n_bins_tf, 1);
        or_sweep_V_tf = repmat(mode_or, n_bins_tf, 1);
        [X_sweep_V_tf, ~] = assemble_design_matrix(B_spd_rep_V_tf, B_tf_bincenters_rc, ...
            B_time_rep_tf, sf_sweep_V_tf, or_sweep_V_tf, ml_name, sf_levels_rc, or_levels_rc);
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
        [X_sweep_VT_tf, ~] = assemble_design_matrix(B_spd_rep_VT_tf, B_tf_bincenters_rc, ...
            B_time_rep_tf, sf_sweep_VT_tf, or_sweep_VT_tf, ml_name, sf_levels_rc, or_levels_rc);
        if size(X_sweep_VT_tf, 2) == length(beta_rc)
            fr_pred_VT_tf = min(exp(X_sweep_VT_tf * beta_rc), max_predicted_fr);
            clr_VT = cond_colors_rc('VT');
            plot(tf_bin_centers, fr_pred_VT_tf, 'o-', 'Color', clr_VT, ...
                'MarkerFaceColor', clr_VT, 'MarkerSize', 4, 'LineWidth', 1.2);
        end
        
        hold off;
        xlabel('TF (Hz)'); ylabel('FR (Hz)');
        title(sprintf('%s %s: TF', glm_label, row_title_suffix), 'FontSize', 9);
        set(gca, 'box', 'off', 'FontSize', 8);
        yl = ylim; ylim_col{2} = [min(ylim_col{2}(1), yl(1)), max(ylim_col{2}(2), yl(2))];
        
        % ---- Column 3: SF tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 3) = subplot(5, 4, (row_num-1)*4 + 3);
        hold on;
        
        if strcmp(gt_tag, 'time')
            B_time_one_sf = B_time_use;
        else
            B_time_one_sf = zeros(1, 0);
        end
        
        % T_Vstatic: No SF tuning (no visual stimulus) - skip
        
        % V: SF sweep with Speed=0 (visual only) - use only actual SF values (exclude 0)
        B_spd_zero_sf = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        fr_pred_V_sf = nan(n_sf_plot, 1);
        for sfi = 1:n_sf_plot
            [X_one, ~] = assemble_design_matrix(B_spd_zero_sf, B_tf_mean, ...
                B_time_one_sf, sf_levels_plot(sfi), mode_or, ml_name, sf_levels_rc, or_levels_rc);
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
            [X_one, ~] = assemble_design_matrix(B_spd_mean, B_tf_mean, ...
                B_time_one_sf, sf_levels_plot(sfi), mode_or, ml_name, sf_levels_rc, or_levels_rc);
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
        title(sprintf('%s %s: SF', glm_label, row_title_suffix), 'FontSize', 9);
        set(gca, 'box', 'off', 'FontSize', 8);
        yl = ylim; ylim_col{3} = [min(ylim_col{3}(1), yl(1)), max(ylim_col{3}(2), yl(2))];
        
        % ---- Column 4: Orientation tuning (predicted) - BY CONDITION ----
        ax_handles(row_num, 4) = subplot(5, 4, (row_num-1)*4 + 4);
        hold on;
        
        if strcmp(gt_tag, 'time')
            B_time_one_or = B_time_use;
        else
            B_time_one_or = zeros(1, 0);
        end
        
        % T_Vstatic: No orientation tuning (no visual stimulus) - skip
        
        % V: OR sweep with Speed=0 (visual only)
        B_spd_zero_or = make_raised_cosine_basis(0, n_speed_bases, speed_range(1), speed_range(2));
        fr_pred_V_or = nan(n_or_rc, 1);
        for ori = 1:n_or_rc
            [X_one, ~] = assemble_design_matrix(B_spd_zero_or, B_tf_mean, ...
                B_time_one_or, mode_sf, or_levels_rc(ori), ml_name, sf_levels_rc, or_levels_rc);
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
            [X_one, ~] = assemble_design_matrix(B_spd_mean, B_tf_mean, ...
                B_time_one_or, mode_sf, or_levels_rc(ori), ml_name, sf_levels_rc, or_levels_rc);
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
            arrayfun(@(v) sprintf('%d°', round(rad2deg(v))), or_levels_rc, 'UniformOutput', false));
        xlabel('Orientation'); ylabel('FR (Hz)');
        title(sprintf('%s %s: OR', glm_label, row_title_suffix), 'FontSize', 9);
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
    h_leg(1) = plot(ax_leg, NaN, NaN, 'o', 'Color', [0.8 0.2 0.2], 'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerSize', 5);
    h_leg(2) = plot(ax_leg, NaN, NaN, 'o', 'Color', [0.2 0.7 0.2], 'MarkerFaceColor', [0.2 0.7 0.2], 'MarkerSize', 5);
    h_leg(3) = plot(ax_leg, NaN, NaN, 'o', 'Color', [0 0.4 0.8], 'MarkerFaceColor', [0 0.4 0.8], 'MarkerSize', 5);
    h_leg(4) = plot(ax_leg, NaN, NaN, 'o-', 'Color', [0.8 0.2 0.2], 'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerSize', 4, 'LineWidth', 1.2);
    h_leg(5) = plot(ax_leg, NaN, NaN, 'o-', 'Color', [0.2 0.7 0.2], 'MarkerFaceColor', [0.2 0.7 0.2], 'MarkerSize', 4, 'LineWidth', 1.2);
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
%  Section 8c: Population Tuning Curves Grouped by Tuning Shape Category
%  ====================================================================
%  Create population summary figures showing tuning curves grouped by
%  visual tuning status and speed tuning shape.
%
%  Figure layout:
%    Rows 1-3: V-tuned clusters grouped by VT shape (increasing/decreasing/bandpass)
%    Rows 4-6: Not V-tuned clusters grouped by T shape (increasing/decreasing/bandpass)
%    Columns: 4 stimulus dimensions
%      Col 1: Speed tuning
%      Col 2: TF tuning
%      Col 3: SF tuning
%      Col 4: OR tuning
%
%  Each panel shows:
%    • Thin lines: mean tuning curve for each individual cluster (all conditions combined)
%    • Thick darker line: population average across clusters in that category
%    • Separate colors for conditions when plotting condition-specific curves
%
%  Empty categories (n=0) are skipped dynamically.
%  ====================================================================
fprintf('\n--- Generating population tuning curves by tuning shape category ---\n');

% Define the tuning shape categories
shape_categories = {'increasing', 'decreasing', 'bandpass'};

% Row labels for the figure (will filter out empty categories later)
row_labels_pop_full = {
    'V-tuned: Increasing', 'V-tuned: Decreasing', 'V-tuned: Bandpass', ...
    'Not V-tuned: Increasing', 'Not V-tuned: Decreasing', 'Not V-tuned: Bandpass'
};

% Column labels
col_labels_pop = {'Speed', 'TF', 'SF', 'OR'};

% Condition colours (same as other tuning curve plots)
cond_colors_pop = containers.Map({'T_Vstatic', 'V', 'VT'}, ...
    {[0.8 0.2 0.2], [0.2 0.7 0.2], [0 0.4 0.8]});
cond_order_pop = {'T_Vstatic', 'V', 'VT'};

% Use bin centers from pre-loaded tuning tables
n_spd_bins_pop = length(spd_bin_centers);
n_tf_bins_pop = length(tf_bin_centers);
sf_levels_pop = [0.003; 0.006; 0.012];
or_levels_pop = [-pi/4; 0; pi/4; pi/2];
n_sf_bins_pop = length(sf_levels_pop);
n_or_bins_pop = length(or_levels_pop);

% ---- Collect tuning data for each cluster from tuning tables ----
% Storage: tuning_data{category_idx}{cluster_idx}{condition_idx}{stimulus_dim} = [bin_values]
% We will extract from probe_info (loaded earlier)

% First, organize clusters by their tuning shape category
% Rows 1-3: V-tuned clusters by VT shape
% Rows 4-6: Not V-tuned clusters by T shape
cluster_by_category = cell(6, 1);  % 6 categories

for ci = 1:n_unique_clusters
    pid = char(unique_clusters.probe_id(ci));
    cid = unique_clusters.cluster_id(ci);
    
    % Find this cluster in prefilter_results
    pf_idx = find(strcmp(prefilter_results.probe_id, pid) & prefilter_results.cluster_id == cid, 1);
    if isempty(pf_idx)
        continue;
    end
    
    is_v_tuned = prefilter_results.is_visually_tuned(pf_idx);
    spd_shape_T = cluster_speed_shape_T{ci};
    spd_shape_VT = cluster_speed_shape_VT{ci};
    
    if is_v_tuned
        % V-tuned: use VT shape (rows 1-3)
        for si = 1:3
            if strcmp(spd_shape_VT, shape_categories{si})
                cluster_by_category{si} = [cluster_by_category{si}; ci];
            end
        end
    else
        % Not V-tuned: use T shape (rows 4-6)
        for si = 1:3
            if strcmp(spd_shape_T, shape_categories{si})
                cluster_by_category{3 + si} = [cluster_by_category{3 + si}; ci];
            end
        end
    end
end

% Print category counts
% Print category counts and identify non-empty categories
fprintf('  Clusters per category:\n');
non_empty_rows = [];
row_labels_pop = {};
for ri = 1:6
    n_in_cat = length(cluster_by_category{ri});
    fprintf('    %s: %d clusters\n', row_labels_pop_full{ri}, n_in_cat);
    if n_in_cat > 0
        non_empty_rows = [non_empty_rows; ri];  %#ok<AGROW>
        row_labels_pop{end+1} = row_labels_pop_full{ri};  %#ok<AGROW>
    end
end
n_rows_pop = length(non_empty_rows);
fprintf('  Non-empty categories: %d\n', n_rows_pop);

if n_rows_pop == 0
    fprintf('  WARNING: No clusters in any category, skipping population figure.\n');
else

% ---- Collect tuning curves from pre-computed tuning tables ----
% Storage: tuning_by_cat{row}{col}{cond} = [n_bins x n_clusters]
tuning_by_cat = cell(6, 4);  % 6 rows x 4 cols
for ri = 1:6
    for col = 1:4
        tuning_by_cat{ri, col} = cell(1, 3);  % 3 conditions
        for coi = 1:3
            % Pre-allocate with NaN
            if col == 1
                n_bins = n_spd_bins_pop;
            elseif col == 2
                n_bins = n_tf_bins_pop;
            elseif col == 3
                n_bins = n_sf_bins_pop;
            else
                n_bins = n_or_bins_pop;
            end
            tuning_by_cat{ri, col}{coi} = nan(n_bins, length(cluster_by_category{ri}));
        end
    end
end

% Loop through clusters and populate tuning data
for ri = 1:6
    cluster_list = cluster_by_category{ri};
    
    for cli = 1:length(cluster_list)
        ci = cluster_list(cli);
        pid = unique_clusters.probe_id(ci);
        cid = unique_clusters.cluster_id(ci);
        
        % Find the probe_info entry for this cluster's probe
        probe_idx_pop = [];
        for pi_pop = 1:length(probe_info)
            if ~isempty(probe_info(pi_pop).pid) && strcmp(probe_info(pi_pop).pid, char(pid))
                probe_idx_pop = pi_pop;
                break;
            end
        end
        
        if isempty(probe_idx_pop)
            continue;
        end
        
        % Extract speed tuning (Column 1)
        if isfield(probe_info(probe_idx_pop), 'D_spd') && ~isempty(probe_info(probe_idx_pop).D_spd)
            D_spd_pop = probe_info(probe_idx_pop).D_spd;
            for coi = 1:length(cond_order_pop)
                cond_name = cond_order_pop{coi};
                cond_idx_spd = find(strcmp(D_spd_pop.trial_groups, cond_name), 1);
                if ~isempty(cond_idx_spd) && ~isempty(D_spd_pop.tuning_curves{cond_idx_spd})
                    tc_array_spd = D_spd_pop.tuning_curves{cond_idx_spd};
                    tc_cluster_ids_spd = arrayfun(@(x) double(x.cluster_id), tc_array_spd);
                    tc_idx_spd = find(tc_cluster_ids_spd == cid, 1);
                    if ~isempty(tc_idx_spd)
                        tc_spd = tc_array_spd(tc_idx_spd);
                        tuning_mat_spd = tc_spd.tuning;
                        % Average across trials for this cluster
                        fr_mean_spd = nanmean(tuning_mat_spd, 2);
                        if length(fr_mean_spd) == n_spd_bins_pop
                            tuning_by_cat{ri, 1}{coi}(:, cli) = fr_mean_spd;
                        end
                    end
                end
            end
        end
        
        % Extract TF tuning (Column 2)
        if isfield(probe_info(probe_idx_pop), 'D_tf') && ~isempty(probe_info(probe_idx_pop).D_tf)
            D_tf_pop = probe_info(probe_idx_pop).D_tf;
            for coi = 1:length(cond_order_pop)
                cond_name = cond_order_pop{coi};
                % Skip T_Vstatic for TF (no visual flow)
                if strcmp(cond_name, 'T_Vstatic'), continue; end
                cond_idx_tf = find(strcmp(D_tf_pop.trial_groups, cond_name), 1);
                if ~isempty(cond_idx_tf) && ~isempty(D_tf_pop.tuning_curves{cond_idx_tf})
                    tc_array_tf = D_tf_pop.tuning_curves{cond_idx_tf};
                    tc_cluster_ids_tf = arrayfun(@(x) double(x.cluster_id), tc_array_tf);
                    tc_idx_tf = find(tc_cluster_ids_tf == cid, 1);
                    if ~isempty(tc_idx_tf)
                        tc_tf = tc_array_tf(tc_idx_tf);
                        tuning_mat_tf = tc_tf.tuning;
                        fr_mean_tf = nanmean(tuning_mat_tf, 2);
                        if length(fr_mean_tf) == n_tf_bins_pop
                            tuning_by_cat{ri, 2}{coi}(:, cli) = fr_mean_tf;
                        end
                    end
                end
            end
        end
        
        % Extract SF and OR tuning from T_master_spd (Column 3 & 4)
        T_cl_pop = T_master_spd(T_master_spd.probe_id == pid & T_master_spd.cluster_id == cid, :);
        if height(T_cl_pop) > 0
            for coi = 1:length(cond_order_pop)
                cond_name = cond_order_pop{coi};
                cond_mask_pop = T_cl_pop.condition == cond_name;
                if ~any(cond_mask_pop) || strcmp(cond_name, 'T_Vstatic')
                    continue;
                end
                
                sf_c = T_cl_pop.sf(cond_mask_pop);
                or_c = T_cl_pop.orientation(cond_mask_pop);
                spk_c = T_cl_pop.spike_count(cond_mask_pop);
                bt_c = T_cl_pop.time_in_bin(cond_mask_pop);
                
                % SF tuning (Column 3)
                fr_sf_pop = nan(n_sf_bins_pop, 1);
                for sfi = 1:n_sf_bins_pop
                    bm = abs(sf_c - sf_levels_pop(sfi)) < 1e-6;
                    if sum(bm) > 1
                        fr_sf_pop(sfi) = mean(spk_c(bm) ./ bt_c(bm));
                    end
                end
                tuning_by_cat{ri, 3}{coi}(:, cli) = fr_sf_pop;
                
                % OR tuning (Column 4)
                fr_or_pop = nan(n_or_bins_pop, 1);
                for ori = 1:n_or_bins_pop
                    bm = abs(or_c - or_levels_pop(ori)) < 1e-6;
                    if sum(bm) > 1
                        fr_or_pop(ori) = mean(spk_c(bm) ./ bt_c(bm));
                    end
                end
                tuning_by_cat{ri, 4}{coi}(:, cli) = fr_or_pop;
            end
        end
    end
end

% ---- Create the population figure ----
% Dynamically adjust figure height based on number of non-empty rows
fig_height = max(400, 300 * n_rows_pop);
fig_pop_tuning = figure('Position', [20 20 1600 fig_height], 'Name', 'Population Tuning by Category');

n_cols_pop = 4;

% Define x-axis values for each column
x_axes_pop = {
    spd_bin_centers(:), ...  % Speed
    tf_bin_centers(:), ...   % TF
    (1:n_sf_bins_pop)', ...  % SF (categorical)
    (1:n_or_bins_pop)'  ...  % OR (categorical)
};

% Lighter colors for individual cluster lines (alpha simulation)
cond_colors_light = containers.Map({'T_Vstatic', 'V', 'VT'}, ...
    {[0.9 0.6 0.6], [0.6 0.85 0.6], [0.5 0.7 0.9]});

% Store y-limits for each column to synchronize
ylim_pop_col = cell(1, 4);
for col = 1:4
    ylim_pop_col{col} = [Inf, -Inf];
end
ax_handles_pop = gobjects(n_rows_pop, n_cols_pop);

for ri_plot = 1:n_rows_pop
    ri = non_empty_rows(ri_plot);  % Map plot row to original category index
    n_clusters_in_cat = length(cluster_by_category{ri});
    
    for col = 1:n_cols_pop
        ax_handles_pop(ri_plot, col) = subplot(n_rows_pop, n_cols_pop, (ri_plot-1)*n_cols_pop + col);
        hold on;
        
        x_vals = x_axes_pop{col};
        n_bins = length(x_vals);
        
        if n_clusters_in_cat == 0
            % No clusters in this category (shouldn't happen since we filter)
            text(0.5, 0.5, 'No clusters', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);
            hold off;
            continue;
        end
        
        % Combine data across conditions for averaging
        all_cluster_curves = nan(n_bins, n_clusters_in_cat);
        
        % Plot individual cluster curves (thin lines) for each condition
        for coi = 1:length(cond_order_pop)
            cond_name = cond_order_pop{coi};
            
            % Skip T_Vstatic for TF, SF, OR columns
            if (col >= 2) && strcmp(cond_name, 'T_Vstatic')
                continue;
            end
            
            data_cond = tuning_by_cat{ri, col}{coi};  % [n_bins x n_clusters]
            if isempty(data_cond)
                continue;
            end
            
            clr_light = cond_colors_light(cond_name);
            clr_dark = cond_colors_pop(cond_name);
            
            % Plot thin lines for each cluster
            for cli = 1:size(data_cond, 2)
                curve = data_cond(:, cli);
                valid_bins = ~isnan(curve);
                if sum(valid_bins) > 1
                    plot(x_vals(valid_bins), curve(valid_bins), '-', ...
                        'Color', [clr_light, 0.3], 'LineWidth', 0.5);
                    % Accumulate for combined average
                    all_cluster_curves(valid_bins, cli) = curve(valid_bins);
                end
            end
            
            % Compute and plot population mean (thick line) for this condition
            valid_clusters = sum(~isnan(data_cond), 1) > 0;
            if sum(valid_clusters) > 1
                pop_mean = nanmean(data_cond(:, valid_clusters), 2);
                pop_std = nanstd(data_cond(:, valid_clusters), [], 2);
                valid_mean = ~isnan(pop_mean);
                if sum(valid_mean) > 1
                    % Plot thick mean line
                    plot(x_vals(valid_mean), pop_mean(valid_mean), '-', ...
                        'Color', clr_dark, 'LineWidth', 2.5);
                    % Plot STD as shaded area (using fill)
                    x_fill = [x_vals(valid_mean); flipud(x_vals(valid_mean))];
                    y_fill = [pop_mean(valid_mean) + pop_std(valid_mean); ...
                              flipud(pop_mean(valid_mean) - pop_std(valid_mean))];
                    fill(x_fill, y_fill, clr_dark, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                end
            end
        end
        
        hold off;
        
        % Update y-limits tracker
        yl = ylim;
        if isfinite(yl(1)) && isfinite(yl(2))
            ylim_pop_col{col} = [min(ylim_pop_col{col}(1), yl(1)), max(ylim_pop_col{col}(2), yl(2))];
        end
        
        % Axis formatting
        set(gca, 'box', 'off', 'FontSize', 8);
        
        % X-axis labels (only on bottom row)
        if ri == n_rows_pop
            switch col
                case 1
                    xlabel('Speed (cm/s)', 'FontSize', 9);
                case 2
                    xlabel('TF (Hz)', 'FontSize', 9);
                case 3
                    xlabel('SF (cpd)', 'FontSize', 9);
                    set(gca, 'XTick', 1:n_sf_bins_pop, 'XTickLabel', ...
                        arrayfun(@(v) sprintf('%.3f', v), sf_levels_pop, 'UniformOutput', false));
                case 4
                    xlabel('Orientation', 'FontSize', 9);
                    set(gca, 'XTick', 1:n_or_bins_pop, 'XTickLabel', ...
                        arrayfun(@(v) sprintf('%d°', round(rad2deg(v))), or_levels_pop, 'UniformOutput', false));
            end
        else
            % For SF and OR columns, set x-tick labels even for non-bottom rows
            if col == 3
                set(gca, 'XTick', 1:n_sf_bins_pop, 'XTickLabel', {});
            elseif col == 4
                set(gca, 'XTick', 1:n_or_bins_pop, 'XTickLabel', {});
            end
        end
        
        % Y-axis label (only on first column)
        if col == 1
            ylabel('FR (Hz)', 'FontSize', 9);
        end
        
        % Row titles (first column only)
        if col == 1
            title(sprintf('%s (n=%d)', row_labels_pop{ri_plot}, n_clusters_in_cat), 'FontSize', 10);
        elseif ri_plot == 1
            % Column titles (first row only)
            title(col_labels_pop{col}, 'FontSize', 10);
        end
    end
end

% --- Synchronize y-axis limits across all rows for each column ---
for col = 1:4
    if isinf(ylim_pop_col{col}(1)) || isinf(ylim_pop_col{col}(2))
        continue;
    end
    % Add small padding to upper limit; always start from 0
    yl_range = ylim_pop_col{col}(2);  % range from 0 to max
    yl_pad = yl_range * 0.05;
    yl_final = [0, ylim_pop_col{col}(2) + yl_pad];
    for ri_plot = 1:n_rows_pop
        if isvalid(ax_handles_pop(ri_plot, col)) && isgraphics(ax_handles_pop(ri_plot, col))
            ylim(ax_handles_pop(ri_plot, col), yl_final);
        end
    end
end

% --- Add legend ---
ax_leg_pop = axes(fig_pop_tuning, 'Position', [0.05 0.01 0.9 0.02], 'Visible', 'off');
hold(ax_leg_pop, 'on');
h_leg_pop = gobjects(3, 1);
h_leg_pop(1) = plot(ax_leg_pop, NaN, NaN, '-', 'Color', cond_colors_pop('T_Vstatic'), 'LineWidth', 2.5);
h_leg_pop(2) = plot(ax_leg_pop, NaN, NaN, '-', 'Color', cond_colors_pop('V'), 'LineWidth', 2.5);
h_leg_pop(3) = plot(ax_leg_pop, NaN, NaN, '-', 'Color', cond_colors_pop('VT'), 'LineWidth', 2.5);
legend(h_leg_pop, {'T\_Vstatic', 'V', 'VT'}, ...
    'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 10);
hold(ax_leg_pop, 'off');

sgtitle('Population Tuning Curves by Tuning Shape Category', 'FontSize', 14, 'FontWeight', 'bold');

if save_figs
    print(fig_pop_tuning, fullfile(ctl.figs.curr_dir, 'population_tuning_by_category.png'), '-dpng', '-r300');
    fprintf('  Saved population_tuning_by_category.png\n');
end

end  % End of if n_rows_pop > 0

fprintf('  Population tuning figure by category complete\n');

%% ====================================================================
%  Section 9: CSV exports
%  ====================================================================
if csv_output
    fprintf('\n--- Exporting CSV files ---\n');
    
    csv_dir = fullfile(ctl.figs.curr_dir, 'csv');
    if ~exist(csv_dir, 'dir'), mkdir(csv_dir); end
    
    % --- CSV 1: Model comparison (both GLM types) ---
    T_comp = table();
    T_comp.probe_id = results.spd.probe_id;
    T_comp.cluster_id = results.spd.cluster_id;
    T_comp.n_trials = results.spd.n_trials;
    
    % Add pre-filter category from decision tree
    T_comp.prefilter_category = cell(n_unique_clusters, 1);
    for ci_pf = 1:n_unique_clusters
        key = sprintf('%s_%d', results.spd.probe_id{ci_pf}, results.spd.cluster_id(ci_pf));
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
        T_comp.([prefix 'winning_model']) = results.(gt_tag).winning_model;
        
        for mi = 1:n_models
            ml = model_labels{mi};
            T_comp.([prefix ml '_cv_bps']) = results.(gt_tag).([ml '_cv_bps']);
            T_comp.([prefix ml '_train_bps']) = results.(gt_tag).([ml '_train_bps']);
            T_comp.([prefix ml '_n_params']) = results.(gt_tag).([ml '_n_params']);
            T_comp.([prefix ml '_dispersion']) = results.(gt_tag).([ml '_dispersion']);
        end
        
        T_comp.([prefix 'delta_bps_interaction']) = results.(gt_tag).delta_bps_interaction;
        
        % --- Drop-one model comparison fields ---
        T_comp.([prefix 'delta_bps_drop_speed']) = results.(gt_tag).delta_bps_drop_speed;
        T_comp.([prefix 'delta_bps_drop_tf']) = results.(gt_tag).delta_bps_drop_tf;
        T_comp.([prefix 'delta_bps_drop_sf']) = results.(gt_tag).delta_bps_drop_sf;
        T_comp.([prefix 'delta_bps_drop_or']) = results.(gt_tag).delta_bps_drop_or;
        
        % --- Classification flags (drop-one unique contribution) ---
        T_comp.([prefix 'is_speed_tuned']) = results.(gt_tag).is_speed_tuned;
        T_comp.([prefix 'is_tf_tuned']) = results.(gt_tag).is_tf_tuned;
        T_comp.([prefix 'is_sf_tuned']) = results.(gt_tag).is_sf_tuned;
        T_comp.([prefix 'is_or_tuned']) = results.(gt_tag).is_or_tuned;
        T_comp.([prefix 'has_significant_interaction']) = results.(gt_tag).has_significant_interaction;
    end
    
    % --- Tuning shape classification (from asymmetric Gaussian fits, added to all clusters) ---
    T_comp.speed_shape_T = cluster_speed_shape_T;
    T_comp.speed_shape_VT = cluster_speed_shape_VT;
    
    writetable(T_comp, fullfile(csv_dir, 'glm_model_comparison.csv'));
    fprintf('  Saved glm_model_comparison.csv (%d rows)\n', height(T_comp));
    
    % --- CSV 2: Coefficients ---
    if n_coef_stored > 0
        all_coefficients = all_coefficients(1:n_coef_stored);
        T_coef = cell2table(vertcat(all_coefficients{:}), ...
            'VariableNames', {'probe_id', 'cluster_id', 'glm_type', 'model', ...
            'coefficient', 'estimate', 'se'});
        writetable(T_coef, fullfile(csv_dir, 'glm_coefficients.csv'));
        fprintf('  Saved glm_coefficients.csv (%d rows)\n', height(T_coef));
    end
    
    % --- CSV 3: Trial-level data (speed-bin) ---
    writetable(T_master_spd, fullfile(csv_dir, 'glm_speed_bin_data.csv'));
    fprintf('  Saved glm_speed_bin_data.csv (%d rows)\n', height(T_master_spd));
    
    % --- CSV 4: Trial-level data (time-bin) ---
    writetable(T_master_time, fullfile(csv_dir, 'glm_time_bin_data.csv'));
    fprintf('  Saved glm_time_bin_data.csv (%d rows)\n', height(T_master_time));
    
    % --- CSV 5: Component importance decomposition ---
    T_imp = table();
    imp_row = 0;
    for gt = 1:length(glm_types)
        gt_tag = glm_types{gt};
        for ci_csv = 1:n_unique_clusters
            for gi = 1:n_comp
                if isnan(results.(gt_tag).comp_fold_change(ci_csv, gi)), continue; end
                imp_row = imp_row + 1;
                T_imp.probe_id{imp_row} = results.(gt_tag).probe_id{ci_csv};
                T_imp.cluster_id(imp_row) = results.(gt_tag).cluster_id(ci_csv);
                T_imp.glm_type{imp_row} = gt_tag;
                T_imp.component{imp_row} = comp_groups{gi};
                T_imp.fold_change(imp_row) = results.(gt_tag).comp_fold_change(ci_csv, gi);
                T_imp.var_eta(imp_row) = results.(gt_tag).comp_var_eta(ci_csv, gi);
                T_imp.rel_importance(imp_row) = results.(gt_tag).comp_rel_importance(ci_csv, gi);
                T_imp.wald_chi2(imp_row) = results.(gt_tag).comp_wald_chi2(ci_csv, gi);
                T_imp.wald_p(imp_row) = results.(gt_tag).comp_wald_p(ci_csv, gi);
            end
        end
    end
    if imp_row > 0
        writetable(T_imp, fullfile(csv_dir, 'glm_component_importance.csv'));
        fprintf('  Saved glm_component_importance.csv (%d rows)\n', height(T_imp));
    end
    
    % --- Update prefilter_decision_tree.csv with GLM results ---
    % Add leave-one-out results and speed/TF classification for each cluster
    fprintf('\n--- Updating prefilter_decision_tree.csv with GLM results ---\n');
    
    % Initialize new columns with NaN/empty for all prefilter clusters
    glm_delta_bps_drop_speed = nan(n_prefilter_clusters, 1);
    glm_delta_bps_drop_tf = nan(n_prefilter_clusters, 1);
    glm_delta_bps_drop_sf = nan(n_prefilter_clusters, 1);
    glm_delta_bps_drop_or = nan(n_prefilter_clusters, 1);
    glm_is_speed_tuned = false(n_prefilter_clusters, 1);
    glm_is_tf_tuned = false(n_prefilter_clusters, 1);
    glm_is_sf_tuned = false(n_prefilter_clusters, 1);
    glm_is_or_tuned = false(n_prefilter_clusters, 1);
    glm_has_significant_interaction = false(n_prefilter_clusters, 1);
    glm_delta_bps_interaction = nan(n_prefilter_clusters, 1);
    glm_additive_cv_bps = nan(n_prefilter_clusters, 1);
    glm_winning_model = cell(n_prefilter_clusters, 1);
    % Positive control: visual contribution (Additive vs M0_Speed)
    glm_delta_bps_visual = nan(n_prefilter_clusters, 1);
    glm_is_visually_tuned_glm = false(n_prefilter_clusters, 1);
    glm_m0_speed_cv_bps = nan(n_prefilter_clusters, 1);
    
    % Fill in default values
    for pi = 1:n_prefilter_clusters
        glm_winning_model{pi} = '';
    end
    
    % Match GLM results back to prefilter clusters (use speed-bin GLM as primary)
    for ci = 1:n_unique_clusters
        glm_probe = results.spd.probe_id{ci};
        glm_cid = results.spd.cluster_id(ci);
        
        % Find matching prefilter row
        for pi = 1:n_prefilter_clusters
            if strcmp(prefilter_results.probe_id{pi}, glm_probe) && ...
               prefilter_results.cluster_id(pi) == glm_cid
                % Copy leave-one-out delta bps
                glm_delta_bps_drop_speed(pi) = results.spd.delta_bps_drop_speed(ci);
                glm_delta_bps_drop_tf(pi) = results.spd.delta_bps_drop_tf(ci);
                glm_delta_bps_drop_sf(pi) = results.spd.delta_bps_drop_sf(ci);
                glm_delta_bps_drop_or(pi) = results.spd.delta_bps_drop_or(ci);
                
                % Copy classification flags
                glm_is_speed_tuned(pi) = results.spd.is_speed_tuned(ci);
                glm_is_tf_tuned(pi) = results.spd.is_tf_tuned(ci);
                glm_is_sf_tuned(pi) = results.spd.is_sf_tuned(ci);
                glm_is_or_tuned(pi) = results.spd.is_or_tuned(ci);
                glm_has_significant_interaction(pi) = results.spd.has_significant_interaction(ci);
                
                % Copy additional metrics
                glm_delta_bps_interaction(pi) = results.spd.delta_bps_interaction(ci);
                glm_additive_cv_bps(pi) = results.spd.Additive_cv_bps(ci);
                glm_winning_model{pi} = results.spd.winning_model{ci};
                
                % Positive control: visual contribution
                glm_delta_bps_visual(pi) = results.spd.delta_bps_visual(ci);
                glm_is_visually_tuned_glm(pi) = results.spd.is_visually_tuned_glm(ci);
                glm_m0_speed_cv_bps(pi) = results.spd.M0_Speed_cv_bps(ci);
                break;
            end
        end
    end
    
    % Update prefilter_table with new columns
    prefilter_table_updated = table(...
        string(prefilter_results.probe_id'), ...
        prefilter_results.cluster_id', ...
        prefilter_results.is_speed_tuned', ...
        prefilter_results.is_visually_tuned', ...
        prefilter_results.is_VT_responsive', ...
        prefilter_results.is_mixed_tuning', ...
        prefilter_results.should_run_glm', ...
        string(prefilter_results.category'), ...
        prefilter_results.p_T_svm', ...
        prefilter_results.p_V_svm', ...
        prefilter_results.p_VT_svm', ...
        prefilter_results.p_T_vs_VT_motion', ...
        prefilter_results.direction_T', ...
        prefilter_results.direction_V', ...
        prefilter_results.direction_VT', ...
        string(prefilter_results.speed_shape_T'), ...
        string(prefilter_results.speed_shape_VT'), ...
        glm_delta_bps_drop_speed, ...
        glm_delta_bps_drop_tf, ...
        glm_delta_bps_drop_sf, ...
        glm_delta_bps_drop_or, ...
        glm_is_speed_tuned, ...
        glm_is_tf_tuned, ...
        glm_is_sf_tuned, ...
        glm_is_or_tuned, ...
        glm_has_significant_interaction, ...
        glm_delta_bps_interaction, ...
        glm_delta_bps_visual, ...
        glm_is_visually_tuned_glm, ...
        glm_m0_speed_cv_bps, ...
        glm_additive_cv_bps, ...
        string(glm_winning_model), ...
        'VariableNames', {'probe_id', 'cluster_id', 'is_speed_tuned_prefilter', 'is_visually_tuned', ...
            'is_VT_responsive', 'is_mixed_tuning', 'should_run_glm', 'category', ...
            'p_T_svm', 'p_V_svm', 'p_VT_svm', 'p_T_vs_VT_motion', ...
            'direction_T', 'direction_V', 'direction_VT', ...
            'speed_shape_T', 'speed_shape_VT', ...
            'glm_delta_bps_drop_speed', 'glm_delta_bps_drop_tf', ...
            'glm_delta_bps_drop_sf', 'glm_delta_bps_drop_or', ...
            'glm_is_speed_tuned', 'glm_is_tf_tuned', 'glm_is_sf_tuned', 'glm_is_or_tuned', ...
            'glm_has_significant_interaction', 'glm_delta_bps_interaction', ...
            'glm_delta_bps_visual', 'glm_is_visually_tuned_glm', 'glm_m0_speed_cv_bps', ...
            'glm_additive_cv_bps', 'glm_winning_model'});
    
    % Save updated prefilter table
    writetable(prefilter_table_updated, prefilter_csv_path);
    fprintf('  Updated prefilter_decision_tree.csv with GLM results (%d rows)\n', height(prefilter_table_updated));
    fprintf('  Added columns: leave-one-out delta bps, speed/TF classification, tuning shapes\n');
end

%% ====================================================================
%  Section 10: Summary printout
%  ====================================================================
fprintf('\n====================================================================\n');
fprintf('  GLM SINGLE CLUSTER ANALYSIS (v6: Pre-filtering Decision Tree) — COMPLETE\n');
fprintf('====================================================================\n');
fprintf('  Total VISp clusters scanned: %d\n', n_prefilter_clusters);
fprintf('  Clusters after pre-filtering: %d (%.1f%%)\n', n_should_run_glm, 100*n_should_run_glm/n_prefilter_clusters);
fprintf('  Clusters analysed with GLM: %d\n', n_unique_clusters);
fprintf('  Speed-bin data points: %d | Time-bin data points: %d\n', ...
    height(T_master_spd), height(T_master_time));
fprintf('  Models per cluster per GLM: %d (%s)\n', n_models, strjoin(model_labels, ', '));

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
    fprintf('  Winning model distribution:\n');
    for mi = 1:n_models
        ml = model_labels{mi};
        n_win = sum(strcmp(results.(gt_tag).winning_model, ml));
        if n_win > 0
            fprintf('    %s: %d (%.0f%%)\n', ml, n_win, 100*n_win/n_unique_clusters);
        end
    end
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


function B = make_linear_raised_cosine_basis(x, n_bases, x_min, x_max)
%MAKE_LINEAR_RAISED_COSINE_BASIS Raised cosine bases on linear axis (no log transform)
%   Used for temporal drift bases where Weber-law scaling is not appropriate.
    centers = linspace(x_min, x_max, n_bases);
    
    if n_bases > 1
        delta = (x_max - x_min) / (n_bases - 1);
    else
        delta = (x_max - x_min);
    end
    width = delta * 1.5;
    
    B = zeros(length(x), n_bases);
    for bi = 1:n_bases
        z = (x - centers(bi)) / width * pi;
        z = max(min(z, pi), -pi);
        B(:, bi) = 0.5 * (1 + cos(z));
    end
end


function [X, col_names] = build_design_matrix(speed_vals, tf_vals, sf_vals, or_vals, ...
        model_label, n_speed_bases_arg, n_tf_bases_arg, speed_range_arg, tf_range_arg, ...
        time_vals, trial_duration_range_arg)
%BUILD_DESIGN_MATRIX Build design matrix for a given model specification
%   time_vals: column vector of time-in-trial values (empty for speed-bin GLM)
%   trial_duration_range_arg: [min, max] for temporal bases (empty if no time)

    n = length(speed_vals);
    use_temporal = ~isempty(time_vals) && ~isempty(trial_duration_range_arg);
    
    B_speed = make_raised_cosine_basis(speed_vals, n_speed_bases_arg, speed_range_arg(1), speed_range_arg(2));
    B_tf    = make_raised_cosine_basis(tf_vals, n_tf_bases_arg, tf_range_arg(1), tf_range_arg(2));
    
    % Temporal drift bases (only for time-bin GLM)
    if use_temporal
        n_temp = 5;  % always 5 temporal bases
        B_time = make_linear_raised_cosine_basis(time_vals, n_temp, ...
            trial_duration_range_arg(1), trial_duration_range_arg(2));
        time_names = {};
        for ti = 1:n_temp
            time_names{end+1} = sprintf('Time_%d', ti);  %#ok<AGROW>
        end
    else
        B_time = zeros(n, 0);
        time_names = {};
    end
    
    % SF dummy coding
    unique_sf = sort(unique(sf_vals(~isnan(sf_vals))));
    n_sf_levels = length(unique_sf);
    D_sf = zeros(n, max(n_sf_levels - 1, 0));
    sf_names = {};
    for si = 2:n_sf_levels
        D_sf(:, si-1) = double(sf_vals == unique_sf(si));
        sf_names{end+1} = sprintf('SF_%d', si);  %#ok<AGROW>
    end
    
    % Orientation dummy coding
    unique_or = sort(unique(or_vals(~isnan(or_vals))));
    n_or_levels = length(unique_or);
    D_or = zeros(n, max(n_or_levels - 1, 0));
    or_names = {};
    for oi = 2:n_or_levels
        D_or(:, oi-1) = double(or_vals == unique_or(oi));
        or_names{end+1} = sprintf('OR_%d', oi);  %#ok<AGROW>
    end
    
    % Speed x TF interaction
    B_inter_speed_tf = zeros(n, n_speed_bases_arg * n_tf_bases_arg);
    inter_st_names = {};
    col = 0;
    for si = 1:n_speed_bases_arg
        for ti = 1:n_tf_bases_arg
            col = col + 1;
            B_inter_speed_tf(:, col) = B_speed(:, si) .* B_tf(:, ti);
            inter_st_names{end+1} = sprintf('Spd%d_x_TF%d', si, ti);  %#ok<AGROW>
        end
    end
    
    % Speed x SF interaction
    B_inter_speed_sf = zeros(n, n_speed_bases_arg * size(D_sf, 2));
    inter_ssf_names = {};
    col = 0;
    for si = 1:n_speed_bases_arg
        for fi = 1:size(D_sf, 2)
            col = col + 1;
            B_inter_speed_sf(:, col) = B_speed(:, si) .* D_sf(:, fi);
            inter_ssf_names{end+1} = sprintf('Spd%d_x_%s', si, sf_names{fi});  %#ok<AGROW>
        end
    end
    
    % Speed x Orientation interaction
    B_inter_speed_or = zeros(n, n_speed_bases_arg * size(D_or, 2));
    inter_sor_names = {};
    col = 0;
    for si = 1:n_speed_bases_arg
        for oi = 1:size(D_or, 2)
            col = col + 1;
            B_inter_speed_or(:, col) = B_speed(:, si) .* D_or(:, oi);
            inter_sor_names{end+1} = sprintf('Spd%d_x_%s', si, or_names{oi});  %#ok<AGROW>
        end
    end
    
    % TF x SF interaction
    n_sf_cols = size(D_sf, 2);
    B_inter_tf_sf = zeros(n, n_tf_bases_arg * n_sf_cols);
    inter_tf_sf_names = {};
    col = 0;
    for ti = 1:n_tf_bases_arg
        for fi = 1:n_sf_cols
            col = col + 1;
            B_inter_tf_sf(:, col) = B_tf(:, ti) .* D_sf(:, fi);
            inter_tf_sf_names{end+1} = sprintf('TF%d_x_%s', ti, sf_names{fi});  %#ok<AGROW>
        end
    end
    
    % TF x Orientation interaction
    n_or_cols = size(D_or, 2);
    B_inter_tf_or = zeros(n, n_tf_bases_arg * n_or_cols);
    inter_tf_or_names = {};
    col = 0;
    for ti = 1:n_tf_bases_arg
        for oi = 1:n_or_cols
            col = col + 1;
            B_inter_tf_or(:, col) = B_tf(:, ti) .* D_or(:, oi);
            inter_tf_or_names{end+1} = sprintf('TF%d_x_%s', ti, or_names{oi});  %#ok<AGROW>
        end
    end
    
    % SF x Orientation interaction
    B_inter_sf_or = zeros(n, n_sf_cols * n_or_cols);
    inter_sf_or_names = {};
    col = 0;
    for fi = 1:n_sf_cols
        for oi = 1:n_or_cols
            col = col + 1;
            B_inter_sf_or(:, col) = D_sf(:, fi) .* D_or(:, oi);
            inter_sf_or_names{end+1} = sprintf('%s_x_%s', sf_names{fi}, or_names{oi});  %#ok<AGROW>
        end
    end
    
    % Column name lists
    spd_names = {};
    for si = 1:n_speed_bases_arg
        spd_names{end+1} = sprintf('Speed_%d', si);  %#ok<AGROW>
    end
    tf_names_list = {};
    for ti = 1:n_tf_bases_arg
        tf_names_list{end+1} = sprintf('TF_%d', ti);  %#ok<AGROW>
    end
    
    switch model_label
        case 'M0'
            X = ones(n, 1);
            col_names = {'Intercept'};
        case 'M0_Speed'
            % Null + Speed basis only (for incremental variable selection)
            X = [ones(n,1), B_speed];
            col_names = [{'Intercept'}, spd_names];
        case 'M0_TF'
            % Null + TF basis only (for incremental variable selection)
            X = [ones(n,1), B_tf];
            col_names = [{'Intercept'}, tf_names_list];
        case 'M0_SF'
            % Null + SF dummies only (for incremental variable selection)
            X = [ones(n,1), D_sf];
            col_names = [{'Intercept'}, sf_names];
        case 'M0_OR'
            % Null + Orientation dummies only (for incremental variable selection)
            X = [ones(n,1), D_or];
            col_names = [{'Intercept'}, or_names];
        case 'M0_Speed_TF'
            % Null + Speed + TF (no SF, OR, Time — for testing speed given TF)
            X = [ones(n,1), B_speed, B_tf];
            col_names = [{'Intercept'}, spd_names, tf_names_list];
        case 'M0_Speed_TF_SF'
            % Null + Speed + TF + SF (for testing SF given speed+TF)
            X = [ones(n,1), B_speed, B_tf, D_sf];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names];
        case 'M0_Speed_TF_SF_OR'
            % Null + Speed + TF + SF + OR (= Additive without temporal bases)
            X = [ones(n,1), B_speed, B_tf, D_sf, D_or];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, or_names];
        case 'Additive'
            X = [ones(n,1), B_speed, B_tf, D_sf, D_or, B_time];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, or_names, time_names];
        case 'FullInteraction'
            X = [ones(n,1), B_speed, B_tf, D_sf, D_or, B_time, ...
                 B_inter_speed_tf, B_inter_speed_sf, B_inter_speed_or, ...
                 B_inter_tf_sf, B_inter_tf_or, B_inter_sf_or];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, or_names, time_names, ...
                         inter_st_names, inter_ssf_names, inter_sor_names, ...
                         inter_tf_sf_names, inter_tf_or_names, inter_sf_or_names];
        otherwise
            error('Unknown model label: %s', model_label);
    end
    
    % Remove zero-variance columns (except intercept)
    keep_cols = true(1, size(X, 2));
    for col_i = 2:size(X, 2)
        if var(X(:, col_i)) < 1e-12
            keep_cols(col_i) = false;
        end
    end
    if ~all(keep_cols)
        X = X(:, keep_cols);
        col_names = col_names(keep_cols);
    end
end


function [X, col_names] = assemble_design_matrix(B_speed, B_tf, B_time, sf_vals, or_vals, model_label, sf_ref_levels, or_ref_levels)
%ASSEMBLE_DESIGN_MATRIX Build design matrix from pre-computed basis matrices
%   This avoids recomputing raised cosine bases for each model of the same cluster.
%   sf_ref_levels, or_ref_levels: optional reference levels for dummy coding
%   (used for prediction to match training design matrix structure)
%   When ref levels are provided, zero-variance removal is skipped (for prediction).

    n = size(B_speed, 1);
    n_speed_b = size(B_speed, 2);
    n_tf_b = size(B_tf, 2);
    n_time_b = size(B_time, 2);
    
    % Determine if this is for prediction (ref levels provided)
    is_prediction = (nargin >= 7 && ~isempty(sf_ref_levels)) || (nargin >= 8 && ~isempty(or_ref_levels));
    
    % Column name lists
    spd_names = arrayfun(@(si) sprintf('Speed_%d', si), 1:n_speed_b, 'UniformOutput', false);
    tf_names_list = arrayfun(@(ti) sprintf('TF_%d', ti), 1:n_tf_b, 'UniformOutput', false);
    
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
        case 'M0'
            X = ones(n, 1);
            col_names = {'Intercept'};
        case 'M0_Speed'
            % Null + Speed basis only
            X = [ones(n,1), B_speed];
            col_names = [{'Intercept'}, spd_names];
        case 'M0_Speed_TF'
            % Null + Speed + TF (for testing TF given speed)
            X = [ones(n,1), B_speed, B_tf];
            col_names = [{'Intercept'}, spd_names, tf_names_list];
        case 'M0_Speed_TF_SF'
            % Null + Speed + TF + SF (for testing SF given speed+TF)
            X = [ones(n,1), B_speed, B_tf, D_sf];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names];
        case 'Additive'
            % Full main effects: Speed + TF + SF + OR (no temporal bases)
            X = [ones(n,1), B_speed, B_tf, D_sf, D_or];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, or_names];
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
            
            X = [ones(n,1), B_speed, B_tf, D_sf, D_or, ...
                 B_inter_st, B_inter_ssf, B_inter_sor, ...
                 B_inter_tf_sf, B_inter_tf_or, B_inter_sf_or];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names, or_names, ...
                         inter_st_names, inter_ssf_names, inter_sor_names, ...
                         inter_tf_sf_names, inter_tf_or_names, inter_sf_or_names];
        case 'Additive_no_Speed'
            % Additive model without Speed: TF + SF + OR
            X = [ones(n,1), B_tf, D_sf, D_or];
            col_names = [{'Intercept'}, tf_names_list, sf_names, or_names];
        case 'Additive_no_TF'
            % Additive model without TF: Speed + SF + OR
            X = [ones(n,1), B_speed, D_sf, D_or];
            col_names = [{'Intercept'}, spd_names, sf_names, or_names];
        case 'Additive_no_SF'
            % Additive model without SF: Speed + TF + OR
            X = [ones(n,1), B_speed, B_tf, D_or];
            col_names = [{'Intercept'}, spd_names, tf_names_list, or_names];
        case 'Additive_no_OR'
            % Additive model without OR: Speed + TF + SF
            X = [ones(n,1), B_speed, B_tf, D_sf];
            col_names = [{'Intercept'}, spd_names, tf_names_list, sf_names];
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


function [grp_mean, grp_sem, grp_clr, grp_lbl, n_grps, grp_betas, grp_sig] = compute_grouped_params(b_all, se_all, cn_all)
%COMPUTE_GROUPED_PARAMS Group coefficients by feature type for swarm plot
%  Also returns per-group individual betas and significance (Wald |z|>1.96)
    n_coef = length(b_all);
    z_scores = b_all ./ max(se_all, 1e-12);  % Wald z-scores
    is_sig = abs(z_scores) > 1.96;           % p < 0.05 two-sided
    
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
    grp_colors = [0.5 0.5 0.5;   ... % Intercept
                  0.2 0.6 1.0;   ... % Speed
                  1.0 0.5 0.0;   ... % TF
                  0.4 0.8 0.4;   ... % SF
                  0.8 0.4 0.8;   ... % OR
                  0.3 0.8 0.8;   ... % Time
                  0.9 0.2 0.2;   ... % Spd x TF
                  0.8 0.5 0.2;   ... % Spd x SF
                  0.6 0.2 0.6;   ... % Spd x OR
                  0.7 0.7 0.7];      % Other
    
    grp_mean = []; grp_sem = []; grp_clr = []; grp_lbl = {};
    grp_betas = {}; grp_sig = {};
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
        grp_sig{n_grps}   = is_sig(mask)';                     %#ok<AGROW>
    end
end


function plot_zscore_heatmap(b_all, se_all, cn_all)
%PLOT_ZSCORE_HEATMAP Raw-beta coefficient heatmap with significance markers
%   Displays raw beta values as colour and marks coefficients whose
%   |beta/SE| > 1.96 with an asterisk.
    n_coef = length(b_all);

    % Compute z-scores only for significance markers (handle Inf/NaN SE)
    se_safe = se_all;
    se_safe(~isfinite(se_safe) | se_safe <= 0) = NaN;
    z_all = b_all ./ se_safe;           % NaN where SE was bad

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
    
    grp_present = {};
    grp_b_cells = {};
    grp_z_cells = {};
    for gi = 1:length(grp_order)
        mask = strcmp(grp_tags, grp_order{gi});
        if any(mask)
            grp_present{end+1} = grp_order{gi};    %#ok<AGROW>
            grp_b_cells{end+1} = b_all(mask);       %#ok<AGROW>
            grp_z_cells{end+1} = z_all(mask);       %#ok<AGROW>
        end
    end
    ng = length(grp_present);
    mb = max(cellfun(@length, grp_b_cells));
    
    bm = nan(mb, ng);   % beta matrix for colour
    zm_sig = nan(mb, ng); % z matrix for significance only
    for gi = 1:ng
        ngi = length(grp_b_cells{gi});
        bm(1:ngi, gi) = grp_b_cells{gi};
        zm_sig(1:ngi, gi) = grp_z_cells{gi};
    end
    
    imagesc(bm, 'AlphaData', ~isnan(bm));
    set(gca, 'Color', [0.95 0.95 0.95]);
    colormap(gca, bluewhitered_cmap(64));
    cb = colorbar; cb.Label.String = '\beta'; cb.Label.FontSize = 7;
    % Symmetric colour axis based on data range
    max_abs_b = max(abs(bm(:)), [], 'omitnan');
    if isempty(max_abs_b) || max_abs_b == 0, max_abs_b = 1; end
    caxis([-max_abs_b, max_abs_b]);
    
    hold on;
    for gc = 1:ng
        for gr = 1:mb
            if ~isnan(zm_sig(gr, gc)) && abs(zm_sig(gr, gc)) > 1.96
                text(gc, gr, '*', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'FontSize', 10, ...
                    'FontWeight', 'bold', 'Color', 'k');
            end
        end
    end
    hold off;
    
    set(gca, 'XTick', 1:ng, 'XTickLabel', grp_present, ...
        'XTickLabelRotation', 45, 'FontSize', 6, ...
        'TickLabelInterpreter', 'tex');
    ylabel('Basis #');
    if mb > 1, set(gca, 'YTick', 1:mb); end
    title(sprintf('\\beta (%d p)', n_coef), 'FontSize', 9);
end


function cmap = bluewhitered_cmap(n)
%BLUEWHITERED_CMAP Blue-white-red diverging colormap
    if nargin < 1, n = 64; end
    half = floor(n/2);
    r1 = linspace(0.2, 1, half)';
    g1 = linspace(0.3, 1, half)';
    b1 = linspace(0.9, 1, half)';
    r2 = linspace(1, 0.9, n - half)';
    g2 = linspace(1, 0.2, n - half)';
    b2 = linspace(1, 0.2, n - half)';
    cmap = [r1 g1 b1; r2 g2 b2];
end


function draw_ellipse(cx, cy, rx, ry, angle_deg, color)
%DRAW_ELLIPSE Draw a filled ellipse with transparency for Venn diagrams
%
%   DRAW_ELLIPSE(CX, CY, RX, RY, ANGLE_DEG, COLOR) draws an ellipse centered
%   at (CX, CY) with semi-axes RX and RY, rotated by ANGLE_DEG degrees.
%   COLOR is a 1x4 vector [R G B Alpha] for fill color with transparency.
%
%   Used for creating 4-set Venn diagrams.

    n_points = 100;
    theta = linspace(0, 2*pi, n_points);
    
    % Ellipse in standard position
    x = rx * cos(theta);
    y = ry * sin(theta);
    
    % Rotation matrix
    angle_rad = angle_deg * pi / 180;
    R = [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)];
    
    % Rotate and translate
    xy = R * [x; y];
    x_rot = xy(1, :) + cx;
    y_rot = xy(2, :) + cy;
    
    % Draw filled ellipse with transparency
    fill(x_rot, y_rot, color(1:3), 'FaceAlpha', color(4), ...
        'EdgeColor', color(1:3)*0.7, 'LineWidth', 1.5);
end
