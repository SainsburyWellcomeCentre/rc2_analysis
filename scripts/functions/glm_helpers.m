classdef glm_helpers
%GLM_HELPERS Helper functions for GLM single cluster analysis
%   Static methods for basis functions, design matrix construction,
%   model fitting, cross-validation, and visualization utilities.
%
%   Usage:
%       B = glm_helpers.make_raised_cosine_basis(x, n_bases, x_min, x_max);
%       result = glm_helpers.fit_poisson_glm(X, y, offset, lambda_ridge);
%
%   Categories:
%       DRY Helpers:            compute_raised_cosine_bumps, make_basis_names,
%                               make_dummy_columns, make_pairwise_interaction,
%                               make_categorical_dummies, make_basis_categorical_interaction,
%                               remove_zero_variance_cols, run_selection_phase
%       Basis Functions:        make_raised_cosine_basis, make_onset_kernel_basis
%       Design Matrix:          assemble_design_matrix, assemble_design_matrix_selected,
%                               build_design_matrix_from_colnames, predict_design_matrix
%       Model Fitting:          fit_poisson_glm, cross_validate_glm, forward_select_model
%       Visualization Helpers:  compute_grouped_params, extract_basis_label, format_radians
%       Parsing/Utilities:      ternary, parse_cloud_name, get_component_column

    methods (Static)
        
        function out = ternary(cond_val, val_true, val_false)
        %TERNARY Inline conditional helper
        %   out = glm_helpers.ternary(condition, value_if_true, value_if_false)
            if cond_val
                out = val_true;
            else
                out = val_false;
            end
        end
        
        
        function B = compute_raised_cosine_bumps(x, centers, width)
        %COMPUTE_RAISED_COSINE_BUMPS Core raised cosine bump computation (DRY helper)
        %   B = glm_helpers.compute_raised_cosine_bumps(x, centers, width)
        %
        %   Shared by make_raised_cosine_basis and make_onset_kernel_basis.
        %
        %   Inputs:
        %       x       - Input values [n x 1]
        %       centers - Center positions for each basis [1 x n_bases]
        %       width   - Width parameter for cosine bumps
        %
        %   Output:
        %       B - Basis matrix [n x n_bases]
            if isempty(x) || isempty(centers)
                B = zeros(length(x), length(centers));
                return;
            end

            x = double(x(:));
            centers = double(centers(:)');
            width = double(width);

            if ~isfinite(width) || width == 0
                width = 1;
            end

            n = length(x);
            n_bases = length(centers);
            B = zeros(n, n_bases);
            for bi = 1:n_bases
                z = (x - centers(bi)) / width * pi;
                z = max(min(z, pi), -pi);
                B(:, bi) = 0.5 * (1 + builtin('cos', z));
            end
        end
        
        
        function names = make_basis_names(prefix, n_bases)
        %MAKE_BASIS_NAMES Generate column names for basis functions (DRY helper)
        %   names = glm_helpers.make_basis_names(prefix, n_bases)
        %
        %   Examples:
        %       make_basis_names('Speed', 5) -> {'Speed_1', 'Speed_2', ..., 'Speed_5'}
            names = arrayfun(@(i) sprintf('%s_%d', prefix, i), 1:n_bases, 'UniformOutput', false);
        end
        
        
        function [D, names, levels] = make_dummy_columns(vals, ref_levels, name_prefix, name_format)
        %MAKE_DUMMY_COLUMNS Create dummy-coded columns for categorical variable (DRY helper)
        %   [D, names, levels] = glm_helpers.make_dummy_columns(vals, ref_levels, prefix, format)
        %
        %   Creates dummy variables with first level as reference (excluded).
        %
        %   Inputs:
        %       vals        - Categorical values [n x 1]
        %       ref_levels  - Reference levels (empty = auto-detect from vals)
        %       name_prefix - Prefix for column names (e.g., 'SF', 'OR')
        %       name_format - Format string for level index (e.g., '%d')
        %                     Uses INDEX starting from 2 (since level 1 is reference)
        %
        %   Outputs:
        %       D      - Dummy matrix [n x (n_levels-1)]
        %       names  - Column names (e.g., 'SF_2', 'SF_3' for indices)
        %       levels - Non-reference levels used (actual values for matching)
            n = length(vals);
            
            % Determine unique levels
            if ~isempty(ref_levels)
                unique_levels = ref_levels(:)';
            else
                unique_levels = sort(unique(vals(~isnan(vals))));
            end
            
            % Exclude reference (first) level
            levels = unique_levels(2:end);
            n_dummies = length(levels);
            
            D = zeros(n, max(n_dummies, 0));
            names = cell(1, n_dummies);
            for di = 1:n_dummies
                D(:, di) = double(vals == levels(di));
                % Use index (di+1) because level 1 is reference, so dummies are levels 2, 3, ...
                names{di} = sprintf([name_prefix '_' name_format], di + 1);
            end
        end
        
        
        function [B_inter, inter_names] = make_pairwise_interaction(B1, B2, names1, names2, name_fmt)
        %MAKE_PAIRWISE_INTERACTION Create all pairwise interaction columns (DRY helper)
        %   [B, names] = glm_helpers.make_pairwise_interaction(B1, B2, names1, names2, fmt)
        %
        %   Creates interaction columns as element-wise products of all pairs.
        %
        %   Inputs:
        %       B1, B2       - Basis/dummy matrices [n x n1], [n x n2]
        %       names1, names2 - Column name components
        %       name_fmt     - Format string (e.g., '%s_x_%s')
        %
        %   Outputs:
        %       B_inter     - Interaction matrix [n x (n1*n2)]
        %       inter_names - Column names
            n = size(B1, 1);
            n1 = size(B1, 2);
            n2 = size(B2, 2);
            
            B_inter = zeros(n, n1 * n2);
            inter_names = cell(1, n1 * n2);
            col = 0;
            for i1 = 1:n1
                for i2 = 1:n2
                    col = col + 1;
                    B_inter(:, col) = B1(:, i1) .* B2(:, i2);
                    inter_names{col} = sprintf(name_fmt, names1{i1}, names2{i2});
                end
            end
        end
        
        
        function [D, names] = make_categorical_dummies(vals, levels, name_prefix, name_format)
        %MAKE_CATEGORICAL_DUMMIES Create dummy matrix from categorical values (DRY helper)
        %   [D, names] = glm_helpers.make_categorical_dummies(vals, levels, prefix, format)
        %
        %   Creates dummy columns for each level, with NaN handling.
        %   Unlike make_dummy_columns, this does NOT exclude a reference level.
        %
        %   Inputs:
        %       vals        - Categorical values [n x 1]
        %       levels      - Levels to create dummies for
        %       name_prefix - Prefix for column names (e.g., 'SF', 'OR')
        %       name_format - Printf format for level value (e.g., '%.4f', '%.3f')
        %
        %   Outputs:
        %       D     - Dummy matrix [n x n_levels]
        %       names - Cell array of column names
            n = length(vals);
            n_levels = length(levels);
            D = zeros(n, n_levels);
            names = cell(1, n_levels);
            for li = 1:n_levels
                D(:, li) = double(vals == levels(li));
                D(isnan(vals), li) = 0;
                names{li} = sprintf([name_prefix '_' name_format], levels(li));
            end
        end
        
        
        function [B_inter, names] = make_basis_categorical_interaction(B, D, basis_prefix, cat_names, name_fmt)
        %MAKE_BASIS_CATEGORICAL_INTERACTION Create basis x categorical dummy interactions (DRY helper)
        %   [B_inter, names] = glm_helpers.make_basis_categorical_interaction(B, D, prefix, cat_names, fmt)
        %
        %   Creates interaction columns: B(:,i) .* D(:,j) for all i,j pairs.
        %
        %   Inputs:
        %       B           - Basis matrix [n x n_bases]
        %       D           - Dummy matrix [n x n_dummies]
        %       basis_prefix - Prefix for basis (e.g., 'Spd', 'TF')
        %       cat_names   - Column names for dummies (e.g., {'SF_0.0060', 'SF_0.0120'})
        %       name_fmt    - Format string (e.g., '%s%d_x_%s')
        %
        %   Outputs:
        %       B_inter - Interaction matrix [n x (n_bases * n_dummies)]
        %       names   - Cell array of column names
            n = size(B, 1);
            n_bases = size(B, 2);
            n_dummies = size(D, 2);
            
            B_inter = zeros(n, n_bases * n_dummies);
            names = cell(1, n_bases * n_dummies);
            col = 0;
            for bi = 1:n_bases
                for di = 1:n_dummies
                    col = col + 1;
                    B_inter(:, col) = B(:, bi) .* D(:, di);
                    names{col} = sprintf(name_fmt, basis_prefix, bi, cat_names{di});
                end
            end
        end
        
        
        function [sf_val, vx_val, or_val] = parse_cloud_name(cname, sf_keys_arg, sf_values_map_arg, or_keys_arg, or_values_map_arg)
        %PARSE_CLOUD_NAME Extract SF, VX, orientation from a cloud name string
        %   [sf, vx, or] = glm_helpers.parse_cloud_name(name, sf_keys, sf_map, or_keys, or_map)
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


        function [mc_sequence, cloud_names] = load_motion_cloud_metadata(mc_sequence_path, mc_folders_path)
        %LOAD_MOTION_CLOUD_METADATA Load presentation sequence and cloud names
        %   [mc_sequence, cloud_names] = glm_helpers.load_motion_cloud_metadata(seq_path, folders_path)
            mc_sequence = [];
            cloud_names = {};

            P = load(mc_sequence_path);
            if isfield(P, 'presentation_sequence')
                mc_sequence = P.presentation_sequence;
            else
                fns_seq = fieldnames(P);
                for i = 1:numel(fns_seq)
                    v = P.(fns_seq{i});
                    if isnumeric(v) && isvector(v)
                        mc_sequence = v;
                        break;
                    end
                end
            end

            S = load(mc_folders_path);
            fns_fold = fieldnames(S);
            for i = 1:numel(fns_fold)
                v = S.(fns_fold{i});
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
        end


        function trial_stim = build_trial_stim_lookup(mc_sequence, cloud_names, exclude_patterns, ...
                sf_keys_arg, sf_values_map_arg, or_keys_arg, or_values_map_arg, ...
                batch_patterns, batch_gains)
        %BUILD_TRIAL_STIM_LOOKUP Build trial-level stimulus lookup struct
        %   trial_stim = glm_helpers.build_trial_stim_lookup(...)
        %
        %   Returns struct array with fields:
        %       sf, batch_gain, or, excluded
            contains_any = @(str, subs) any(cellfun(@(s) contains(str, s), subs));

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
                    [sf_v, ~, or_v] = glm_helpers.parse_cloud_name(cname, sf_keys_arg, sf_values_map_arg, or_keys_arg, or_values_map_arg);
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
        end


        function [ok, template] = get_trial_speed_template(aligned_obj, n_samples)
        %GET_TRIAL_SPEED_TEMPLATE Returns a fixed-length speed profile template per trial
        %   [ok, template] = glm_helpers.get_trial_speed_template(aligned_obj, n_samples)
            ok = false;
            template = nan(1, n_samples);

            tr_probe_t = aligned_obj.probe_t;
            tr_vel = aligned_obj.velocity;
            tr_mmask = aligned_obj.motion_mask;

            motion_idx = find(tr_mmask);
            if isempty(motion_idx)
                return;
            end
            i0 = motion_idx(1);
            i1 = motion_idx(end);
            if i1 <= i0
                return;
            end

            t_seg = tr_probe_t(i0:i1);
            v_seg = abs(tr_vel(i0:i1));

            if length(t_seg) < 5 || range(t_seg) <= 0
                return;
            end

            if exist('smoothdata', 'file')
                v_seg = smoothdata(v_seg, 'movmedian', 5);
            end

            t_norm = (t_seg - t_seg(1)) ./ (t_seg(end) - t_seg(1));
            t_target = linspace(0, 1, n_samples);
            template = interp1(t_norm, v_seg, t_target, 'linear', 'extrap');

            if any(~isfinite(template))
                return;
            end

            ok = true;
        end
        
        
        function B = make_raised_cosine_basis(x, n_bases, x_min, x_max)
        %MAKE_RAISED_COSINE_BASIS Raised cosine bases on log-shifted axis (Weber-law scaling)
        %   B = glm_helpers.make_raised_cosine_basis(x, n_bases, x_min, x_max)
        %
        %   Inputs:
        %       x        - Input values (e.g., speed or temporal frequency)
        %       n_bases  - Number of basis functions
        %       x_min    - Minimum value for basis placement
        %       x_max    - Maximum value for basis placement
        %
        %   Output:
        %       B - Basis matrix [length(x) x n_bases]
            epsilon = 1e-6;

            if isempty(x)
                B = zeros(0, n_bases);
                return;
            end

            x = double(x(:));
            x_min = double(x_min);
            x_max = double(x_max);

            support_min = min([x; x_min; x_max], [], 'omitnan');
            if ~isfinite(support_min)
                support_min = 0;
            end
            shift = max(0, -support_min + epsilon);

            x_shifted = x + shift;
            x_shifted(~isfinite(x_shifted) | x_shifted <= 0) = epsilon;

            min_shifted = x_min + shift;
            max_shifted = x_max + shift;
            if ~isfinite(min_shifted) || min_shifted <= 0
                min_shifted = epsilon;
            end
            if ~isfinite(max_shifted) || max_shifted <= 0
                max_shifted = min_shifted + 1;
            end
            if max_shifted <= min_shifted
                max_shifted = min_shifted + 1;
            end

            log_x = builtin('log', x_shifted);
            log_min = builtin('log', min_shifted);
            log_max = builtin('log', max_shifted);
            
            centers = linspace(log_min, log_max, n_bases);
            width = glm_helpers.ternary(n_bases > 1, ...
                (log_max - log_min) / (n_bases - 1) * 1.5, ...
                (log_max - log_min) * 1.5);
            
            B = glm_helpers.compute_raised_cosine_bumps(log_x, centers, width);
        end
        
        
        function B = make_onset_kernel_basis(t_since_onset, n_bases, t_max)
        %MAKE_ONSET_KERNEL_BASIS Causal raised cosine bases for onset dynamics (Park et al. 2014 style)
        %   B = glm_helpers.make_onset_kernel_basis(t_since_onset, n_bases, t_max)
        %
        %   Creates raised cosine bases that are NON-ZERO only for t >= 0 (causal).
        %   For t < 0 (stationary periods before motion onset), all bases are zero.
        %   This captures the transient response at motion onset and adaptation.
        %
        %   Inputs:
        %       t_since_onset - time since motion onset (negative for stationary, positive for motion)
        %       n_bases       - number of raised cosine bases
        %       t_max         - maximum time covered by the kernel (e.g., 2.0 s)
        %
        %   Output:
        %       B - Basis matrix [length(t_since_onset) x n_bases]
        %
        %   Following Park et al. 2014: bases are raised cosine bumps on linear axis,
        %   spanning from 0 to t_max, with overlapping coverage.

            n = length(t_since_onset);
            B = zeros(n, n_bases);
            
            centers = linspace(0, t_max, n_bases);
            width = glm_helpers.ternary(n_bases > 1, t_max / (n_bases - 1), t_max) * 1.5;
            
            % Mask: only compute for t >= 0 (motion period)
            motion_mask = t_since_onset >= 0;
            t_motion = t_since_onset(motion_mask);
            
            if ~isempty(t_motion)
                B(motion_mask, :) = glm_helpers.compute_raised_cosine_bumps(t_motion, centers, width);
            end
            % B remains zero for t < 0 (stationary periods)
        end
        
        
        function [X, col_names] = assemble_design_matrix(B_speed, B_tf, B_onset, sf_vals, or_vals, model_label, sf_ref_levels, or_ref_levels)
        %ASSEMBLE_DESIGN_MATRIX Build design matrix from pre-computed basis matrices
        %   [X, col_names] = glm_helpers.assemble_design_matrix(B_speed, B_tf, B_onset, ...
        %                       sf_vals, or_vals, model_label, sf_ref_levels, or_ref_levels)
        %
        %   This avoids recomputing raised cosine bases for each model of the same cluster.
        %
        %   Inputs:
        %       B_speed, B_tf, B_onset - Pre-computed basis matrices
        %       sf_vals, or_vals       - Spatial frequency and orientation values
        %       model_label            - One of: 'Null', 'M0', 'M0_Speed', 'M0_Speed_TF',
        %                                'M0_Speed_TF_SF', 'Additive', 'FullInteraction',
        %                                'Additive_no_Speed', 'Additive_no_TF',
        %                                'Additive_no_SF', 'Additive_no_OR'
        %       sf_ref_levels, or_ref_levels - Optional reference levels for dummy coding
        %                                      (used for prediction to match training structure)
        %
        %   Outputs:
        %       X         - Design matrix [n_obs x n_params]
        %       col_names - Cell array of column names

            if nargin < 7, sf_ref_levels = []; end
            if nargin < 8, or_ref_levels = []; end

            n = size(B_speed, 1);
            n_speed_b = size(B_speed, 2);
            n_tf_b = size(B_tf, 2);
            n_onset_b = size(B_onset, 2);
            
            % Determine if this is for prediction (ref levels provided)
            is_prediction = ~isempty(sf_ref_levels) || ~isempty(or_ref_levels);
            
            % Column name lists using DRY helper
            spd_names = glm_helpers.make_basis_names('Speed', n_speed_b);
            tf_names_list = glm_helpers.make_basis_names('TF', n_tf_b);
            onset_names = glm_helpers.make_basis_names('Onset', n_onset_b);
            
            % SF dummy coding using DRY helper
            [D_sf, sf_names] = glm_helpers.make_dummy_columns(sf_vals, sf_ref_levels, 'SF', '%d');
            
            % Orientation dummy coding using DRY helper
            [D_or, or_names] = glm_helpers.make_dummy_columns(or_vals, or_ref_levels, 'OR', '%d');
            
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
                    % Create short names for interaction naming
                    spd_short = arrayfun(@(i) sprintf('Spd%d', i), 1:n_speed_b, 'UniformOutput', false);
                    tf_short = arrayfun(@(i) sprintf('TF%d', i), 1:n_tf_b, 'UniformOutput', false);
                    
                    % Build all pairwise interactions using DRY helper
                    [B_inter_st, inter_st_names] = glm_helpers.make_pairwise_interaction(...
                        B_speed, B_tf, spd_short, tf_short, '%s_x_%s');
                    [B_inter_ssf, inter_ssf_names] = glm_helpers.make_pairwise_interaction(...
                        B_speed, D_sf, spd_short, sf_names, '%s_x_%s');
                    [B_inter_sor, inter_sor_names] = glm_helpers.make_pairwise_interaction(...
                        B_speed, D_or, spd_short, or_names, '%s_x_%s');
                    [B_inter_tf_sf, inter_tf_sf_names] = glm_helpers.make_pairwise_interaction(...
                        B_tf, D_sf, tf_short, sf_names, '%s_x_%s');
                    [B_inter_tf_or, inter_tf_or_names] = glm_helpers.make_pairwise_interaction(...
                        B_tf, D_or, tf_short, or_names, '%s_x_%s');
                    [B_inter_sf_or, inter_sf_or_names] = glm_helpers.make_pairwise_interaction(...
                        D_sf, D_or, sf_names, or_names, '%s_x_%s');
                    
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
                [X, col_names] = glm_helpers.remove_zero_variance_cols(X, col_names, 1e-12);
            end
        end
        
        
        function result = fit_poisson_glm(X, y, offset, lambda_ridge)
        %FIT_POISSON_GLM Fit Poisson GLM with log link via IRLS
        %   result = glm_helpers.fit_poisson_glm(X, y, offset, lambda_ridge)
        %
        %   Inputs:
        %       X            - Design matrix [n x p]
        %       y            - Response (spike counts) [n x 1]
        %       offset       - Offset term (log bin width) [n x 1]
        %       lambda_ridge - Ridge regularization parameter (default: 0)
        %
        %   Output:
        %       result - Struct with fields:
        %           beta, se, log_likelihood, aic, bic, deviance, n_params, n_obs,
        %           predicted_count, pearson_residuals, dispersion, converged, lambda_ridge
            if nargin < 4, lambda_ridge = 0; end
            
            % Ensure y is double (histcounts returns uint32)
            y = double(y);
            
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
        %   [cv_ll, cv_bps, fold_ll, cv_pred] = glm_helpers.cross_validate_glm(X, y, offset, fold_ids, lambda)
        %
        %   cv_predicted_count returns mu = E[Y] (expected spike count), not firing rate.
        %   To get firing rate lambda(t), divide by exposure time: lambda = mu / t
        %
        %   Inputs:
        %       X            - Design matrix
        %       y            - Spike counts
        %       offset       - Offset term (log bin width)
        %       fold_ids     - CV fold assignments
        %       lambda_ridge - Ridge parameter (default: 0)
        %
        %   Outputs:
        %       cv_ll             - Total cross-validated log-likelihood
        %       cv_bits_per_spike - Bits per spike (information metric)
        %       fold_ll           - Log-likelihood per fold
        %       cv_predicted_count - Cross-validated predictions
            if nargin < 5, lambda_ridge = 0; end
            
            % Ensure y is double (histcounts returns uint32)
            y = double(y);
            
            folds = unique(fold_ids);
            K = length(folds);
            
            total_ll = 0;
            total_spikes = 0;
            fold_ll = zeros(K, 1);
            cv_predicted_count = nan(size(y));
            
            for fi = 1:K
                test_idx = fold_ids == folds(fi);
                train_idx = ~test_idx;
                
                res = glm_helpers.fit_poisson_glm(X(train_idx,:), y(train_idx), offset(train_idx), lambda_ridge);
                
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
        %   [vars, history, final_bps, null_bps] = glm_helpers.forward_select_model(...)
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
        %       threshold              - Minimum delta bps for inclusion
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
            [X_null, ~] = glm_helpers.assemble_design_matrix_selected(B_speed, B_tf, B_onset, ...
                sf_v, or_v, {}, sf_ref_levels, or_ref_levels);
            [~, null_cv_bps] = glm_helpers.cross_validate_glm(X_null, y, offset, fold_ids, 0);
            
            selected_vars = {};
            current_cv_bps = null_cv_bps;
            selection_history = struct('round', {}, 'phase', {}, 'tested', {}, 'delta_bps', {}, ...
                'best_candidate', {}, 'added', {}, 'cv_bps_after', {});
            
            round_num = 0;
            
            % =====================================================================
            % PHASE 1: Main effects only
            % =====================================================================
            fprintf('    --- Phase 1: Main Effects ---\n');
            [selected_vars, ~, round_num, current_cv_bps, phase1_history] = ...
                glm_helpers.run_selection_phase(1, 'main', main_effects, selected_vars, ...
                    current_cv_bps, round_num, B_speed, B_tf, B_onset, ...
                    sf_v, or_v, y, offset, fold_ids, threshold, sf_ref_levels, or_ref_levels);
            for hi = 1:length(phase1_history)
                selection_history(end+1) = phase1_history(hi); %#ok<AGROW>
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
                if all(ismember(parents, selected_vars))
                    eligible_interactions{end+1} = int_name; %#ok<AGROW>
                end
            end
            
            if isempty(eligible_interactions)
                fprintf('    No eligible interactions (need both parent main effects selected)\n');
            else
                fprintf('    Eligible interactions: %s\n', strjoin(eligible_interactions, ', '));
                [selected_vars, ~, round_num, current_cv_bps, phase2_history] = ...
                    glm_helpers.run_selection_phase(2, 'interaction', eligible_interactions, selected_vars, ...
                        current_cv_bps, round_num, B_speed, B_tf, B_onset, ...
                        sf_v, or_v, y, offset, fold_ids, threshold, sf_ref_levels, or_ref_levels);
                for hi = 1:length(phase2_history)
                    selection_history(end+1) = phase2_history(hi); %#ok<AGROW>
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
        %   [X, col_names] = glm_helpers.assemble_design_matrix_selected(...)
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

            if nargin < 7, sf_ref_levels = []; end
            if nargin < 8, or_ref_levels = []; end

            n = size(B_speed, 1);
            n_speed_bases = size(B_speed, 2);
            n_tf_bases = size(B_tf, 2);
            n_onset_bases = size(B_onset, 2);
            
            % Determine SF levels for dummy coding
            % When sf_ref_levels is provided (prediction mode), use it as the complete set of levels
            % and exclude the first level (reference). This matches assemble_design_matrix behavior.
            if ~isempty(sf_ref_levels)
                unique_sf = sf_ref_levels(:)';
                sf_levels = unique_sf(2:end);  % Exclude first (reference) level
            else
                sf_unique = unique(sf_vals(~isnan(sf_vals)));
                sf_levels = sf_unique(2:end);  % First level is reference
            end
            
            % Determine OR levels for dummy coding (same logic as SF)
            if ~isempty(or_ref_levels)
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
            has_speed = ismember('Speed', selected_vars);
            has_tf = ismember('TF', selected_vars);
            has_sf = ismember('SF', selected_vars);
            has_or = ismember('OR', selected_vars);
            
            has_speed_x_tf = ismember('Speed_x_TF', selected_vars);
            has_speed_x_sf = ismember('Speed_x_SF', selected_vars);
            has_speed_x_or = ismember('Speed_x_OR', selected_vars);
            has_tf_x_sf = ismember('TF_x_SF', selected_vars);
            has_tf_x_or = ismember('TF_x_OR', selected_vars);
            has_sf_x_or = ismember('SF_x_OR', selected_vars);
            
            % Pre-compute dummy matrices for SF and OR (used by main effects and interactions)
            [D_sf, sf_names] = glm_helpers.make_categorical_dummies(sf_vals, sf_levels, 'SF', '%.4f');
            [D_or, or_names] = glm_helpers.make_categorical_dummies(or_vals, or_levels, 'OR', '%.3f');
            
            % Short names for interactions (preserve original naming convention without prefix underscore)
            sf_int_names = arrayfun(@(v) sprintf('SF%.4f', v), sf_levels, 'UniformOutput', false);
            or_int_names = arrayfun(@(v) sprintf('OR%.3f', v), or_levels, 'UniformOutput', false);
            
            % --- Add main effects ---
            if has_speed
                X = [X, B_speed];
                col_names = [col_names, glm_helpers.make_basis_names('Speed', n_speed_bases)];
            end
            
            if has_tf
                X = [X, B_tf];
                col_names = [col_names, glm_helpers.make_basis_names('TF', n_tf_bases)];
            end
            
            if has_sf
                X = [X, D_sf];
                col_names = [col_names, sf_names];
            end
            
            if has_or
                X = [X, D_or];
                col_names = [col_names, or_names];
            end
            
            % --- Add interactions ---
            % For interactions, we create the interaction columns directly
            % (products of bases/dummies), without requiring main effects to be selected
            
            if has_speed_x_tf
                spd_names = arrayfun(@(i) sprintf('Spd%d', i), 1:n_speed_bases, 'UniformOutput', false);
                tf_names_short = arrayfun(@(i) sprintf('TF%d', i), 1:n_tf_bases, 'UniformOutput', false);
                [B_int, int_names] = glm_helpers.make_pairwise_interaction(B_speed, B_tf, spd_names, tf_names_short, '%s_x_%s');
                X = [X, B_int];
                col_names = [col_names, int_names];
            end
            
            if has_speed_x_sf
                [B_int, int_names] = glm_helpers.make_basis_categorical_interaction(B_speed, D_sf, 'Spd', sf_int_names, '%s%d_x_%s');
                X = [X, B_int];
                col_names = [col_names, int_names];
            end
            
            if has_speed_x_or
                [B_int, int_names] = glm_helpers.make_basis_categorical_interaction(B_speed, D_or, 'Spd', or_int_names, '%s%d_x_%s');
                X = [X, B_int];
                col_names = [col_names, int_names];
            end
            
            if has_tf_x_sf
                [B_int, int_names] = glm_helpers.make_basis_categorical_interaction(B_tf, D_sf, 'TF', sf_int_names, '%s%d_x_%s');
                X = [X, B_int];
                col_names = [col_names, int_names];
            end
            
            if has_tf_x_or
                [B_int, int_names] = glm_helpers.make_basis_categorical_interaction(B_tf, D_or, 'TF', or_int_names, '%s%d_x_%s');
                X = [X, B_int];
                col_names = [col_names, int_names];
            end
            
            if has_sf_x_or
                [B_int, int_names] = glm_helpers.make_pairwise_interaction(D_sf, D_or, sf_int_names, or_int_names, '%s_x_%s');
                X = [X, B_int];
                col_names = [col_names, int_names];
            end
            
            % Remove zero-variance columns
            [X, col_names] = glm_helpers.remove_zero_variance_cols(X, col_names, 1e-10);
        end
        
        
        function [grp_mean, grp_sem, grp_clr, grp_lbl, n_grps, grp_betas, grp_cn] = compute_grouped_params(b_all, se_all, cn_all)
        %COMPUTE_GROUPED_PARAMS Group coefficients by feature type for swarm plot
        %   [mean, sem, clr, lbl, n, betas, cn] = glm_helpers.compute_grouped_params(b, se, cn)
        %
        %   Returns per-group individual betas and column names for basis labeling
        %
        %   Inputs:
        %       b_all  - All beta coefficients
        %       se_all - Standard errors (unused but kept for API)
        %       cn_all - Column names
        %
        %   Outputs:
        %       grp_mean  - Mean beta per group
        %       grp_sem   - SEM of betas per group
        %       grp_clr   - Color per group [n_groups x 3]
        %       grp_lbl   - Label per group
        %       n_grps    - Number of groups
        %       grp_betas - Cell array of individual betas per group
        %       grp_cn    - Cell array of column names per group
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
            % Colors matched to model barplot color scheme
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
        %   lbl = glm_helpers.extract_basis_label(col_name)
        %
        %   Examples:
        %       'Speed_3' -> '3'
        %       'TF_2' -> '2'
        %       'Spd2_x_TF4' -> '2-4'
        %       'Spd1_x_SF_2' -> '1-2'
        %       'Intercept' -> ''
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
        %   str = glm_helpers.format_radians(val)
        %
        %   Converts common angles to nice pi notation (e.g., -pi/4, 0, pi/4, pi/2)
            
            tol = 1e-6;
            if abs(val) < tol
                str = '0';
            elseif abs(val - pi/2) < tol
                str = 'pi/2';
            elseif abs(val + pi/2) < tol
                str = '-pi/2';
            elseif abs(val - pi/4) < tol
                str = 'pi/4';
            elseif abs(val + pi/4) < tol
                str = '-pi/4';
            elseif abs(val - pi) < tol
                str = 'pi';
            elseif abs(val + pi) < tol
                str = '-pi';
            elseif abs(val - 3*pi/4) < tol
                str = '3pi/4';
            elseif abs(val + 3*pi/4) < tol
                str = '-3pi/4';
            else
                % General case: show decimal radians
                str = sprintf('%.2f', val);
            end
        end
        
        
        function X = build_design_matrix_from_colnames(col_names, B_speed, B_tf, B_onset, sf_val, or_val, sf_all_levels, or_all_levels)
        %BUILD_DESIGN_MATRIX_FROM_COLNAMES Build design matrix matching stored column names
        %   X = glm_helpers.build_design_matrix_from_colnames(col_names, B_speed, B_tf, ...
        %           B_onset, sf_val, or_val, sf_all_levels, or_all_levels)
        %
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
                        comp1 = glm_helpers.get_component_column(parts{1}, B_speed, B_tf, sf_val, or_val, sf_dummy_levels, or_dummy_levels);
                        % Get second component
                        comp2 = glm_helpers.get_component_column(parts{2}, B_speed, B_tf, sf_val, or_val, sf_dummy_levels, or_dummy_levels);
                        X(:, ci) = comp1 .* comp2;
                    end
                end
            end
        end
        
        
        function col = get_component_column(name, B_speed, B_tf, sf_val, or_val, sf_levels, or_levels)
        %GET_COMPONENT_COLUMN Get a single column for interaction term component
        %   col = glm_helpers.get_component_column(name, B_speed, B_tf, sf_val, or_val, sf_levels, or_levels)
        %
        %   Helper function for build_design_matrix_from_colnames to parse
        %   individual components of interaction terms.
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
        
        
        function [X, col_names] = remove_zero_variance_cols(X, col_names, tol)
        %REMOVE_ZERO_VARIANCE_COLS Remove columns with near-zero variance
        %   [X, col_names] = glm_helpers.remove_zero_variance_cols(X, col_names, tol)
        %
        %   Always keeps the first column (intercept). Removes columns whose
        %   variance is below the tolerance threshold.
        %
        %   Inputs:
        %       X         - Design matrix [n x p]
        %       col_names - Cell array of column names
        %       tol       - Variance threshold (default: 1e-10)
        %
        %   Outputs:
        %       X         - Pruned design matrix
        %       col_names - Pruned column names
            if nargin < 3, tol = 1e-10; end
            col_vars = var(X, 0, 1);
            keep_cols = col_vars > tol | (1:size(X,2)) == 1;  % Always keep intercept
            if ~all(keep_cols)
                X = X(:, keep_cols);
                col_names = col_names(keep_cols);
            end
        end
        
        
        function X = predict_design_matrix(model_label, col_names_train, B_speed, B_tf, ...
                B_onset, sf_val, or_val, sf_levels, or_levels)
        %PREDICT_DESIGN_MATRIX Build prediction design matrix for any model type
        %   X = glm_helpers.predict_design_matrix(model_label, col_names_train, ...)
        %
        %   Dispatches between build_design_matrix_from_colnames (for 'Selected' model)
        %   and assemble_design_matrix (for named models like 'Null', 'Additive', etc.).
        %
        %   This eliminates the repeated if/else pattern in prediction code.
        %
        %   Inputs:
        %       model_label      - Model name ('Selected', 'Null', 'Additive', etc.)
        %       col_names_train  - Column names from training (used only for 'Selected')
        %       B_speed, B_tf, B_onset - Basis matrices for prediction data
        %       sf_val, or_val   - SF and OR values for prediction data
        %       sf_levels        - All SF levels used during training
        %       or_levels        - All OR levels used during training
        %
        %   Output:
        %       X - Design matrix matching the trained model structure
            if strcmp(model_label, 'Selected')
                X = glm_helpers.build_design_matrix_from_colnames(col_names_train, ...
                    B_speed, B_tf, B_onset, sf_val, or_val, sf_levels, or_levels);
            else
                [X, ~] = glm_helpers.assemble_design_matrix(B_speed, B_tf, ...
                    B_onset, sf_val, or_val, model_label, sf_levels, or_levels);
            end
        end
        
        
        function [selected_vars, remaining, round_num, current_cv_bps, phase_history] = ...
                run_selection_phase(phase_num, phase_label, remaining, selected_vars, ...
                    current_cv_bps, round_num, B_speed, B_tf, B_onset, ...
                    sf_v, or_v, y, offset, fold_ids, threshold, sf_ref_levels, or_ref_levels)
        %RUN_SELECTION_PHASE Core greedy-forward loop for one phase of selection
        %   Shared by Phase 1 (main effects) and Phase 2 (interactions) in
        %   forward_select_model, eliminating code duplication.
        %
        %   Iteratively tests each remaining candidate, adds the best if its
        %   delta CV bits/spike exceeds the threshold, and repeats until no
        %   candidate passes.
        %
        %   Inputs:
        %       phase_num, phase_label - Phase identifier (1/'main' or 2/'interaction')
        %       remaining              - Cell array of candidate variable names
        %       selected_vars          - Currently selected variables
        %       current_cv_bps         - Current CV bits/spike
        %       round_num              - Current round counter
        %       B_speed, B_tf, B_onset - Pre-computed basis matrices
        %       sf_v, or_v             - Categorical variable vectors
        %       y                      - Spike counts
        %       offset                 - Log of bin width
        %       fold_ids               - CV fold assignments
        %       threshold              - Minimum delta bps for inclusion
        %       sf_ref_levels, or_ref_levels - Reference levels for dummy coding
        %
        %   Outputs:
        %       selected_vars  - Updated selected variables
        %       remaining      - Updated remaining candidates
        %       round_num      - Updated round counter
        %       current_cv_bps - Updated CV bits/spike
        %       phase_history  - Struct array with round details for this phase

            phase_history = struct('round', {}, 'phase', {}, 'tested', {}, 'delta_bps', {}, ...
                'best_candidate', {}, 'added', {}, 'cv_bps_after', {});
            
            while ~isempty(remaining)
                round_num = round_num + 1;
                
                best_delta = -Inf;
                best_candidate = '';
                tested_results = struct('candidate', {}, 'delta_bps', {}, 'cv_bps', {});
                
                for ci = 1:length(remaining)
                    candidate = remaining{ci};
                    
                    % Build model with current selected + this candidate
                    test_vars = [selected_vars, {candidate}];
                    [X_test, ~] = glm_helpers.assemble_design_matrix_selected(B_speed, B_tf, B_onset, ...
                        sf_v, or_v, test_vars, sf_ref_levels, or_ref_levels);
                    
                    % Check if design matrix is valid
                    if size(X_test, 2) >= length(y)
                        tested_results(ci).candidate = candidate;
                        tested_results(ci).delta_bps = -Inf;
                        tested_results(ci).cv_bps = -Inf;
                        continue;
                    end
                    
                    % Cross-validate
                    [~, test_cv_bps] = glm_helpers.cross_validate_glm(X_test, y, offset, fold_ids, 0);
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
                entry_idx = length(phase_history) + 1;
                phase_history(entry_idx).round = round_num;
                phase_history(entry_idx).phase = phase_num;
                phase_history(entry_idx).tested = tested_results;
                phase_history(entry_idx).best_candidate = best_candidate;
                phase_history(entry_idx).delta_bps = best_delta;
                
                % Diagnostic logging
                fprintf('    Round %d (%s): ', round_num, phase_label);
                for ti = 1:length(tested_results)
                    fprintf('%s=%.4f ', tested_results(ti).candidate, tested_results(ti).delta_bps);
                end
                fprintf('\n');
                fprintf('    Best: %s (delta bps=%.4f), threshold=%.4f\n', best_candidate, best_delta, threshold);
                
                % Warn about extreme negative deltas (numerical instability)
                if best_delta < -1
                    fprintf('    WARNING: Large negative delta (%.2f) suggests numerical instability\n', best_delta);
                end
                
                % Add best candidate if it exceeds threshold
                if best_delta > threshold && ~isempty(best_candidate)
                    fprintf('    -> ADDED %s\n', best_candidate);
                    selected_vars{end+1} = best_candidate; %#ok<AGROW>
                    remaining = setdiff(remaining, {best_candidate});
                    % Update current_cv_bps
                    for ti = 1:length(tested_results)
                        if strcmp(tested_results(ti).candidate, best_candidate)
                            current_cv_bps = tested_results(ti).cv_bps;
                            break;
                        end
                    end
                    phase_history(entry_idx).added = true;
                    phase_history(entry_idx).cv_bps_after = current_cv_bps;
                else
                    phase_history(entry_idx).added = false;
                    phase_history(entry_idx).cv_bps_after = current_cv_bps;
                    break;  % Stop this phase: no candidate improved enough
                end
            end
        end
        
    end
end
