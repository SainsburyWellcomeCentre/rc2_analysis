classdef ModelSelectionTuning < handle
% ModelSelectionTuning Class for model selection and significance testing
%
%   This class implements a two-step approach:
%     1. Fit multiple candidate models to data and select best using R² on mean firing rates
%     2. Generate null distribution by shuffling data for the selected model
%
%   Candidate Models:
%     - Linear:              b0 + b1*x
%     - Quadratic:           b0 + b1*x + b2*x^2
%     - Cubic:               b0 + b1*x + b2*x^2 + b3*x^3
%     - Gaussian:            A*exp(-(x-mu)^2/(2*sigma^2)) + baseline
%     - Asymmetric Gaussian: R_max*exp(-(x-x_max)^2/sigma_x) where sigma_x differs for x < x_max and x >= x_max
%     - Sigmoid:             A/(1 + exp(-k*(x-x0))) + baseline
%
%   ModelSelectionTuning Properties:
%         n_reps          - number of shuffles to perform
%         x               - # bins x 1 vector giving bin centers
%         tuning          - # bins x # trials giving firing rate per bin per trial
%         n_bins          - # bins
%         n_trials        - # trials
%         
%         models          - struct array with all fitted models
%         best_model      - name of best model (highest R² on mean)
%         best_model_info - struct with fit details of best model
%         
%         bic_shuff       - n_reps x 1 vector, BIC values for shuffled data
%         fit_metric_shuff - n_reps x 1 vector, R² values for shuffled data
%         shuff_tuning    - # bins x n_reps containing average for each shuffle
%         shuff_sd        - # bins x n_reps containing SD for each shuffle
%         shuff_n         - # bins x n_reps containing n trials averaged for each shuffle
%         
%         p               - p-value from the shuffling procedure
%
%   ModelSelectionTuning Methods:
%         fit_all_models         - fit all candidate models to data
%         select_best_model      - select model with highest R² on mean firing rates
%         get_shuffled_null      - perform shuffling for best model
%         get_summary            - return properties as a structure
%
%   See also: AsymmetricGaussianFit, ShuffleTuning, TuningTable, TFTuningTable

    properties
        
        n_reps = 1000
        force_linear = false  % if true, only fit linear model (no model selection)
        x
        tuning
        n_bins
        n_trials
        
        models              % struct array with all model fits
        best_model          % name of best model
        best_model_info     % detailed info about best model
        
        bic_shuff
        fit_metric_shuff
        shuff_tuning
        shuff_sd
        shuff_n
        
        p
        p_lm                % p-value from linear model (for backward compatibility)
    end
    
    
    
    methods
        
        function obj = ModelSelectionTuning(tuning, x, varargin)
        % ModelSelectionTuning Constructor
        %
        %   ModelSelectionTuning(FIRING_RATE, X)
        %   creates an object for model selection and significance testing.
        %
        %   ModelSelectionTuning(FIRING_RATE, X, 'force_linear', true)
        %   forces linear model only (no model selection), for backward
        %   compatibility with ShuffleTuning and simple velocity/acceleration tuning.
        %
        %   FIRING_RATE is a # bins x # trials matrix with the
        %   firing rates of a cluster for a set of trials. X is a # bins x 1 
        %   vector with the centers of the bins.
        %
        %   The constructor:
        %     1. Fits all candidate models to the data (or only linear if force_linear=true)
        %     2. Selects the best model using R² on mean firing rates (or uses linear if force_linear=true)
        %     3. Generates null distribution by shuffling for the best model (stores R²)
        %     4. Computes p-value by comparing true R² to shuffled R² distribution
        
            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'force_linear', false, @islogical);
            parse(p, varargin{:});
            
            obj.force_linear = p.Results.force_linear;
            obj.x = x(:); % n bins x 1
            obj.tuning = tuning; % n bins x n trials
            obj.n_bins = size(obj.tuning, 1);
            obj.n_trials = size(obj.tuning, 2);
            
            % Step 1: Fit all models (or only linear)
            obj.fit_all_models();
            
            % Step 2: Select best model (or use linear)
            obj.select_best_model();
            
            % Step 3: Generate null distribution for best model
            obj.get_shuffled_null();
            
            % Step 4: Compute p-value using R²
            if isnan(obj.best_model_info.rsq)
                obj.p = nan;
                obj.p_lm = nan;
            else
                % One-tailed test: is the true R² better than shuffled?
                % For R², higher is better
                obj.p = sum(obj.fit_metric_shuff >= obj.best_model_info.rsq) / obj.n_reps;
                
                % For backward compatibility, also compute p_lm for linear model
                if strcmp(obj.best_model, 'linear')
                    obj.p_lm = obj.p;
                else
                    obj.p_lm = nan;
                end
            end
        end
        
        
        
        function fit_all_models(obj)
        %%fit_all_models Fit all candidate models to the data
        %
        %   fit_all_models() fits linear, Gaussian, asymmetric Gaussian,
        %   and sigmoid models to the data and stores results.
        %   If force_linear is true, only fits the linear model.
        
            x_ = repmat(obj.x, 1, obj.n_trials);
            y_ = obj.tuning(:);
            
            % Remove NaN values
            valid = ~isnan(y_);
            x_valid = x_(valid);
            y_valid = y_(valid);
            
            n_valid = length(y_valid);
            
            % Initialize models struct
            obj.models = struct();
            
            % 1. Linear model (always fit)
            obj.models.linear = obj.fit_polynomial(x_valid, y_valid, n_valid, 1);
            
            % If force_linear is true, skip other models
            if obj.force_linear
                return;
            end
            
            % 2. Quadratic model
            obj.models.quadratic = obj.fit_polynomial(x_valid, y_valid, n_valid, 2);
            
            % 3. Cubic model
            obj.models.cubic = obj.fit_polynomial(x_valid, y_valid, n_valid, 3);
            
            % 4. Gaussian model
            obj.models.gaussian = obj.fit_gaussian(x_valid, y_valid, n_valid);
            
            % 5. Asymmetric Gaussian model
            obj.models.asymmetric_gaussian = obj.fit_asymmetric_gaussian(x_valid, y_valid, n_valid);
            
            % 6. Sigmoid model
            obj.models.sigmoid = obj.fit_sigmoid(x_valid, y_valid, n_valid);
        end
        
        
        
        function select_best_model(obj)
        %%select_best_model Select model with highest R² on mean firing rates
        %
        %   select_best_model() compares R² values calculated on mean firing rates
        %   across bins and selects the model with the highest R².
        %   If force_linear is true, simply uses the linear model.
        
            if obj.force_linear
                % Use linear model without selection
                obj.best_model = 'linear';
                obj.best_model_info = obj.models.linear;
                return;
            end
            
            % Calculate mean firing rate for each bin
            mean_firing_rates = nanmean(obj.tuning, 2);
            
            model_names = fieldnames(obj.models);
            rsq_on_mean_values = zeros(length(model_names), 1);
            
            % Calculate R² on mean for each model
            for i = 1:length(model_names)
                model = obj.models.(model_names{i});
                
                % Get fitted values for the mean
                switch model.name
                    case {'linear', 'quadratic', 'cubic'}
                        y_fit = polyval(model.beta, obj.x);
                        
                    case 'gaussian'
                        % params = [amplitude, mu, sigma, baseline]
                        params = model.beta;
                        y_fit = params(1) * exp(-(obj.x - params(2)).^2 / (2*params(3)^2)) + params(4);
                        
                    case 'asymmetric_gaussian'
                        % params = [R_max, x_max, sigma_minus, sigma_plus]
                        params = model.beta;
                        y_fit = ModelSelectionTuning.asymmetric_gaussian_static(params, obj.x);
                        
                    case 'sigmoid'
                        % params = [amplitude, k, x0, baseline]
                        params = model.beta;
                        y_fit = params(1) ./ (1 + exp(-params(2) * (obj.x - params(3)))) + params(4);
                        
                    otherwise
                        rsq_on_mean_values(i) = -inf;
                        continue;
                end
                
                % Calculate R² on mean firing rates
                SS_res = sum((mean_firing_rates - y_fit).^2);
                SS_tot = sum((mean_firing_rates - nanmean(mean_firing_rates)).^2);
                rsq_on_mean_values(i) = 1 - (SS_res / SS_tot);
            end
            
            % Select model with maximum R² on mean (higher is better)
            [~, max_idx] = max(rsq_on_mean_values);
            obj.best_model = model_names{max_idx};
            obj.best_model_info = obj.models.(obj.best_model);
            
            % Store the R² on mean in the best model info
            obj.best_model_info.rsq_on_mean = rsq_on_mean_values(max_idx);
        end
        
        
        
        function get_shuffled_null(obj)
        %%get_shuffled_null Generate null distribution for best model
        %
        %   get_shuffled_null() shuffles the tuning matrix and refits
        %   the best model n_reps times to create a null distribution of R².
        %   Uses parallel computing for efficiency (uses existing pool).
        
            rng(1);  % For reproducibility
            
            % Use existing parallel pool (don't create/restart here)
            % Pool should be started once at script level
            
            x_ = repmat(obj.x, 1, obj.n_trials);
            
            % Preallocate arrays
            bic_shuff = nan(1, obj.n_reps);
            fit_metric_shuff = nan(1, obj.n_reps);
            shuff_tuning = nan(obj.n_bins, obj.n_reps);
            shuff_sd = nan(obj.n_bins, obj.n_reps);
            shuff_n = nan(obj.n_bins, obj.n_reps);
            
            % Pre-extract data for parallel access
            tuning_data = obj.tuning;
            x_data = obj.x;
            best_model_name = obj.best_model;
            
            % Track failures for error reporting
            fit_failures = cell(1, obj.n_reps);
            
            parfor rand_i = 1:obj.n_reps
                try
                    % Shuffle by randomly sampling from all values
                    I = randi(numel(tuning_data), size(tuning_data));
                    new_tuning = tuning_data(I);
                    
                    % Flatten for fitting
                    x_rep = repmat(x_data, 1, size(new_tuning, 2));
                    y_shuff = new_tuning(:);
                    valid = ~isnan(y_shuff);
                    x_valid = x_rep(valid);
                    y_valid = y_shuff(valid);
                    n_valid = length(y_valid);
                    
                    % Fit the best model to shuffled data
                    switch best_model_name
                        case 'linear'
                            fit_result = ModelSelectionTuning.fit_polynomial(x_valid, y_valid, n_valid, 1);
                        case 'quadratic'
                            fit_result = ModelSelectionTuning.fit_polynomial(x_valid, y_valid, n_valid, 2);
                        case 'cubic'
                            fit_result = ModelSelectionTuning.fit_polynomial(x_valid, y_valid, n_valid, 3);
                        case 'gaussian'
                            fit_result = ModelSelectionTuning.fit_gaussian(x_valid, y_valid, n_valid);
                        case 'asymmetric_gaussian'
                            fit_result = ModelSelectionTuning.fit_asymmetric_gaussian(x_valid, y_valid, n_valid);
                        case 'sigmoid'
                            fit_result = ModelSelectionTuning.fit_sigmoid(x_valid, y_valid, n_valid);
                        otherwise
                            error('Unknown model type: %s', best_model_name);
                    end
                    
                    bic_shuff(rand_i) = fit_result.bic;
                    fit_metric_shuff(rand_i) = fit_result.rsq;
                    
                    % Store shuffled tuning statistics
                    shuff_tuning(:, rand_i) = nanmean(new_tuning, 2);
                    shuff_sd(:, rand_i) = nanstd(new_tuning, [], 2);
                    shuff_n(:, rand_i) = sum(~isnan(new_tuning), 2);
                    
                catch ME
                    % Store error information
                    fit_failures{rand_i} = sprintf('Shuffle %d failed: %s', rand_i, ME.message);
                    % Fill with NaN to continue
                    bic_shuff(rand_i) = nan;
                    fit_metric_shuff(rand_i) = nan;
                    shuff_tuning(:, rand_i) = nan;
                    shuff_sd(:, rand_i) = nan;
                    shuff_n(:, rand_i) = nan;
                end
            end
            
            % Report failures if any occurred
            failures = fit_failures(~cellfun(@isempty, fit_failures));
            if ~isempty(failures)
                warning('ModelSelectionTuning:ShuffleFailures', ...
                    'Shuffling encountered %d failures (out of %d):\n%s', ...
                    length(failures), obj.n_reps, strjoin(failures(1:min(5, length(failures))), '\n'));
                if length(failures) > 5
                    warning('ModelSelectionTuning:ShuffleFailures', ...
                        '... and %d more failures (showing first 5 only)', length(failures) - 5);
                end
            end
            
            % Assign results to object properties after parfor
            obj.bic_shuff = bic_shuff;
            obj.fit_metric_shuff = fit_metric_shuff;
            obj.shuff_tuning = shuff_tuning;
            obj.shuff_sd = shuff_sd;
            obj.shuff_n = shuff_n;
        end
        
        
        
        function val = get_summary(obj)
        %%get_summary Return essential properties in a structure to save
        %
        %   STRUCT = get_summary()
        %   Returns a structure with all relevant properties.
        %   For backward compatibility, when force_linear=true, also includes
        %   ShuffleTuning-compatible fields (rsq, beta, r, rsq_shuff, etc.)
        
            val.n_reps = obj.n_reps;
            val.best_model = obj.best_model;
            val.best_model_info = obj.best_model_info;
            val.all_models = obj.models;
            val.bic_shuff = obj.bic_shuff;
            val.fit_metric_shuff = obj.fit_metric_shuff;
            val.shuff_tuning = obj.shuff_tuning;
            val.shuff_sd = obj.shuff_sd;
            val.shuff_n = obj.shuff_n;
            val.p = obj.p;
            
            % Backward compatibility: add ShuffleTuning-like fields for linear model
            if strcmp(obj.best_model, 'linear')
                val.rsq = obj.best_model_info.rsq;
                val.beta = obj.best_model_info.beta;
                val.r = obj.best_model_info.r;
                val.rsq_shuff = obj.fit_metric_shuff;  % Same as fit_metric for linear
                val.p_lm = obj.p_lm;
            end
            
            % Backward compatibility: add AsymmetricGaussianFit-like fields for asymmetric_gaussian model
            if strcmp(obj.best_model, 'asymmetric_gaussian')
                params = obj.best_model_info.beta;
                val.R_max = params(1);
                val.x_max = params(2);
                val.sigma_minus = params(3);
                val.sigma_plus = params(4);
                val.rsq = obj.best_model_info.rsq;
                val.rsq_shuff = obj.fit_metric_shuff;
                % Store parameters from all shuffles for plotting
                val.params_shuff = nan(obj.n_reps, 4);  % Not stored in current implementation
            end
        end
    end
    
    
    
    methods (Static = true)
        
        function model_info = fit_polynomial(x, y, n, degree)
        %%fit_polynomial Fit polynomial model of specified degree
        %
        %   MODEL = fit_polynomial(X, Y, N, DEGREE)
        
            beta = polyfit(x, y, degree);
            y_pred = polyval(beta, x);
            
            % Compute fit metrics
            residuals = y - y_pred;
            rss = sum(residuals.^2);
            tss = sum((y - mean(y)).^2);
            rsq = 1 - rss/tss;
            
            % Log-likelihood
            sigma2 = rss / n;
            if sigma2 <= 0
                log_likelihood = -inf;
            else
                log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
            end
            
            % BIC
            k = degree + 1;  % number of parameters
            bic = -2*log_likelihood + k*log(n);
            
            % Correlation
            r = corr(x, y);
            
            % Store results
            switch degree
                case 1
                    name = 'linear';
                case 2
                    name = 'quadratic';
                case 3
                    name = 'cubic';
                otherwise
                    name = sprintf('poly%d', degree);
            end
            
            model_info.name = name;
            model_info.beta = beta;
            model_info.rsq = rsq;
            model_info.r = r;
            model_info.bic = bic;
            model_info.log_likelihood = log_likelihood;
            model_info.n_params = k;
        end
        
        
        
        function model_info = fit_gaussian(x, y, n)
        %%fit_gaussian Fit Gaussian model: A·exp(-(x-μ)²/(2σ²)) + baseline
        %
        %   MODEL = fit_gaussian(X, Y, N)
        
            % Initialize parameters from data
            [max_y, max_idx] = max(y);
            min_y = min(y);
            amplitude_init = max_y - min_y;
            mu_init = x(max_idx);
            sigma_init = (max(x) - min(x)) / 4;  % quarter of range
            baseline_init = min_y;
            
            % Define Gaussian function
            gauss_fun = @(params, x) params(1) * exp(-(x - params(2)).^2 / (2*params(3)^2)) + params(4);
            
            % Initial guess and bounds
            params_init = [amplitude_init, mu_init, sigma_init, baseline_init];
            lb = [0, min(x), 0.01, -inf];  % Lower bounds
            ub = [inf, max(x), inf, inf];  % Upper bounds
            
            % Fit using nonlinear least squares
            try
                options = optimoptions('lsqcurvefit', 'Display', 'off');
                params = lsqcurvefit(gauss_fun, params_init, x, y, lb, ub, options);
                y_pred = gauss_fun(params, x);
                
                % Compute fit metrics
                residuals = y - y_pred;
                rss = sum(residuals.^2);
                tss = sum((y - mean(y)).^2);
                rsq = 1 - rss/tss;
                
                % Log-likelihood
                sigma2 = rss / n;
                if sigma2 <= 0
                    log_likelihood = -inf;
                else
                    log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
                end
                
                % BIC
                k = 4;  % 4 parameters (amplitude, mu, sigma, baseline)
                bic = -2*log_likelihood + k*log(n);
                
            catch
                % If fit fails, return poor metrics
                params = params_init;
                rsq = -inf;
                bic = inf;
                log_likelihood = -inf;
            end
            
            % Store results
            model_info.name = 'gaussian';
            model_info.beta = params;  % [amplitude, mu, sigma, baseline]
            model_info.rsq = rsq;
            model_info.bic = bic;
            model_info.log_likelihood = log_likelihood;
            model_info.n_params = 4;
        end
        
        
        
        function model_info = fit_asymmetric_gaussian(x, y, n)
        %%fit_asymmetric_gaussian Fit asymmetric Gaussian: R_max·exp(-(x-x_max)²/σ_x)
        %
        %   MODEL = fit_asymmetric_gaussian(X, Y, N)
        %   where σ_x = sigma_minus when x < x_max, sigma_plus when x ≥ x_max
        
            % Initialize parameters from data
            [R_max_init, idx_max] = max(y);
            x_max_init = x(idx_max);
            sigma_init = abs(x(end) - x(1)) / 4;
            
            % Define asymmetric Gaussian function
            asym_gauss_fun = @(params, x) ModelSelectionTuning.asymmetric_gaussian_static(params, x);
            
            % Parameter bounds: [R_max, x_max, sigma_minus, sigma_plus]
            params_init = [R_max_init, x_max_init, sigma_init, sigma_init];
            lb = [0, min(x), 0.001, 0.001];
            ub = [max(y) * 2, max(x), abs(x(end) - x(1)) * 2, abs(x(end) - x(1)) * 2];
            
            % Fit using nonlinear least squares
            try
                options = optimoptions('lsqcurvefit', 'Display', 'off');
                params = lsqcurvefit(asym_gauss_fun, params_init, x, y, lb, ub, options);
                y_pred = asym_gauss_fun(params, x);
                
                % Compute fit metrics
                residuals = y - y_pred;
                rss = sum(residuals.^2);
                tss = sum((y - mean(y)).^2);
                rsq = 1 - rss/tss;
                
                % Log-likelihood
                sigma2 = rss / n;
                if sigma2 <= 0
                    log_likelihood = -inf;
                else
                    log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
                end
                
                % BIC
                k = 4;  % 4 parameters (R_max, x_max, sigma_minus, sigma_plus)
                bic = -2*log_likelihood + k*log(n);
                
            catch
                % If fit fails, return poor metrics
                params = params_init;
                rsq = -inf;
                bic = inf;
                log_likelihood = -inf;
            end
            
            % Store results
            model_info.name = 'asymmetric_gaussian';
            model_info.beta = params;  % [R_max, x_max, sigma_minus, sigma_plus]
            model_info.rsq = rsq;
            model_info.bic = bic;
            model_info.log_likelihood = log_likelihood;
            model_info.n_params = 4;
        end
        
        
        
        function model_info = fit_sigmoid(x, y, n)
        %%fit_sigmoid Fit sigmoid model: A/(1 + exp(-k·(x-x₀))) + baseline
        %
        %   MODEL = fit_sigmoid(X, Y, N)
        
            % Initialize parameters from data
            max_y = max(y);
            min_y = min(y);
            amplitude_init = max_y - min_y;
            baseline_init = min_y;
            x0_init = mean(x);  % midpoint
            k_init = 1 / std(x);  % steepness
            
            % Define sigmoid function
            sigmoid_fun = @(params, x) params(1) ./ (1 + exp(-params(2) * (x - params(3)))) + params(4);
            
            % Initial guess and bounds
            params_init = [amplitude_init, k_init, x0_init, baseline_init];
            lb = [0, 0, min(x), -inf];
            ub = [inf, inf, max(x), inf];
            
            % Fit using nonlinear least squares
            try
                options = optimoptions('lsqcurvefit', 'Display', 'off');
                params = lsqcurvefit(sigmoid_fun, params_init, x, y, lb, ub, options);
                y_pred = sigmoid_fun(params, x);
                
                % Compute fit metrics
                residuals = y - y_pred;
                rss = sum(residuals.^2);
                tss = sum((y - mean(y)).^2);
                rsq = 1 - rss/tss;
                
                % Log-likelihood
                sigma2 = rss / n;
                if sigma2 <= 0
                    log_likelihood = -inf;
                else
                    log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
                end
                
                % BIC
                k = 4;  % 4 parameters
                bic = -2*log_likelihood + k*log(n);
                
            catch
                % If fit fails, return poor metrics
                params = params_init;
                rsq = -inf;
                bic = inf;
            end
            
            % Store results
            model_info.name = 'sigmoid';
            model_info.beta = params;  % [amplitude, k, x0, baseline]
            model_info.rsq = rsq;
            model_info.bic = bic;
            model_info.log_likelihood = log_likelihood;
            model_info.n_params = 4;
        end
        
        
        
        function y = asymmetric_gaussian_static(params, x)
        %%asymmetric_gaussian_static Static asymmetric Gaussian function
        %
        %   Y = asymmetric_gaussian_static(PARAMS, X)
        %   PARAMS = [R_max, x_max, sigma_minus, sigma_plus]
        %   
        %   The function is:
        %       R = R_max * exp(-(x - x_max)^2 / sigma_x)
        %   where sigma_x = sigma_minus when x < x_max
        %         sigma_x = sigma_plus when x >= x_max
        
            R_max = params(1);
            x_max = params(2);
            sigma_minus = params(3);
            sigma_plus = params(4);
            
            % Initialize output
            y = zeros(size(x));
            
            % Left side (x < x_max): use sigma_minus
            left_idx = x < x_max;
            y(left_idx) = R_max * exp(-(x(left_idx) - x_max).^2 / sigma_minus);
            
            % Right side (x >= x_max): use sigma_plus
            right_idx = x >= x_max;
            y(right_idx) = R_max * exp(-(x(right_idx) - x_max).^2 / sigma_plus);
        end
    end
end
