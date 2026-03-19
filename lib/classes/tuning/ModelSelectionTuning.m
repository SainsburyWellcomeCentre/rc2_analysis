classdef ModelSelectionTuning < handle
% ModelSelectionTuning Class for model selection and significance testing
%
%   This class implements a two-step approach:
%     1. Fit multiple candidate models to data and select best using BIC
%     2. Generate null distribution by shuffling data for the selected model
%
%   Candidate Models:
%     - Linear:    β₀ + β₁·x
%     - Quadratic: β₀ + β₁·x + β₂·x²
%     - Cubic:     β₀ + β₁·x + β₂·x² + β₃·x³
%     - Gaussian:  A·exp(-(x-μ)²/(2σ²)) + baseline
%     - ReLU:      max(0, β₀ + β₁·x)
%     - Sigmoid:   A/(1 + exp(-k·(x-x₀))) + baseline
%
%   ModelSelectionTuning Properties:
%         n_reps          - number of shuffles to perform
%         x               - # bins x 1 vector giving bin centers
%         tuning          - # bins x # trials giving firing rate per bin per trial
%         n_bins          - # bins
%         n_trials        - # trials
%         
%         models          - struct array with all fitted models
%         best_model      - name of best model (lowest BIC)
%         best_model_info - struct with fit details of best model
%         
%         bic_shuff       - n_reps x 1 vector, BIC values for shuffled data
%         fit_metric_shuff - n_reps x 1 vector, fit metric (log-likelihood) for shuffled
%         shuff_tuning    - # bins x n_reps containing average for each shuffle
%         shuff_sd        - # bins x n_reps containing SD for each shuffle
%         shuff_n         - # bins x n_reps containing n trials averaged for each shuffle
%         
%         p               - p-value from the shuffling procedure
%
%   ModelSelectionTuning Methods:
%         fit_all_models         - fit all candidate models to data
%         select_best_model      - select model with lowest BIC
%         get_shuffled_null      - perform shuffling for best model
%         get_summary            - return properties as a structure
%
%   See also: ShuffleTuning, TuningTable, VelocityBins, VelocityTuningCurve

    properties
        
        n_reps = 1000
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
    end
    
    
    
    methods
        
        function obj = ModelSelectionTuning(tuning, x)
        % ModelSelectionTuning Constructor
        %
        %   ModelSelectionTuning(FIRING_RATE, X)
        %   creates an object for model selection and significance testing.
        %
        %   FIRING_RATE is a # bins x # trials matrix with the
        %   firing rates of a cluster for a set of trials. X is a # bins x 1 
        %   vector with the centers of the bins.
        %
        %   The constructor:
        %     1. Fits all candidate models to the data
        %     2. Selects the best model using BIC
        %     3. Generates null distribution by shuffling for the best model
        %     4. Computes p-value by comparing true fit to null distribution
        
            obj.x = x(:); % n bins x 1
            obj.tuning = tuning; % n bins x n trials
            obj.n_bins = size(obj.tuning, 1);
            obj.n_trials = size(obj.tuning, 2);
            
            % Step 1: Fit all models
            obj.fit_all_models();
            
            % Step 2: Select best model
            obj.select_best_model();
            
            % Step 3: Generate null distribution for best model
            obj.get_shuffled_null();
            
            % Step 4: Compute p-value
            if isnan(obj.best_model_info.fit_metric)
                obj.p = nan;
            else
                % One-tailed test: is the true fit better than shuffled?
                obj.p = sum(obj.fit_metric_shuff >= obj.best_model_info.fit_metric) / obj.n_reps;
            end
        end
        
        
        
        function fit_all_models(obj)
        %%fit_all_models Fit all candidate models to the data
        %
        %   fit_all_models() fits flat, linear, quadratic, cubic, Gaussian,
        %   ReLU, and sigmoid models to the data and stores results.
        
            x_ = repmat(obj.x, 1, obj.n_trials);
            y_ = obj.tuning(:);
            
            % Remove NaN values
            valid = ~isnan(y_);
            x_valid = x_(valid);
            y_valid = y_(valid);
            
            n_valid = length(y_valid);
            
            % Initialize models struct
            obj.models = struct();
            
            % 1. Linear model
            obj.models.linear = obj.fit_polynomial(x_valid, y_valid, n_valid, 1);
            
            % 2. Quadratic model
            obj.models.quadratic = obj.fit_polynomial(x_valid, y_valid, n_valid, 2);
            
            % 3. Cubic model
            obj.models.cubic = obj.fit_polynomial(x_valid, y_valid, n_valid, 3);
            
            % 4. Gaussian model
            obj.models.gaussian = obj.fit_gaussian(x_valid, y_valid, n_valid);
            
            % 5. ReLU model
            obj.models.relu = obj.fit_relu(x_valid, y_valid, n_valid);
            
            % 6. Sigmoid model
            obj.models.sigmoid = obj.fit_sigmoid(x_valid, y_valid, n_valid);
        end
        
        
        
        function select_best_model(obj)
        %%select_best_model Select model with lowest BIC
        %
        %   select_best_model() compares BIC values and selects the best model.
        
            model_names = fieldnames(obj.models);
            bic_values = zeros(length(model_names), 1);
            
            for i = 1:length(model_names)
                bic_values(i) = obj.models.(model_names{i}).bic;
            end
            
            % Select model with minimum BIC (lower is better)
            [min_bic, min_idx] = min(bic_values);
            obj.best_model = model_names{min_idx};
            obj.best_model_info = obj.models.(obj.best_model);
            
            % Store BIC difference from linear model for reference
            obj.best_model_info.delta_bic_from_linear = min_bic - obj.models.linear.bic;
        end
        
        
        
        function get_shuffled_null(obj)
        %%get_shuffled_null Generate null distribution for best model
        %
        %   get_shuffled_null() shuffles the tuning matrix and refits
        %   the best model n_reps times to create a null distribution.
        
            rng(1);  % For reproducibility
            
            x_ = repmat(obj.x, 1, obj.n_trials);
            
            % Preallocate arrays
            obj.bic_shuff = nan(1, obj.n_reps);
            obj.fit_metric_shuff = nan(1, obj.n_reps);
            obj.shuff_tuning = nan(obj.n_bins, obj.n_reps);
            obj.shuff_sd = nan(obj.n_bins, obj.n_reps);
            obj.shuff_n = nan(obj.n_bins, obj.n_reps);
            
            for rand_i = 1:obj.n_reps
                % Shuffle by randomly sampling from all values
                I = randi(numel(obj.tuning), size(obj.tuning));
                new_tuning = obj.tuning(I);
                
                % Flatten for fitting
                y_shuff = new_tuning(:);
                valid = ~isnan(y_shuff);
                x_valid = x_(valid);
                y_valid = y_shuff(valid);
                n_valid = length(y_valid);
                
                % Fit the best model to shuffled data
                switch obj.best_model
                    case 'linear'
                        fit_result = obj.fit_polynomial(x_valid, y_valid, n_valid, 1);
                    case 'quadratic'
                        fit_result = obj.fit_polynomial(x_valid, y_valid, n_valid, 2);
                    case 'cubic'
                        fit_result = obj.fit_polynomial(x_valid, y_valid, n_valid, 3);
                    case 'gaussian'
                        fit_result = obj.fit_gaussian(x_valid, y_valid, n_valid);
                    case 'relu'
                        fit_result = obj.fit_relu(x_valid, y_valid, n_valid);
                    case 'sigmoid'
                        fit_result = obj.fit_sigmoid(x_valid, y_valid, n_valid);
                    otherwise
                        error('Unknown model type: %s', obj.best_model);
                end
                
                obj.bic_shuff(rand_i) = fit_result.bic;
                obj.fit_metric_shuff(rand_i) = fit_result.fit_metric;
                
                % Store shuffled tuning statistics
                obj.shuff_tuning(:, rand_i) = nanmean(new_tuning, 2);
                obj.shuff_sd(:, rand_i) = nanstd(new_tuning, [], 2);
                obj.shuff_n(:, rand_i) = sum(~isnan(new_tuning), 2);
            end
        end
        
        
        
        function val = get_summary(obj)
        %%get_summary Return essential properties in a structure to save
        %
        %   STRUCT = get_summary()
        
            val.n_reps = obj.n_reps;
            val.best_model = obj.best_model;
            val.best_model_info = obj.best_model_info;
            val.all_models = obj.models;
            val.bic_shuff = obj.bic_shuff;
            val.fit_metric_shuff = obj.fit_metric_shuff;
            val.shuff_tuning = obj.shuff_tuning;
            val.p = obj.p;
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
            log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
            
            % BIC
            k = degree + 1;  % number of parameters
            bic = -2*log_likelihood + k*log(n);
            
            % Correlation
            r = corr(polyval(beta, x), y);  % correlation between predicted and actual
            model_info.fit_metric = r;
            
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
            model_info.log_likelihood = log_likelihood;
            model_info.bic = bic;
            model_info.n_params = k;
            model_info.fit_metric = r;
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
                log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
                
                % BIC
                k = 4;  % 4 parameters (amplitude, mu, sigma, baseline)
                bic = -2*log_likelihood + k*log(n);
                %r
                r = corr(y_pred(:), y(:));
                model_info.r = r;
                model_info.fit_metric = r;
                
            catch
                % If fit fails, return poor metrics
                params = params_init;
                rsq = -inf;
                bic = inf;
                log_likelihood = -inf;
                r = -inf;  %
            end
            
            % Store results
            model_info.name = 'gaussian';
            model_info.beta = params;  % [amplitude, mu, sigma, baseline]
            model_info.rsq = rsq;
            model_info.log_likelihood = log_likelihood;
            model_info.bic = bic;
            model_info.n_params = 4;
            model_info.fit_metric = r;
        end
        
        
        
        function model_info = fit_relu(x, y, n)
        %%fit_relu Fit ReLU model: max(0, β₀ + β₁·x)
        %
        %   MODEL = fit_relu(X, Y, N)
        
            % Define ReLU function
            relu_fun = @(params, x) max(0, params(1) + params(2) * x);
            
            % Initial guess from linear fit
            beta_linear = polyfit(x, y, 1);
            params_init = [beta_linear(2), beta_linear(1)];  % [intercept, slope]
            
            % Fit using nonlinear least squares
            try
                options = optimoptions('lsqcurvefit', 'Display', 'off');
                params = lsqcurvefit(@(p, x) arrayfun(@(xi) max(0, p(1) + p(2)*xi), x), ...
                                     params_init, x, y, [], [], options);
                y_pred = arrayfun(@(xi) max(0, params(1) + params(2)*xi), x);
                
                % Compute fit metrics
                residuals = y - y_pred;
                rss = sum(residuals.^2);
                tss = sum((y - mean(y)).^2);
                rsq = 1 - rss/tss;
                
                % Log-likelihood
                sigma2 = rss / n;
                log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
                
                % BIC
                k = 2;  % 2 parameters (intercept, slope)
                bic = -2*log_likelihood + k*log(n);

                % r
                r = corr(y_pred(:), y(:));
                model_info.r = r;
                model_info.fit_metric = r;
                
            catch
                % If fit fails, return poor metrics
                params = params_init;
                rsq = -inf;
                bic = inf;
                log_likelihood = -inf;
                r = -inf;
            end
            
            % Store results
            model_info.name = 'relu';
            model_info.beta = params;  % [intercept, slope]
            model_info.rsq = rsq;
            model_info.log_likelihood = log_likelihood;
            model_info.bic = bic;
            model_info.n_params = 2;
            model_info.fit_metric = r;
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
                log_likelihood = -0.5 * n * (log(2*pi*sigma2) + 1);
                
                % BIC
                k = 4;  % 4 parameters
                bic = -2*log_likelihood + k*log(n);
                %r
                r = corr(y_pred(:), y(:));
                model_info.r = r;
                model_info.fit_metric = r;
                
            catch
                % If fit fails, return poor metrics
                params = params_init;
                rsq = -inf;
                bic = inf;
                log_likelihood = -inf;
                r = -inf;
            end
            
            % Store results
            model_info.name = 'sigmoid';
            model_info.beta = params;  % [amplitude, k, x0, baseline]
            model_info.rsq = rsq;
            model_info.log_likelihood = log_likelihood;
            model_info.bic = bic;
            model_info.n_params = 4;
            model_info.fit_metric = r;
        end
    end
end
