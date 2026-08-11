classdef AsymmetricGaussianFit < handle
% AsymmetricGaussianFit Class for fitting asymmetric Gaussian tuning curves
%
%   Fits temporal frequency tuning data with an asymmetric Gaussian function:
%       R = R_max * exp(-(x - x_max)^2 / sigma_x)
%   where sigma_x = sigma_minus when x < x_max, and sigma_x = sigma_plus when x >= x_max
%
%   This is the descriptive function from Angelucci et al. papers for characterizing
%   spatial and temporal frequency tuning curves.
%
%   AsymmetricGaussianFit Properties:
%       n_reps          - number of shuffles to perform (default: 1000)
%       x               - # bins x 1 vector giving TF bin centers
%       tuning          - # bins x # trials giving firing rate in a TF bin for each trial
%       n_bins          - # bins
%       n_trials        - # trials
%       
%       % Fit parameters for true data
%       R_max           - maximum response (peak firing rate)
%       x_max           - preferred temporal frequency
%       sigma_minus     - width parameter for x < x_max
%       sigma_plus      - width parameter for x >= x_max
%       rsq             - R^2 goodness of fit
%       rmse            - root mean squared error
%       
%       % Shuffled fit results
%       rsq_shuff       - n_reps x 1 vector, R^2 values for shuffled data
%       params_shuff    - n_reps x 4 matrix, [R_max, x_max, sigma_minus, sigma_plus] for each shuffle
%       shuff_tuning    - # bins x n_reps containing average for each shuffle
%       shuff_sd        - # bins x n_reps containing SD for each shuffle
%       shuff_n         - # bins x n_reps containing n trials averaged for each shuffle
%       
%       p               - p-value from the shuffling procedure
%
%   AsymmetricGaussianFit Methods:
%       fit_data            - fit the asymmetric Gaussian to the data
%       get_shuffled_fits   - perform shuffling and fit each shuffle
%       get_summary         - return properties as a structure
%       predict             - predict firing rate at given x values
%
%   See also: ShuffleTuning, TFTuningTable

    properties
        
        n_reps = 1000
        x
        tuning
        n_bins
        n_trials
        
        % True fit parameters
        R_max
        x_max
        sigma_minus
        sigma_plus
        rsq
        rmse
        
        % Shuffled results
        rsq_shuff
        params_shuff
        shuff_tuning
        shuff_sd
        shuff_n
        
        p
    end
    
    
    
    methods
        
        function obj = AsymmetricGaussianFit(tuning, x)
        % AsymmetricGaussianFit
        %
        %   AsymmetricGaussianFit(FIRING_RATE, TF)
        %   creates an object for fitting asymmetric Gaussian tuning curves
        %   and shuffling the data to determine significance of tuning.
        %
        %   FIRING_RATE is a # TF bins x # trials matrix with the
        %   firing rates of a cluster to a TF bin for a set of trials.
        %   TF is a # TF bins x 1 vector with the centers of the TF bins
        %   in which the firing rates were computed.
        %
        %   An asymmetric Gaussian fit is performed on the data. Then the
        %   FIRING_RATE matrix is completely shuffled, and the fitting
        %   performed multiple times.
        %   
        %   A p-value is calculated by using the fraction of shuffled R^2 value
        %   above the true R^2 value.
        
            obj.x = x(:); % n bins x 1
            obj.tuning = tuning; % n bins x n trials
            obj.n_bins = size(obj.tuning, 1);
            obj.n_trials = size(obj.tuning, 2);
            
            % Fit to data
            obj.fit_data();
            
            % Perform shuffling
            obj.get_shuffled_fits();
            
            % Calculate p-value
            if isnan(obj.rsq)
                obj.p = nan;
            else
                p_up = sum(obj.rsq_shuff > obj.rsq) / obj.n_reps;
                obj.p = p_up;
            end
        end
        
        
        
        function fit_data(obj)
        %%fit_data Fit the asymmetric Gaussian to the data
        %
        %   fit_data() fits the asymmetric Gaussian function to the mean
        %   tuning curve and stores parameters in object properties.
        
            % Compute mean tuning curve
            mean_tuning = nanmean(obj.tuning, 2);
            
            % Check if we have valid data
            if all(isnan(mean_tuning)) || sum(~isnan(mean_tuning)) < 4
                obj.R_max = nan;
                obj.x_max = nan;
                obj.sigma_minus = nan;
                obj.sigma_plus = nan;
                obj.rsq = nan;
                obj.rmse = nan;
                return;
            end
            
            % Initial parameter guesses
            [R_max_init, idx_max] = max(mean_tuning);
            x_max_init = obj.x(idx_max);
            
            % Estimate initial sigma from half-width at half-maximum
            half_max = R_max_init / 2;
            sigma_init = abs(obj.x(end) - obj.x(1)) / 4; % rough estimate
            
            % Parameter bounds: [R_max, x_max, sigma_minus, sigma_plus]
            lb = [0, min(obj.x), 0.001, 0.001];
            ub = [max(mean_tuning) * 2, max(obj.x), abs(obj.x(end) - obj.x(1)) * 2, abs(obj.x(end) - obj.x(1)) * 2];
            x0 = [R_max_init, x_max_init, sigma_init, sigma_init];
            
            % Fit using lsqcurvefit (nonlinear least squares)
            options = optimoptions('lsqcurvefit', 'Display', 'off');
            
            try
                [params, ~, residual, ~, ~] = lsqcurvefit(@obj.asymmetric_gaussian, x0, obj.x, mean_tuning, lb, ub, options);
                
                obj.R_max = params(1);
                obj.x_max = params(2);
                obj.sigma_minus = params(3);
                obj.sigma_plus = params(4);
                
                % Calculate goodness of fit
                y_fit = obj.asymmetric_gaussian(params, obj.x);
                obj.rmse = sqrt(mean(residual.^2));
                
                % Calculate R^2
                SS_res = sum((mean_tuning - y_fit).^2);
                SS_tot = sum((mean_tuning - nanmean(mean_tuning)).^2);
                obj.rsq = 1 - SS_res / SS_tot;
                
            catch ME
                warning('AsymmetricGaussianFit:fit_data', 'Fitting failed: %s', ME.message);
                obj.R_max = nan;
                obj.x_max = nan;
                obj.sigma_minus = nan;
                obj.sigma_plus = nan;
                obj.rsq = nan;
                obj.rmse = nan;
            end
        end
        
        
        
        function get_shuffled_fits(obj)
        %%get_shuffled_fits Perform the shuffling and fit each shuffle
        %
        %   get_shuffled_fits() takes the data in `tuning` and shuffles
        %   with replacement `n_reps` times. Stores the results in
        %   properties of the object. Parallelized across shuffles.
        
            rng(1);
            
            % Preallocate arrays
            rsq_shuff = nan(1, obj.n_reps);
            params_shuff = nan(obj.n_reps, 4);
            shuff_tuning = nan(obj.n_bins, obj.n_reps);
            shuff_sd = nan(obj.n_bins, obj.n_reps);
            shuff_n = nan(obj.n_bins, obj.n_reps);
            
            % Pre-extract data for parallel access
            tuning_data = obj.tuning;
            x_data = obj.x;
            n_bins = obj.n_bins;
            n_reps = obj.n_reps;
            
            % Parallelize the shuffle loop
            parfor rand_i = 1 : n_reps
                % Shuffle the data
                I = randi(numel(tuning_data), size(tuning_data));
                new_tuning = tuning_data(I);
                
                % Compute mean of shuffled data
                mean_shuffled = nanmean(new_tuning, 2);
                
                % Fit the shuffled data using static method
                [rsq_shuffle, params_shuffle] = AsymmetricGaussianFit.fit_shuffled_data_static(mean_shuffled, x_data, n_bins);
                
                rsq_shuff(rand_i) = rsq_shuffle;
                params_shuff(rand_i, :) = params_shuffle;
                
                shuff_tuning(:, rand_i) = mean_shuffled;
                shuff_sd(:, rand_i) = nanstd(new_tuning, [], 2);
                shuff_n(:, rand_i) = sum(~isnan(new_tuning), 2);
            end
            
            % Assign results to object properties after parfor
            obj.rsq_shuff = rsq_shuff;
            obj.params_shuff = params_shuff;
            obj.shuff_tuning = shuff_tuning;
            obj.shuff_sd = shuff_sd;
            obj.shuff_n = shuff_n;
        end
        
        
        
        function [rsq, params] = fit_shuffled_data(obj, mean_tuning)
        %%fit_shuffled_data Fit a single shuffled dataset
        %
        %   [RSQ, PARAMS] = fit_shuffled_data(MEAN_TUNING)
        
            % Check if we have valid data
            if all(isnan(mean_tuning)) || sum(~isnan(mean_tuning)) < 4
                rsq = nan;
                params = [nan, nan, nan, nan];
                return;
            end
            
            % Initial parameter guesses
            [R_max_init, idx_max] = max(mean_tuning);
            x_max_init = obj.x(idx_max);
            sigma_init = abs(obj.x(end) - obj.x(1)) / 4;
            
            % Parameter bounds
            lb = [0, min(obj.x), 0.001, 0.001];
            ub = [max(mean_tuning) * 2, max(obj.x), abs(obj.x(end) - obj.x(1)) * 2, abs(obj.x(end) - obj.x(1)) * 2];
            x0 = [R_max_init, x_max_init, sigma_init, sigma_init];
            
            % Fit using lsqcurvefit
            options = optimoptions('lsqcurvefit', 'Display', 'off');
            
            try
                [params, ~, ~, ~, ~] = lsqcurvefit(@obj.asymmetric_gaussian, x0, obj.x, mean_tuning, lb, ub, options);
                
                % Calculate R^2
                y_fit = obj.asymmetric_gaussian(params, obj.x);
                SS_res = sum((mean_tuning - y_fit).^2);
                SS_tot = sum((mean_tuning - nanmean(mean_tuning)).^2);
                rsq = 1 - SS_res / SS_tot;
                
            catch
                rsq = nan;
                params = [nan, nan, nan, nan];
            end
        end
        
        
        
        function y = predict(obj, x_vals)
        %%predict Predict firing rate at given x values using fitted parameters
        %
        %   Y = predict(X_VALS) returns the predicted firing rates at
        %   the x values in X_VALS using the fitted parameters.
        
            params = [obj.R_max, obj.x_max, obj.sigma_minus, obj.sigma_plus];
            y = obj.asymmetric_gaussian(params, x_vals);
        end
        
        
        
        function val = get_summary(obj)
        %%get_summary Return the essential properties in a structure to save
        %
        %   STRUCT = get_summary()
        
            props = {'n_reps', ...
                     'R_max', ...
                     'x_max', ...
                     'sigma_minus', ...
                     'sigma_plus', ...
                     'rsq', ...
                     'rmse', ...
                     'rsq_shuff', ...
                     'params_shuff', ...
                     'shuff_tuning', ...
                     'p'};
                 
            for ii = 1 : length(props)
                val.(props{ii}) = obj.(props{ii});
            end
        end
    end
    
    
    
    methods (Static = true)
        
        function y = asymmetric_gaussian(params, x)
        %%asymmetric_gaussian The asymmetric Gaussian function
        %
        %   Y = asymmetric_gaussian(PARAMS, X)
        %   
        %   PARAMS is a 4-element vector: [R_max, x_max, sigma_minus, sigma_plus]
        %   X is the stimulus parameter (temporal frequency)
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
        
        
        
        function [rsq, params] = fit_shuffled_data_static(mean_tuning, x, n_bins)
        %%fit_shuffled_data_static Static version of fit_shuffled_data for parfor
        %
        %   [RSQ, PARAMS] = fit_shuffled_data_static(MEAN_TUNING, X, N_BINS)
        %   Static method that can be called from parfor loop
        
            % Check if we have valid data
            if all(isnan(mean_tuning)) || sum(~isnan(mean_tuning)) < 4
                rsq = nan;
                params = [nan, nan, nan, nan];
                return;
            end
            
            % Initial parameter guesses
            [R_max_init, idx_max] = max(mean_tuning);
            x_max_init = x(idx_max);
            sigma_init = abs(x(end) - x(1)) / 4;
            
            % Parameter bounds
            lb = [0, min(x), 0.001, 0.001];
            ub = [max(mean_tuning) * 2, max(x), abs(x(end) - x(1)) * 2, abs(x(end) - x(1)) * 2];
            x0 = [R_max_init, x_max_init, sigma_init, sigma_init];
            
            % Fit using lsqcurvefit
            options = optimoptions('lsqcurvefit', 'Display', 'off');
            
            try
                [params, ~, ~, ~, ~] = lsqcurvefit(@AsymmetricGaussianFit.asymmetric_gaussian, x0, x, mean_tuning, lb, ub, options);
                
                % Calculate R^2
                y_fit = AsymmetricGaussianFit.asymmetric_gaussian(params, x);
                SS_res = sum((mean_tuning - y_fit).^2);
                SS_tot = sum((mean_tuning - nanmean(mean_tuning)).^2);
                rsq = 1 - SS_res / SS_tot;
                
            catch
                rsq = nan;
                params = [nan, nan, nan, nan];
            end
        end
    end
end
