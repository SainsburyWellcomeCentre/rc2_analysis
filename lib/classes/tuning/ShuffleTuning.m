classdef ShuffleTuning < handle
% ShuffleTuning Class for getting tuning information
%
%   ShuffleTuning Properties:
%         n_reps        - number of shuffles to perform
%         x             - # bins x 1 vector giving velocity of bin centers
%         tuning        - # bins x # trials giving firing rate in a velocity bin for each trial
%         n_bins        - # bins
%         n_trials      - # trials
%         rsq           - the R^2 of a linear fit to the data
%         beta          - the parameters of the linear fit
%         r             - the r value of the fit
%         
%         rsq_shuff     - n_reps x 1 vector, R^2 values for shuffled data
%         beta_shuff    - n_reps x 2 vector, parameters of each shuffled fit
%         r_shuff       - n_reps x 1 vector, correlation 
%         shuff_tuning  - # bins x n_reps containing average for each shuffle
%         shuff_sd      - # bins x n_reps containing SD for each shuffle
%         shuff_n       - # bins x n_reps containing n trials averaged for each shuffle
%         
%         p             - p-value from the shuffling procedure
%         p_lm          - p-value from the linear fit
%
%   ShuffleTuning Methods:
%         get_shuffled_rsq - perform the shuffling and get R^2 for each shuffle
%         get_summary      - return properties of this class as a structure
%         get_rsq          - perform the fitting
%
%   See also: TuningTable, VelocityBins, VelocityTuningCurve

    properties
        
        n_reps = 1000
        x
        tuning
        n_bins
        n_trials
        
        rsq
        beta
        r
        
        rsq_shuff
        beta_shuff
        r_shuff
        shuff_tuning
        shuff_sd
        shuff_n
        
        p
        p_lm
    end
    
    
    
    methods
        
        function obj = ShuffleTuning(tuning, x)
        % ShuffleTuning
        %
        %   ShuffleTuning(FIRING_RATE, VELOCITY)
        %   creates an object for fitting tuning cuves and shuffling the
        %   data to determine significance of tuning.
        %
        %   FIRING_RATE is a # velocity bins x # trials matrix with the
        %   firing rates of a cluster to a velocity bin for a set of
        %   trials. VELOCITY is a #velocity bins x 1 vector with the centers 
        %   of the velocity bins in which the firing rates were computed.
        %
        %   A linear fit is performed on the data. Then the FIRING_RATE
        %   matrix is completely shuffled, and the fitting performed
        %   multiple times.
        %   
        %   A p-value is calculated by using the fraction of shuffled R^2 value
        %   above the true R^2 value.
        
            obj.x = x(:); % n bins x 1
            obj.tuning = tuning; % n bins x n trials
            obj.n_bins = size(obj.tuning, 1);
            obj.n_trials = size(obj.tuning, 2);
            
            % fit to data R^2
            x_ = repmat(obj.x, 1, obj.n_trials);
            lm = fitlm(x_(:), tuning(:));
            obj.p_lm = lm.Coefficients.pValue(2);
            
            [obj.rsq, obj.beta, obj.r] = obj.get_rsq(x_(:), tuning(:));
            
            obj.get_shuffled_rsq();
            
            if isnan(obj.rsq)
                obj.p = nan;
            else
                p_up = sum(obj.rsq_shuff > obj.rsq)/obj.n_reps;
%                 p_down = sum(obj.rsq_shuff < obj.rsq)/obj.n_reps;
                obj.p = p_up; %min(p_up, p_down);
            end
        end
        
        
        
        function get_shuffled_rsq(obj)
        %%get_shuffed_rsq Perform the shuffling and get R^2 for each
        %%shuffle
        %
        %   get_shuffled_rsq() takes the data in `tuning` and shuffles
        %   with replacement `n_reps` times. Stores the results in
        %   properties of the object.
        
            rng(1);
            
            x_ = repmat(obj.x, 1, obj.n_trials);
            
            % preallocate arrays
            obj.rsq_shuff = nan(1, obj.n_reps);
            obj.beta_shuff = nan(obj.n_reps, 2);
            obj.shuff_tuning = nan(obj.n_bins, obj.n_reps);
            obj.shuff_sd = nan(obj.n_bins, obj.n_reps);
            obj.shuff_n = nan(obj.n_bins, obj.n_reps);
            
            for rand_i = 1 : obj.n_reps
                I = randi(numel(obj.tuning), size(obj.tuning));
                new_tuning = obj.tuning(I);
                
                [obj.rsq_shuff(rand_i), obj.beta_shuff(rand_i, :), obj.r_shuff(rand_i, :)] = ...
                    obj.get_rsq(x_(:), new_tuning(:));
                
                obj.shuff_tuning(:, rand_i) = nanmean(new_tuning, 2);
                obj.shuff_sd(:, rand_i) = nanstd(new_tuning, [], 2);
                obj.shuff_n(:, rand_i) = sum(~isnan(new_tuning), 2);
            end
        end
        
        
        function val = get_summary(obj)
        %%get_summary Return the essential properties in a structure to save
        %
        %   STRUCT = get_summary()
        
            props = {'n_reps', ...
                     'rsq', ...
                     'beta', ...
                     'r', ...
                     'rsq_shuff', ...
                     'beta_shuff', ...
                     'r_shuff', ...
                     'shuff_tuning', ...
                     'p', ...
                     'p_lm'};
                 
            for ii = 1 : length(props)
                val.(props{ii}) = obj.(props{ii});
            end
        end
    end
    
    methods (Static = true)
        
        function [rsq, beta, r] = get_rsq(x, tuning)
        %%get_rsq Performs a linear fit to the tuning matrix
        %
        %   [R_sq, BETA, R] = get_rsp(VELOCITY, FIRING_RATE)
        %   takes a # bins x # trials FIRING_RATE matrix and # biins x #
        %   trials VELOCITY matrix (velocities at which FIRING_RATE
        %   measured) and performs linear fit and correlation.
        %
        %   R_sq & BETA are the R^2 and parameters of the lienar fit
        %   and R is the correlation between VELOCITY and FIRING_RATE.
        
            x(isnan(tuning)) = [];
            tuning(isnan(tuning)) = [];
            
            beta = polyfit(x, tuning, 1);
            yfit =  beta(1)*x + beta(2);
            yresid = tuning - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(tuning)-1) * var(tuning);
            rsq = 1 - SSresid/SStotal;
            r = corr(x, tuning);
        end
    end
end
