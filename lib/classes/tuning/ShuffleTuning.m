classdef ShuffleTuning < handle
    
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
        % Shuffle
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
