classdef TFTuningCurvePlot < RC2Axis
% TFTuningCurvePlot Class for plotting temporal frequency tuning curves with model fits
%
%   Supports multiple fit types: linear, Gaussian, asymmetric Gaussian, and sigmoid.
%   Automatically detects the model type from the tuning data structure.
%
%   TFTuningCurvePlot Properties:
%       error_type      - 'sem' (default) or 'sd' determines the error bar to plot
%       print_stats     - true (default) or false, whether to print stats on figure
%       dot_size        - size of dots for line plot (default = 10 points)
%       line_width      - thickness of main line (default = 1 point)
%       n_shuffs_to_plot - number of shuffled fits to plot in background
%       main_col        - colour of the true fit and tuning curve
%
%   TFTuningCurvePlot Methods:
%       plot            - plot the tuning curve and stats
%       color           - change the colour of the tuning curve
%       ylim            - get or set the y-axis limits
%       xlim            - get or set the x-axis limits

    properties
        
        error_type = 'sem'
        print_stats = true
        dot_size = 10
        line_width = 1
        
        n_shuffs_to_plot = 4
        main_col = [0, 0, 0]
    end

    properties (SetAccess = private)
        
        h_line
        h_dots
        h_errorbars
        h_fit
        h_txt
        
        h_line_shuff
        h_dots_shuff
        h_errorbars_shuff
        h_fit_shuff
        
        shuff
        x_fit  % Store x values for fitting curve
    end
    
    
    
    methods
        
        function obj = TFTuningCurvePlot(h_ax)
        %%TFTuningCurvePlot
        %
        %   TFTuningCurvePlot(AXIS_HANDLE)
        
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        function plot(obj, tuning)
        %%plot Plots the tuning curve
        %
        %   plot(TUNING) plots the tuning curve information contained in
        %   the structure TUNING. For more information about the tuning curve
        %   structure see TFTuningTable. The information plotted includes:
        %   the mean firing rate in each TF bin, along with the fitted model.
        %   Several fits to the shuffled data in light grey.
        %   Further, statistics about the tuning are printed.
        %
        %   See also: TFTuningTable, ModelSelectionTuning
            
            obj.shuff = tuning.shuffled;
            fr      = nanmean(tuning.tuning, 2);
            sd      = nanstd(tuning.tuning, [], 2);
            n       = sum(~isnan(tuning.tuning), 2);
            x       = tuning.bin_centers;
            
            % Create fine-grained x values for smooth curves
            obj.x_fit = linspace(min(x), max(x), 100);
        
            % Note: Shuffled fit parameters (beta) are not stored, only metrics (R², BIC).
            % Therefore we cannot plot the shuffled fitted curves. The shuffled distribution
            % is shown in the histogram of R² values instead.
            
            stat_fr = nanmean(tuning.stationary_fr);
            stat_sd = nanstd(tuning.stationary_fr);
            stat_n  = sum(~isnan(tuning.stationary_fr));
            
            % Plot data points with error bars
            obj.h_dots = scatter(obj.h_ax, x, fr, obj.dot_size, obj.main_col, 'fill');
            for i = 1 : length(fr)
                if strcmp(obj.error_type, 'std')
                    y = fr(i) + sd(i) * [-1, 1];
                else
                    y = fr(i) + (sd(i)/sqrt(n(i))) * [-1, 1];
                end
                
                obj.h_errorbars(i) = line(obj.h_ax, x([i, i]), y, 'color', obj.main_col, 'linewidth', obj.line_width);
            end
            
            % Plot the true fit (model-dependent)
            y_fit = obj.evaluate_model(tuning.shuffled, obj.x_fit);
            if ~isempty(y_fit)
                obj.h_fit = line(obj.h_ax, obj.x_fit, y_fit, 'color', obj.main_col, 'linewidth', obj.line_width);
            end
            
            % Plot stationary firing rate at x=0
            scatter(obj.h_ax, 0, stat_fr, obj.dot_size, [0, 0, 0], 'fill')
            
            if obj.print_stats
                
                % Build stats string based on model type
                str = obj.build_stats_string(tuning.shuffled);
                
                if tuning.shuffled.p < 0.05
                    col = 'r';
                else
                    col = 'b';
                end
                
                obj.h_txt = text(obj.h_ax, obj.xmax, obj.ymin, ...
                    str, 'color', col, ...
                    'verticalalignment', 'bottom', 'horizontalalignment', 'right', ...
                    'fontsize', 6);
            end
            
            set(obj.h_ax, 'xlim', [obj.xmin, obj.xmax]);
        end
        
        
        
        function color(obj, val)
        %%color Sets the colour of the true fit and tuning curve
        %
        %   color(VALUE) where VALUE is a MATLAB colour format.
        
            set(obj.h_line, 'color', val);
            for i = 1 : length(obj.h_errorbars)
                set(obj.h_errorbars(i), 'color', val);
            end
        end
        
        
        
        function val = ylim(obj, val)
        %%ylim Get or set the y-axis limits
        %
        %   VALUEOUT = ylim(VALUEIN) similar to RC2Axis.ylim, but also
        %   adjusts the position of the text information.
        %
        %   See also: RC2Axis
        
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'ylim');
            if obj.print_stats
                pos = get(obj.h_txt, 'position');
                set(obj.h_txt, 'position', [pos(1), val(1), 0]);
            end
        end
        
        
        
        function val = xlim(obj, val)
        %%xlim Get or set the x-axis limits
        %
        %   VALUEOUT = xlim(VALUEIN) similar to RC2Axis.xlim, but also
        %   adjusts the fitted curves and text information.
        %
        %   See also: RC2Axis
        
            VariableDefault('val', []);
            val = obj.get_set_limits(val, 'xlim');
            
            if obj.print_stats
                pos = get(obj.h_txt, 'position');
                set(obj.h_txt, 'position', [val(2), pos(2), 0]);
            end
            
            % Update x range for curves
            obj.x_fit = linspace(val(1), val(2), 100);
            
            % Update true fit curve
            y_fit = obj.evaluate_model(obj.shuff, obj.x_fit);
            if ~isempty(y_fit) && ~isempty(obj.h_fit)
                set(obj.h_fit, 'xdata', obj.x_fit, 'ydata', y_fit);
            end
            
            % Note: Shuffled curves are bin averages, not re-evaluated models
        end
        
        
        
        function y = evaluate_model(obj, shuffled, x)
        %%evaluate_model Evaluate the fitted model at given x values
        %
        %   Y = evaluate_model(SHUFFLED, X) evaluates the model at x points
        
            model_type = shuffled.best_model;
            params = shuffled.best_model_info.beta;
            
            % Evaluate model
            switch model_type
                case 'linear'
                    y = polyval(params, x);
                case 'quadratic'
                    y = polyval(params, x);
                case 'cubic'
                    y = polyval(params, x);
                case 'gaussian'
                    % params = [amplitude, mu, sigma, baseline]
                    y = params(1) * exp(-(x - params(2)).^2 / (2*params(3)^2)) + params(4);
                case 'asymmetric_gaussian'
                    % params = [R_max, x_max, sigma_minus, sigma_plus]
                    y = ModelSelectionTuning.asymmetric_gaussian_static(params, x);
                case 'sigmoid'
                    % params = [amplitude, k, x0, baseline]
                    y = params(1) ./ (1 + exp(-params(2) * (x - params(3)))) + params(4);
            end
        end
        
        
        
        function str = build_stats_string(obj, shuffled)
        %%build_stats_string Build statistics string based on model type
        %
        %   STR = build_stats_string(SHUFFLED)
        
            model_type = shuffled.best_model;
            params = shuffled.best_model_info.beta;
            
            % Build string based on model
            str = sprintf('Model: %s\n', model_type);
            
            % Show R² (used for significance testing)
            rsq = shuffled.best_model_info.rsq;
            str = [str, sprintf('R^2 = %.2f\n', rsq)];
            
            switch model_type
                case 'linear'
                    str = [str, sprintf('slope = %.2f\n', params(1))];
                    str = [str, sprintf('intercept = %.2f\n', params(2))];
                case 'quadratic'
                    str = [str, sprintf('a = %.2e\n', params(1))];
                    str = [str, sprintf('b = %.2f\n', params(2))];
                    str = [str, sprintf('c = %.2f\n', params(3))];
                case 'cubic'
                    str = [str, sprintf('a = %.2e\n', params(1))];
                    str = [str, sprintf('b = %.2e\n', params(2))];
                    str = [str, sprintf('c = %.2f\n', params(3))];
                    str = [str, sprintf('d = %.2f\n', params(4))];
                case 'gaussian'
                    str = [str, sprintf('A = %.2f Hz\n', params(1))];
                    str = [str, sprintf('\\mu = %.2f Hz\n', params(2))];
                    str = [str, sprintf('\\sigma = %.2f\n', params(3))];
                case 'asymmetric_gaussian'
                    str = [str, sprintf('R_{max} = %.2f Hz\n', params(1))];
                    str = [str, sprintf('x_{max} = %.2f Hz\n', params(2))];
                    str = [str, sprintf('\\sigma_- = %.2f\n', params(3))];
                    str = [str, sprintf('\\sigma_+ = %.2f\n', params(4))];
                case 'sigmoid'
                    str = [str, sprintf('A = %.2f Hz\n', params(1))];
                    str = [str, sprintf('k = %.2f\n', params(2))];
                    str = [str, sprintf('x_0 = %.2f Hz\n', params(3))];
            end
            
            % Add p-value
            if shuffled.p < 0.05
                str = [str, sprintf('p=%.2e', shuffled.p)];
            else
                str = [str, sprintf('p=%.2f', shuffled.p)];
            end
        end
    end
end
