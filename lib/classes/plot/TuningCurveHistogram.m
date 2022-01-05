classdef TuningCurveHistogram < RC2Axis
% TuningCurveHistogram Class for plotting the distribution of r values from
% the shuffling of the tuning curve data.
%
%  TuningCurveHistogram Properties:
%       dot_col         - color of the dot for the true r value (default = [1, 0, 0])
%       dot_size        - size of the dot for the true r value (default = 10 points)
%       plot_labels     - whether to plot the axis lable (default = true)
%       h_bar           - handle to the bars in the histogram
%
%  TuningCurveHistogram Methods:
%       plot            - plot the shuffled r values

    properties
        
        h_bar
        
        dot_col = [1, 0, 0]
        dot_size = 10
        plot_labels = true;
    end
    
    
    methods
        
        function obj = TuningCurveHistogram(h_ax)
        %%TuningCurveHistogram
        %
        %   TuningCurveHistogram(AXIS_HANDLE)
        
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        function plot(obj, tuning)
        %%plot Plots the distribution of shuffled r values
        %
        %   plot(TUNING) takes a MATLAB structure TUNING with information
        %   about the tuning curves. For more information about the tuning curve
        %   structure see TuningTable.
        %
        %   See also: TuningTable
        
            [n, c] = histcounts(tuning.shuffled.r_shuff, 50);
            obj.h_bar = barh(obj.h_ax, (c(1:end-1)+c(2:end))/2, n, 'facecolor', [0.6, 0.6, 0.6], 'edgecolor', 'none', 'barwidth', 1);
            set(obj.h_ax, 'plotboxaspectratio', [1, 3, 1])
            ylim([-0.5, 0.5]);
            scatter(obj.h_ax, 0, tuning.shuffled.r, obj.dot_size, obj.dot_col, 'fill');
            if obj.plot_labels
                xlabel('count');
                ylabel('r');
            end
        end
    end
end