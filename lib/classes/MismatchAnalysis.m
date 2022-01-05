classdef MismatchAnalysis < handle
% MismatchAnalysis Class for handling analysis of mismatch events
%
%   MismatchAnalysis Properties:
%       method              - 'anova', not currently used
%       n_windows           - number of time windows to split the baseline and response periods into
%       window_t            - size of each window
%       simple_bsl_window_t - unused
%
%   MismatchAnalysis Methods:
%       get_avg_baseline_fr            - get the average firing rate of a cluster in a baseline period for a set of trials
%       get_avg_response_fr            - get the average firing rate of a cluster in a response period for a set of trials
%       get_baseline_fr_mtx            - get the firing rates of a cluster in the windows of the baseline period
%       get_response_fr_mtx            - get the firing rates of a cluster in the windows of the response period
%       get_control_response_fr_mtx    - get the firing rates of a cluster in the windows of the control period
%       get_fr_mtx                     - shared function for getting the firing rate matrices
%       get_fr_simple_baseline         - get the firing rates of a cluster in the baseline periods
%       get_fr_simple_response         - get the firing rates of a cluster in the response periods
%       is_response_significant        - determine whether the response to the mismatch is significant (increase or decrease)
%       is_baseline_normal             - tests whether the trial-to-trial firing rates in the baseline can 
%                                        be considered normally distributed
%       is_response_normal             - tests whether the trial-to-trial firing rates in the response can 
%                                        be considered normally distributed
%       anova                          - perfom a 2-way anova on the firing rates


    properties
        
        method = 'anova'
        n_windows = 4;
        window_t = 0.1;
        
        simple_bsl_window_t = 0.4
    end
    
    
    methods
        
        function obj = MismatchAnalysis()
        %%MismatchAnalysis
        %
        %   MismatchAnalysis()
        %   class of methods handling analysis of mismatch events
        
        end
        
        
        function val = get_avg_baseline_fr(obj, cluster, trials)
        %%get_avg_baseline_fr Get the average firing rate of a cluster in a baseline
        %%period for a set of trials
        %
        %   FIRING_RATE = get_avg_baseline_fr(CLUSTER, TRIALS) for a cell
        %   array of (mismatch) Trial objects, TRIALS, get the average
        %   firing rate in the mismatch baseline window for the cluster
        %   CLUSTER (object of class Cluster).
        %
        %   A single value is returned as FIRING_RATE
        
            if strcmp(obj.method, 'anova')
                mtx = obj.get_baseline_fr_mtx(cluster, trials);
                val = mean(mtx(:));
            elseif strcmp(obj.method, 'ranksum')
                fr = obj.get_fr_simple_baseline(cluster, trials);
                val = median(fr);
            end
        end
        
        
        
        function val = get_avg_response_fr(obj, cluster, trials)
        %%get_avg_response_fr Get the average firing rate of a cluster in a
        %%response period for a set of trials
        %
        %   FIRING_RATE = get_avg_response_fr(CLUSTER, TRIALS) for a cell
        %   array of (mismatch) Trial objects, TRIALS, get the average
        %   firing rate in the mismatch response window for the cluster
        %   CLUSTER (object of class Cluster).
        
            if strcmp(obj.method, 'anova')
                mtx = obj.get_response_fr_mtx(cluster, trials);
                val = mean(mtx(:));
            elseif strcmp(obj.method, 'ranksum')
                fr = obj.get_fr_simple_response(cluster, trials);
                val = median(fr);
            end
        end
        
        
        
        function mtx = get_baseline_fr_mtx(obj, cluster, trials)
        %%get_baseline_fr_mtx Get the firing rates of a cluster in the
        %%windows of the baseline period
        %
        %   FIRING_RATE_MTX = get_baseline_fr_mtx(CLUSTER, TRIALS) for a cell
        %   array of (mismatch) Trial objects, TRIALS, get the firing rates
        %   in the mismatch baseline windows for the cluster 
        %   CLUSTER (object of class Cluster). Firing rates are returned as
        %   a matrix #trials x #windows
        
            offset = -obj.n_windows * obj.window_t;
            mtx = obj.get_fr_mtx(cluster, trials, offset);
        end
        
        
        
        function mtx = get_response_fr_mtx(obj, cluster, trials)
        %%get_response_fr_mtx Get the firing rates of a cluster in the
        %%windows of the response period
        %
        %   FIRING_RATE_MTX = get_response_fr_mtx(CLUSTER, TRIALS) for a cell
        %   array of (mismatch) Trial objects, TRIALS, get the firing rates
        %   in the mismatch response windows for the cluster 
        %   CLUSTER (object of class Cluster). Firing rates are returned as
        %   a matrix #trials x #windows
        
            offset = 0;
            mtx = obj.get_fr_mtx(cluster, trials, offset);
        end
        
        
        
        function mtx = get_control_response_fr_mtx(obj, cluster, trials)
        %%get_control_response_fr_mtx Get the firing rates of a cluster in the
        %%windows of the control period
        %
        %   FIRING_RATE_MTX = get_control_response_fr_mtx(CLUSTER, TRIALS)
        %   for a cell array of (mismatch) Trial objects, TRIALS, get the
        %   firing rates in the mismatch control windows for the cluster  
        %   CLUSTER (object of class Cluster). Firing rates are returned as
        %   a matrix #trials x #windows. The control period comes before
        %   the baseline period.
        
            offset = -2 * obj.n_windows * obj.window_t;
            mtx = obj.get_fr_mtx(cluster, trials, offset);
        end
        
        
        
        function mtx = get_fr_mtx(obj, cluster, trials, offset)
        %%get_fr_mtx Shared function for getting the firing rate matrices
        %
        %   FIRING_RATE_MTX = get_fr_mtx(CLUSTER, TRIALS, OFFSET) for a cell
        %   array of (mismatch) Trial objects, TRIALS, get the firing rates
        %   in the windows starting at offset, OFFSET (specified in
        %   seconds) for the cluster CLUSTER (object of class Cluster).
        %   Firing rates are returned as  a matrix #trials x #windows.
        
            mtx = nan(length(trials), obj.n_windows);
            
            for ii = 1 : length(trials)
                
                mm_onset = trials{ii}.mismatch_onset_t();
                
                limits = mm_onset + offset + [(0:obj.n_windows-1)', (1:obj.n_windows)']*obj.window_t;
                
                for jj = 1 : obj.n_windows
                    
                    mtx(ii, jj) = cluster.fr.get_fr_in_window(limits(jj, :));
                end
            end
        end
        
        
        
        function fr = get_fr_simple_baseline(obj, cluster, trials)
        %%get_fr_simple_baseline Get the firing rates of a cluster in the
        %%baseline periods
        %
        %   FIRING_RATE = get_fr_simple_baseline(CLUSTER, TRIALS) for a cell
        %   array of Trial objects, TRIALS, get the firing rates
        %   in the baseline window (-`simple_bsl_window_t` to 0 before
        %   mismatch onset) for the cluster CLUSTER (object of class
        %   Cluster).  Firing rates are returned as  a vector #trials x 1.
        
            fr = nan(length(trials), 1);
            
            for ii = 1 : length(trials)
                
                mm_onset = trials{ii}.mismatch_onset_t();
                
                T = linspace(-obj.simple_bsl_window_t, 0, round(obj.simple_bsl_window_t*10e3)+1);
                T = mm_onset + T;
                
                fr(ii, 1) = mean(cluster.fr.get_convolution(T));
            end
        end
        
        
        
        function fr = get_fr_simple_response(obj, cluster, trials)
        %%get_fr_simple_response Get the firing rates of a cluster in the
        %%response periods
        %
        %   FIRING_RATE = get_fr_simple_response(CLUSTER, TRIALS) for a cell
        %   array of Trial objects, TRIALS, get the firing rates
        %   in the response window (0 to `simple_bsl_window_t` after
        %   mismatch onset) for the cluster CLUSTER (object of class
        %   Cluster).  Firing rates are returned as  a vector #trials x 1. 
        
            fr = nan(length(trials), 1);
            
            for ii = 1 : length(trials)
                
                mm_onset = trials{ii}.mismatch_onset_t();
                
                T = linspace(0, obj.simple_bsl_window_t, round(obj.simple_bsl_window_t*10e3)+1);
                T = mm_onset + T;
                
                fr(ii, 1) = mean(cluster.fr.get_convolution(T));
            end
        end
        
        
        
        function [sig, p_response, direction] = is_response_significant(obj, cluster, trials)
        %%is_response_significant Determine whether the response to the
        %%mismatch is significant (increase or decrease)
        %
        %   [SIGNIFICANT, P_VALUE, DIRECTION] = is_response_significant(CLUSTER, TRIALS)
        %   for a cell array of (mismatch) Trial objects, TRIALS, determine whether
        %   the response of cluster CLUSTER (object of class Cluster) is
        %   significant.
        %
        %   Outputs:
        %    SIGNIFICANT    - whether the test is significant 1 or not 0.
        %                     significance is determined by the P_VALUE < 0.05
        %    P_VALUE        - p-value of the TEST
        %    DIRECTION      - either 1, 0, or -1
        %                     1  - mismatch response > baseline firing
        %                     0  - no change in firing rate
        %                     -1 - mismatch response < baseline firing
        %
        %   Currently, the test performed is determined by the property
        %   `method` which is either 'anova' or 'ranksum'. If 'anova', the
        %   time around the mismatch onset is split into multiple windows
        %   (`n_windows` property), each of size `window_t` seconds both
        %   before (baseline), after (response) and before the baseline
        %   (control). That is, if `n_windows` is 4, `window_t` is 0.1,
        %   there will be a total of 12 windows, spanning 1.2 seconds
        %   around the mismatch onset from -0.8 seconds before to 0.4
        %   seconds after. The 4 windows between -0.8 and -0.4 are part of
        %   the control, the 4 windows between -0.4 and 0 are part of the
        %   baseline and the 4 windows between 0 and 0.4 are part of the
        %   response.
        %
        %   The firing rate is determined in each of these windows for each
        %   trials, giving 3 matrices of firing (#trials x #windows). A two
        %   way ANOVA is used to compare control to baseline and baseline
        %   to response. If the main effect is significant between the
        %   control and baseline, the p-value is set to NaN and the cluster
        %   is not deemed to be significant. If the main effect is not
        %   significant between control and baseline, the cluster is
        %   significant if the p-value of the main effect between baseline
        %   and response is < 0.05.
        %
        %   See also: FormattedData.is_mismatch_significant
            
            if strcmp(obj.method, 'anova')
                baseline_mtx    = obj.get_baseline_fr_mtx(cluster, trials);
                response_mtx    = obj.get_response_fr_mtx(cluster, trials);
                control_mtx     = obj.get_control_response_fr_mtx(cluster, trials);
                
                p_response      = obj.anova(baseline_mtx', response_mtx');
                p_control       = obj.anova(baseline_mtx', control_mtx');
                
                if p_control < 0.05
                    p_response = nan;
                    sig = false;
                else
                    sig = p_response < 0.05;
                end
                
                if sig
                    if mean(baseline_mtx(:)) < mean(response_mtx(:))
                        direction = 1;
                    else
                        direction = -1;
                    end
                else
                    direction = 0;
                end
                
            elseif strcmp(obj.method, 'ranksum')
                
                baseline_fr = obj.get_fr_simple_baseline(cluster, trials);
                response_fr = obj.get_fr_simple_response(cluster, trials);
                
                p_response = signrank(baseline_fr, response_fr);
                
                sig = p_response < 0.05;
                
                if sig
                    if median(baseline_fr) < median(response_fr)
                        direction = 1;
                    elseif median(baseline_fr) > median(response_fr)
                        direction = -1;
                    elseif median(baseline_fr) == median(response_fr)
                        error('medians are the same')
                    end
                else
                    direction = 0;
                end
            end
        end
        
        
        function sig = is_baseline_normal(obj, cluster, trials)
        %%is_baseline_normal   Tests whether the trial-to-trial firing
        %%rates in the baseline can be considered normally distributed
        %
        %   SIGNIFICANT = is_baseline_normal(CLUSTER, TRIALS)
        %   computes the firing rates in the baseline periods for mismatch
        %   trials in the cell array of Trial objects TRIALS, for cluster
        %   CLUSTER (object of class Cluster). Returns SIGNIFICANT, which
        %   is false if the firing rates reject the null hypothesis of a
        %   Lilliefors test, and true otherwise.
        
            baseline_fr = obj.get_fr_simple_baseline(cluster, trials);
            sig = ~lillietest(baseline_fr);
        end
        
        
        function sig = is_response_normal(obj, cluster, trials)
        %%is_response_normal   Tests whether the trial-to-trial firing
        %%rates in the baseline can be considered normally distributed
        %
        %   SIGNIFICANT = is_response_normal(CLUSTER, TRIALS)
        %   computes the firing rates in the response periods for mismatch
        %   trials in the cell array of Trial objects TRIALS, for cluster
        %   CLUSTER (object of class Cluster). Returns SIGNIFICANT, which
        %   is false if the firing rates reject the null hypothesis of a
        %   Lilliefors test, and true otherwise.   
        
            response_fr = obj.get_fr_simple_response(cluster, trials);
            sig = ~lillietest(response_fr);
        end
    end
    
    
    
    methods (Static = true)
        
        function p = anova(g1, g2)
        %%anova Perfom a 2-way anova on the firing rates
        %
        %   P_VALUE = anova(FIRING_RATES_1, FIRING_RATES_2)
        %   takes two matrices of firing rates of shape #trialsx#windows
        %   and looks to see if there is a main effect between the two
        %   groups.
            
            group_1 = [ones(numel(g1), 1); 2*ones(numel(g2), 1)];
            sg_1 = repmat((1:size(g1, 1))', 1, size(g1, 2));
            sg_2 = repmat((1:size(g2, 1))', 1, size(g2, 2));
            group_2 = [sg_1(:); sg_2(:)];
            
            p = anovan([g1(:); g2(:)], {group_1, group_2}, 'display', 'off');
            p = p(1);
        end
    end
end
