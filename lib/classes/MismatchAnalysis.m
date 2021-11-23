classdef MismatchAnalysis < handle
    
    properties
        
        method = 'anova'
        n_windows = 4;
        window_t = 0.1;
        
        simple_bsl_window_t = 0.4
    end
    
    
    methods
        
        function obj = MismatchAnalysis()
        %%class of methods aiding analysis of the mismatch responses
        end
        
        
        function val = get_avg_baseline_fr(obj, cluster, trials)
        %%get average firing rate during the baseline period
            if strcmp(obj.method, 'anova')
                mtx = obj.get_baseline_fr_mtx(cluster, trials);
                val = mean(mtx(:));
            elseif strcmp(obj.method, 'ranksum')
                fr = obj.get_fr_simple_baseline(cluster, trials);
                val = median(fr);
            end
        end
        
        
        
        function val = get_avg_response_fr(obj, cluster, trials)
        %%get average firing rate during the response period
            if strcmp(obj.method, 'anova')
                mtx = obj.get_response_fr_mtx(cluster, trials);
                val = mean(mtx(:));
            elseif strcmp(obj.method, 'ranksum')
                fr = obj.get_fr_simple_response(cluster, trials);
                val = median(fr);
            end
        end
        
        
        
        function mtx = get_baseline_fr_mtx(obj, cluster, trials)
        %%get matrix of firing rates during the baseline period 
            offset = -obj.n_windows * obj.window_t;
            mtx = obj.get_fr_mtx(cluster, trials, offset);
        end
        
        
        
        function mtx = get_response_fr_mtx(obj, cluster, trials)
        %%get matrix of firing rates during the response period
            offset = 0;
            mtx = obj.get_fr_mtx(cluster, trials, offset);
        end
        
        
        
        function mtx = get_control_response_fr_mtx(obj, cluster, trials)
        %%get matrix of firing rates during the control period
            offset = -2 * obj.n_windows * obj.window_t;
            mtx = obj.get_fr_mtx(cluster, trials, offset);
        end
        
        
        
        function mtx = get_fr_mtx(obj, cluster, trials, offset)
        %%gets matrix of firing rates for the cluster on the trials
        %%provided.
        %   response matrix gives firing rates during each window with
        %   offset
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
            
            fr = nan(length(trials), 1);
            
            for ii = 1 : length(trials)
                
                mm_onset = trials{ii}.mismatch_onset_t();
                
                T = linspace(-obj.simple_bsl_window_t, 0, round(obj.simple_bsl_window_t*10e3)+1);
                T = mm_onset + T;
                
                fr(ii, 1) = mean(cluster.fr.get_convolution(T));
            end
        end
        
        
        
        function fr = get_fr_simple_response(obj, cluster, trials)
            
            fr = nan(length(trials), 1);
            
            for ii = 1 : length(trials)
                
                mm_onset = trials{ii}.mismatch_onset_t();
                
                T = linspace(0, obj.simple_bsl_window_t, round(obj.simple_bsl_window_t*10e3)+1);
                T = mm_onset + T;
                
                fr(ii, 1) = mean(cluster.fr.get_convolution(T));
            end
        end
        
        
        
        function [sig, p_response, direction] = is_response_significant(obj, cluster, trials)
            %%determines whether the cluster significantly responds to the
            %%mismatch based on the responses to the trials in 'trials'
            
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
            
            baseline_fr = obj.get_fr_simple_baseline(cluster, trials);
            sig = ~lillietest(baseline_fr);
        end
        
        
        function sig = is_response_normal(obj, cluster, trials)
            
            response_fr = obj.get_fr_simple_response(cluster, trials);
            sig = ~lillietest(response_fr);
        end
    end
    
    
    
    methods (Static = true)
        
        function p = anova(g1, g2)
            
            group_1 = [ones(numel(g1), 1); 2*ones(numel(g2), 1)];
            sg_1 = repmat((1:size(g1, 1))', 1, size(g1, 2));
            sg_2 = repmat((1:size(g2, 1))', 1, size(g2, 2));
            group_2 = [sg_1(:); sg_2(:)];
            
            p = anovan([g1(:); g2(:)], {group_1, group_2}, 'display', 'off');
            p = p(1);
        end
    end
end
