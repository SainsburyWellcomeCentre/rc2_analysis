classdef MismatchAnalysis < handle
    
    properties
        
        n_windows = 4;
        window_t = 0.1;
    end
    
    
    methods
        
        function obj = MismatchAnalysis()
        %%class of methods aiding analysis of the mismatch responses
        end
        
        
        function val = get_avg_baseline_fr(obj, cluster, trials)
        %%get average firing rate during the baseline period
            mtx = obj.get_baseline_fr_mtx(cluster, trials);
            val = mean(mtx(:));
        end
        
        
        
        function val = get_avg_response_fr(obj, cluster, trials)
        %%get average firing rate during the response period
            mtx = obj.get_response_fr_mtx(cluster, trials);
            val = mean(mtx(:));
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
        
        
        
        function [sig, p_response, direction] = is_response_significant(obj, cluster, trials)
        %%determines whether the cluster significantly responds to the
        %%mismatch based on the responses to the trials in 'trials'
        
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
