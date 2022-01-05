classdef MismatchNov2020Session < RVTSession
% MismatchNov2020Session Class for handling details of a
% mismatch_nov2020 session.
%
%   MismatchNov2020Session Properties:
%       trial_group_ids         - integer IDs for each trial group
%       trial_group_labels      - string IDs for each trial group
%
%   MismatchNov2020Session Methods:
%       get_trial_group_id      - Given a trial, return the trial group ID
%
%   See also: RVTSession  
%
%   TODO: change the name of the experiment group
    
    properties (Constant = true)
        
        trial_group_ids      = 1 : 4
        trial_group_labels    = {'RVT_gain_up', 'RVT_gain_down', 'RV_gain_up', 'RV_gain_down'};
    end
    
    
    
    methods
        
        function obj = MismatchNov2020Session(session)
        % MismatchNov2020Session
        %
        %   MismatchNov2020Session(SESSION) takes a session structure,
        %   SESSION, saved in the formatted data file and creates an object
        %   to handle the data of a 'mismatch_nov2020' session.
        
            obj = obj@RVTSession(session);
            
            for ii = 1 : session.n_trials
                
                % get which group the trial belongs to
                trial_group_id = obj.get_trial_group_id(session.trials(ii));
                
                % create trial object (may have to subclass this at one
                % point)
                obj.trials{ii} = Trial(session.trials(ii), obj);
                
                % set the group id and label
                obj.trials{ii}.trial_group_id = trial_group_id;
                obj.trials{ii}.trial_group_label = obj.get_group_label_from_id(trial_group_id);
            end
        end
        
        
        
        function group_id = get_trial_group_id(~, trial)
        %%get_trial_group_id Given a trial, return the trial group ID
        %
        %   GROUP_ID = get_trial_group_id(TRIAL)
        %   get which trial group ID the given trial in TRIAL, of class
        %   Trial, belongs to.
        
            if strcmp(trial.protocol, 'CoupledMismatch')
                if strcmp(trial.config.gain_direction, 'up')
                    group_id = 1;
                else
                    group_id = 2;
                end
            elseif strcmp(trial.protocol, 'EncoderOnlyMismatch')
                if strcmp(trial.config.gain_direction, 'up')
                    group_id = 3;
                else
                    group_id = 4;
                end
            end
        end
        
        
        
        function trials = get_trials_with_trial_group_label(obj, varargin)
        %%get_trials_with_trial_group_label Return cell array of Trial
        %%objects which are part of a trial group label
        %
        %   TRIALS = get_trials_with_trial_group_label(VARARGIN)
        %   reimplementation of get_trials_with_trial_group_label to remove
        %   some trials which didn't have a mismatch - reason unsolved
        %
        %   See also: RVTSession.get_trials_with_trial_group_label
        %
        %   TODO:   is VARARGIN compatible with the RVTSession
        %           implmentation?
        
            trials = get_trials_with_trial_group_label@RVTSession(obj, varargin{:});
            
            mm_duration = cellfun(@(x)(x.mismatch_duration), trials);
            trials = trials(mm_duration > 0.05);
        end
        
        
        
        
        
        
        function idx = get_svm_table_index(obj, cluster_id, protocol_id)
        %%TODO: UNUSED, REMOVE
        %%Overwrite RVTSession method as it does not work for the mismatch protocol
        %   need better design upstream
        
            trial_type = obj.protocol_type{obj.protocol_ids == protocol_id};
            these_trials = obj.trials_of_type(protocol_id);
            ids_for_these_trials = [these_trials(:).id];
            
            idx = obj.svm_table.cluster_id == cluster_id & ...
                    obj.svm_table.protocol == trial_type & ...
                    ismember(obj.svm_table.trial_id, ids_for_these_trials);
        end
        
        
        
        function [baseline, response, response_ctl] = windowed_mm_responses(obj, cluster, prot_i)
        %%TODO: UNUSED, REMOVE    
            n_windows = 4;
            window_t = 0.1;
            
            cluster_fr = FiringRate(cluster.spike_times);
            trials = obj.trials_of_type(prot_i);
            
            baseline = nan(length(trials), n_windows);
            response = nan(length(trials), n_windows);
            response_ctl = nan(length(trials), n_windows);
            
            for ti = 1 : length(trials)
                
                mm_onset = trials(ti).mismatch_onset_t();
                mm_offset = trials(ti).mismatch_offset_t();
                
                if mm_offset - mm_onset < 0.05
                    continue
                end
                
                rc_lims = mm_onset - 2*n_windows*window_t + [(0:n_windows-1)', (1:n_windows)']*window_t;
                b_lims = mm_onset - n_windows*window_t + [(0:n_windows-1)', (1:n_windows)']*window_t;
                r_lims = mm_onset + [(0:n_windows-1)', (1:n_windows)']*window_t;
                
                for wi = 1 : n_windows
                    
                    response_ctl(ti, wi) = cluster_fr.get_fr_in_window(rc_lims(wi, :));
                    baseline(ti, wi) = cluster_fr.get_fr_in_window(b_lims(wi, :));
                    response(ti, wi) = cluster_fr.get_fr_in_window(r_lims(wi, :));
                end
            end
        end
        
        
        function response_magnitude = average_mismatch_response_by_protocol(obj, cluster, prot_i)
        %%TODO: UNUSED, REMOVE    
            [baseline, response, ~] = obj.windowed_mm_responses(cluster, prot_i);
            
            avg_response = mean(response(:));
            avg_bsl = mean(baseline(:));
            response_magnitude = avg_response - avg_bsl;
        end
        
        
        
        function [running, t] = running_around_mismatch_by_protocol(obj, prot_i, limits)
        %%TODO: UNUSED, REMOVE    
            trials = obj.trials_of_type(prot_i);
            
            n_trials = length(trials);
            n_samples = range(limits) * trials(1).fs;
            
            running = nan(n_samples, n_trials);
            
            for trial_i = 1 : n_trials
                
                mm_onset = trials(trial_i).mismatch_onset_t();
                start_idx = find(trials(trial_i).probe_t > mm_onset + limits(1), 1, 'first');
                full_idx = start_idx + (0:n_samples-1);
                running(:, trial_i) = trials(trial_i).velocity(full_idx);
            end
            
            t = limits(1) + (0:n_samples-1)*(1/trials(1).fs);
        end
        
        
        
        function [spike_rate, t, common_t] = firing_around_mismatch_by_protocol(obj, cluster, prot_i, limits)
        %%TODO: UNUSED, REMOVE    
            cluster_fr = FiringRate(cluster.spike_times);
            trials = obj.trials_of_type(prot_i);
            
            n_trials = length(trials);
            n_samples = range(limits) * trials(1).fs;
            
            spike_rate = nan(n_samples, n_trials);
            
            common_t = limits(1) + (0:n_samples-1)*(1/trials(1).fs);
            
            for trial_i = 1 : n_trials
                
                mm_onset = trials(trial_i).mismatch_onset_t();
                t = common_t + mm_onset;
                
                spike_rate(:, trial_i) = cluster_fr.get_convolution(t);
            end
        end
        
        
        
        function [spike_rate, common_t] = average_firing_around_mismatch_by_protocol(obj, cluster, prot_i, limits)
        %%TODO: UNUSED, REMOVE    
            [spike_rate, ~, common_t] = obj.firing_around_mismatch_by_protocol(cluster, prot_i, limits);
            spike_rate = mean(spike_rate, 2);
        end
        
        
        
        function [running, t] = running_around_mismatch(obj, trial_id, limits)
        %%TODO: UNUSED, REMOVE    
            trial = obj.trial_by_id(trial_id);
            
            n_samples = range(limits) * trial.fs;
            
            mm_onset = trial.mismatch_onset_t();
            start_idx = find(trial.probe_t > mm_onset + limits(1), 1, 'first');
            full_idx = start_idx + (0:n_samples-1);
            running = trial.velocity(full_idx);
            
            t = limits(1) + (0:n_samples-1)*(1/trial.fs);
        end
        
        
        
        function changed_down = trial_changed_velocity(obj, trial_id, window_1, window_2, n_sds)
        %%TODO: UNUSED, REMOVE    
            limits = [min(window_1), max(window_2)];
            
            [running, t] = running_around_mismatch(obj, trial_id, limits);
            
            % determine whether the trial has changed velocity after
            % mismatch
            baseline_idx = t > window_1(1) & t < window_1(2);
            response_idx = t >= window_2(1) & t < window_2(2);
            
            m_before = mean(running(baseline_idx));
            sd = std(running(baseline_idx));
            
            m_after = mean(running(response_idx));
            
            changed_down = m_after < m_before - n_sds*sd;
        end
        
        
        
        function val = get_trial_protocol_id(obj, trial)
        %%TODO: UNUSED, REMOVE    
            idx_prot = strcmp(obj.protocol_type, trial.protocol);
            idx_gain = strcmp(obj.protocol_gain, trial.config.gain_direction);
            
            val = obj.protocol_ids(idx_prot & idx_gain);
        end
    end
end
