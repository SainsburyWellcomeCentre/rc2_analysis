classdef RVTSession < RC2Session
    
    
    properties (SetAccess = protected)
        
        trials = {}
    end
    
    properties
        
        include_200ms = true
        use_aligned_data = true
        min_bout_duration = 2  % seconds
    end
    
    
    
    methods
        
        function obj = RVTSession(session)
            
            obj = obj@RC2Session(session);
            
            for ii = 1 : session.n_trials
                
                % create trial object (may have to subclass this at one
                % point)
                obj.trials{ii} = Trial(session.trials(ii), obj);
                
                % get which group the trial belongs to
                trial_group_id = obj.get_trial_group_id(obj.trials{ii});
                
                % set the group id and label
                obj.trials{ii}.trial_group_id = trial_group_id;
                obj.trials{ii}.trial_group_label = obj.get_group_label_from_id(trial_group_id);
            end
        end
        
        
        
        function label = get_group_label_from_id(obj, trial_group_id)
        %%returns trial group label from the trial group id
            idx = obj.trial_group_ids == trial_group_id;
            label = obj.trial_group_labels{idx};
        end
        
        
        
        function label = get_group_id_from_label(obj, trial_group_label)
        %%returns trial group id from the trial group label
            idx = strcmp(obj.trial_group_labels, trial_group_label);
            label = obj.trial_group_ids(idx);
        end
        
        
        
        function trials = get_trials_with_trial_group_label(obj, trial_group_labels)
        %%returns trials with trial group label (can be cell array of group
        %%labels)
            
            if ~iscell(trial_group_labels)
                trial_group_labels = {trial_group_labels};
            end
            
            trial_array = [obj.trials{:}];
            
            idx = false(1, length(trial_array));
            for ii = 1 : length(trial_group_labels)
                idx = idx | strcmp({trial_array(:).trial_group_label}, trial_group_labels{ii});
            end
            
            trials = obj.trials(idx);
        end
        
        
        
        function trial = get_trial_with_id(obj, trial_id)
        %%returns the trial with id (not to be confused with 'group_id'
            idx = cellfun(@(x)(x.trial_id), obj.trials) == trial_id;
            trial = obj.trials{idx};
        end
        
        
        
        function original_trial = get_original_trial_for_replay(obj, replay_trial_id)
        %%for a trial which is a replay, return the original trial
        
            replay_trial = obj.get_trial_with_id(replay_trial_id);
            
            if ~replay_trial.is_replay
                original_trial = [];
                return
            end
            
            original_trial = obj.get_trial_with_id(replay_trial.original_trial_id);
        end
        
        
        
        function offset = match_replay_to_original(obj, trial_id)
        %%matches a replay trial to its original, and returns the amount
        %%the replay trial is offset
            
            original_trial = obj.get_original_trial_for_replay(trial_id);
            replay_trial = obj.get_trial_with_id(trial_id);
            
            %use the multiplexer output if its available, otherwise use the
            %filtered teensy 
            if ~isempty(original_trial.multiplexer_output)
                ori_trace = original_trial.multiplexer_output;
                rep_trace = replay_trial.multiplexer_output;
            else
                ori_trace = original_trial.filtered_teensy;
                rep_trace = max(replay_trial.stage, 0);
            end
            
            n_ori_samples = length(ori_trace);
            n_rep_samples = length(rep_trace);
            
            % roughly loop over the traces first
            start_idx = 1 : 1000 : max(n_ori_samples, n_rep_samples);
            
            r = nan(1, length(start_idx));
            
            % for each shift
            for i = 1 : length(start_idx)
                
                if n_ori_samples < n_rep_samples
                    m = min(n_ori_samples, n_rep_samples - start_idx(i) + 1);
                    o = ori_trace(1:m);
                    O = rep_trace(start_idx(i) + (0:m-1));
                else
                    m = max(1, n_rep_samples - start_idx(i) + 1);
                    M = min(n_rep_samples, start_idx(i));
                    o = ori_trace(1:m);
                    O = rep_trace(M:end);
                end
                
                r(i) = corr(o, O);
            end
            
            [~, max_corr_idx] = max(r);
            
            % next loop over more carefully
            offset_correction = start_idx(max_corr_idx) - 500;
            start_idx = start_idx(max_corr_idx) + (-500:499);
            start_idx(start_idx<1) = [];
            start_idx(start_idx>max(n_ori_samples, n_rep_samples)) = [];
            
            r = nan(1, length(start_idx));
            for i = 1 : length(start_idx)
                
                if n_ori_samples < n_rep_samples
                    m = min(n_ori_samples, n_rep_samples - start_idx(i) + 1);
                    o = ori_trace(1:m);
                    O = rep_trace(start_idx(i) + (0:m-1));
                else
                    m = n_rep_samples - start_idx(i) + 1;
                    o = ori_trace(1:m);
                    O = rep_trace(start_idx(i):end);
                end
                    
                r(i) = corr(o, O);
            end
            
            [~, offset] = max(r);
            offset = offset_correction + offset - 1;
        end
        
        
        
        function all_bouts = get_motion_bouts_for_trial_group(obj, trial_group_id)
            
            trials = obj.get_trials_with_trial_group_label(trial_group_id); %#ok<*PROPLC>
            
            all_bouts = {};
            for ii = 1 : length(trials)
                
                if obj.use_aligned_data
                    trial = obj.to_aligned(trials{ii});
                else
                    trial = trials{ii};
                end
                
                motion_bouts = trial.motion_bouts(obj.include_200ms);
                all_bouts = [all_bouts, motion_bouts];
            end
            
            valid = cellfun(@(x)(x.duration > obj.min_bout_duration), all_bouts);
            all_bouts = all_bouts(valid);
        end
        
        
        
        function apply_offsets(obj, offsets_tbl)
        %%applies offsets to replay trials from the table of offsets
            for ii = 1 : length(obj.trials)
                idx = obj.trials{ii}.trial_id == offsets_tbl.replay_trial_id;
                if any(idx)
                    obj.trials{ii}.replay_offset = offsets_tbl.sample_offset(idx);
                end
            end
        end
        
        
        
        function out_trial = to_aligned(obj, in_trial)
        %%for a trial in the session find the original trial and create an
        %%AlignedTrial object
        
            original_id = in_trial.original_trial_id;
            
            if ~isempty(original_id)
                original_trial = obj.get_trial_with_id(original_id);
                out_trial = AlignedTrial(in_trial, original_trial);
            else
                out_trial = AlignedTrial(in_trial, in_trial);
            end
        end
        
        
        
        
        
        
        
        
        
        
       
        
        
        
        
        
        
        
        function protocol_list = list_protocols(obj)
            
            protocol_list = arrayfun(@(x)(x.protocol), obj.trials, 'uniformoutput', false);
        end
        
        
        
        function protocol_list = list_replay_protocols(obj)
            
            protocol_list = arrayfun(@(x)(x.replay_of), obj.trials, 'uniformoutput', false);
        end
        
        
        
        function log_fnames = list_log_fnames(obj)
            
             log_fnames = arrayfun(@(x)(x.config.log_fname), obj.trials, 'uniformoutput', false);
        end
        
        
        
        function wave_fnames = list_wave_fnames(obj)
            
             wave_fnames = arrayfun(@(x)(x.config.wave_fname), obj.trials, 'uniformoutput', false);
        end
        
        
        
        
        
        
        
        
        
        function load_stationary_vs_motion_table(obj)
            
            csv_fname = fullfile(obj.config.summary_data_dir, 'stationary_vs_motion_fr', sprintf('%s.csv', obj.probe_fname));
            if exist(csv_fname, 'file')
                obj.svm_table = readsvmtable(csv_fname);
            end
        end
        
        
        
        function load_tuning_table(obj)
            
            mat_fname = fullfile(obj.config.summary_data_dir, 'tuning_table', sprintf('%s.mat', obj.probe_fname));
            if exist(mat_fname, 'file')
                t = load(mat_fname);
                obj.tt_table = t.tuning_table;
            end
        end
        
        
        
        function fr = trial_stationary_fr(obj, cluster_id, protocol_id)
            
            idx = obj.get_svm_table_index(cluster_id, protocol_id);
            fr = obj.svm_table.stationary_firing_rate(idx);
            
        end
        
        
        
        function fr = trial_motion_fr(obj, cluster_id, protocol_id)
            
            idx = obj.get_svm_table_index(cluster_id, protocol_id);
            fr = obj.svm_table.motion_firing_rate(idx);
            
        end
        
        
        
        function idx = get_svm_table_index(obj, cluster_id, protocol_id)
            
            trial_type = obj.protocol_type{obj.protocol_ids == protocol_id};
            replayed_type = obj.protocol_replayed_type{obj.protocol_ids == protocol_id};
            
            if strcmp(replayed_type, 'any') || isempty(replayed_type)
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    strcmp(obj.svm_table.protocol, trial_type);
            else
                idx = obj.svm_table.cluster_id == cluster_id & ...
                    strcmp(obj.svm_table.protocol, trial_type) & ...
                    strcmp(obj.svm_table.replay_of, replayed_type);
            end
        end
        
        
        
        function [fr, std, n, x, shuff, stat_fr, stat_sd, stat_n] = ...
                tuning_curve(obj, cluster_id, protocol_id)
            
            idx = obj.get_tuning_table_index(cluster_id, protocol_id);
            
            tuning = obj.tt_table.tuning{idx};
            timing = obj.tt_table.timing{idx};
            
            bin_edges = obj.tt_table.bin_edges{idx};
            stat_fr = obj.tt_table.stationary_fr{idx};
            
            fr = nanmean(tuning, 2);
            std = nanstd(tuning, [], 2);
            n = sum(~isnan(tuning), 2);
            x = (bin_edges(1:end-1) + bin_edges(2:end))/2;
            
            shuff = ShuffleTuning(tuning, x);
            
%             p_anova = s.p;%anova1(tuning', [], 'off');
%             beta = s.beta(1);
            
            stat_fr = nanmean(stat_fr);
            stat_sd = nanstd(stat_fr);
            stat_n = sum(~isnan(stat_fr));
        end
        
        
        
        function idx = get_tuning_table_index(obj, cluster_id, protocol_id)
            
            trial_type = obj.protocol_type{obj.protocol_ids == protocol_id};
            replayed_type = obj.protocol_replayed_type{obj.protocol_ids == protocol_id};
            
            if strcmp(replayed_type, 'any') || isempty(replayed_type)
                idx = obj.tt_table.cluster_id == cluster_id & ...
                    strcmp(obj.tt_table.protocol, trial_type);
            else
                idx = obj.tt_table.cluster_id == cluster_id & ...
                    strcmp(obj.tt_table.protocol, trial_type) & ...
                    strcmp(obj.tt_table.replay_of, replayed_type);
            end
        end
    end
end