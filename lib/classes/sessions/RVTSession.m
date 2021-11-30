classdef RVTSession < RC2Session
    
    
    properties (SetAccess = protected)
        
        trials = {}
    end
    
    properties
        
        use_aligned_data = true
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
            
            % offset, of large gaps, zero indexed
            offset_correction = max(start_idx(max_corr_idx) - 500, 1) - 1;
            
            % next loop over more carefully
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
        
        
        
        function all_bouts = get_motion_bouts_for_trial_group(obj, trial_group_id, options)
        %   'options' can be supplied to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low
        
            trials = obj.get_trials_with_trial_group_label(trial_group_id); %#ok<*PROPLC>
            
            all_bouts = {};
            for ii = 1 : length(trials)
                
                if obj.use_aligned_data
                    trial = obj.to_aligned(trials{ii});
                else
                    trial = trials{ii};
                end
                
                motion_bouts = trial.motion_bouts(options.include_200ms);
                all_bouts = [all_bouts, motion_bouts];
            end
            
            valid = cellfun(@(x)(x.duration > options.min_bout_duration), all_bouts);
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
    end
end
