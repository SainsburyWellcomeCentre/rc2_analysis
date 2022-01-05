classdef RVTSession < RC2Session
% RVTSession Class for handling details of a session with R, VF or T
% trials
%
%   RVTSession Properties:
%       trials              - cell array with entries of class Trial
%       use_aligned_data    - default: true
%
%   RVTSession Methods:
%       get_group_label_from_id            - return a trial group label from its ID
%       get_group_id_from_label            - return a trial group ID from its label
%       get_trials_with_trial_group_label  - return a cell array of trials with the specified trial label
%       get_trial_with_id                  - return trials with specified ID (not trial group ID, but trial ID)
%       get_original_trial_for_replay      - for a replay trial, return the original trial
%       match_replay_to_original           - match the replay trial to the original trial
%       get_motion_bouts_for_trial_group   - for all trials in the trial group, get all of the motion bouts
%       apply_offsets                      - applies offsets to trials from the cached table of offsets
%       to_aligned                         - for a trial, return the corresponding AlignedTrial
%
%   All subclasses deal with a particular protocol and should define:
%       trial_group_ids
%       trial_group_labels
%
%   Each R/VF/T trial is assigned to a 'group', which is given a label and
%   an ID. The ID is an integer starting at 1, and the label is a short
%   string.  This class does not have default for those properties, so they
%   must be defined by the subclass.
%   
%   See also: FourWayProtocolSession, HeadTiltMiceSession,
%   LocoVestLocoVestSession, LocoVestLocoVestDarknessSession,
%   MismatchNov2020Session, MismatchJul2021Session,
%   MismatchDarknessOct2021Session, PassiveExperiment,
%   PassiveProtocolSession, SparseNoiseSession
    
    properties (SetAccess = protected)
        
        trials = {}
    end
    
    properties
        
        use_aligned_data = true
    end
    
    
    
    methods
        
        function obj = RVTSession(session)
        % RVTSession
        %
        %   RVTSession(SESSION) takes a session structure saved in the
        %   formatted data and creates an object to handle the data.
        
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
        %%get_group_label_from_id Return a trial group label from its ID
        %
        %   LABEL = get_group_label_from_id(GROUP_ID) 
        
            idx = obj.trial_group_ids == trial_group_id;
            label = obj.trial_group_labels{idx};
        end
        
        
        
        function label = get_group_id_from_label(obj, trial_group_label)
        %%get_group_id_from_label Return a trial group ID from its label
        %
        %   GROUP_ID = get_group_id_from_label(LABEL) 
        
            idx = strcmp(obj.trial_group_labels, trial_group_label);
            label = obj.trial_group_ids(idx);
        end
        
        
        
        function trials = get_trials_with_trial_group_label(obj, trial_group_labels)
        %%get_trials_with_trial_group_label Return a cell array of trials with the specified trial label
        %
        %   TRIALS = get_trials_with_trial_group_label(GROUP_LABELS) 
        %   gets the set of trials with the trial groups labels in GROUP_LABELS.
        %   LABELS can be a single string with the trial group label or a
        %   cell array of group labels, in which case all trials satisfying
        %   any of the LABELS is returned.
            
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
        %%get_trial_with_id Return trials with specified ID (not trial group ID, but trial ID)
        %
        %   TRIAL = get_trial_with_id(TRIAL_ID) 
        %   returns the trial with trial ID, TRIAL_ID. TRIAL is the Trial object. The
        %   supplied ID is *not* the trial group ID, but the trial ID in
        %   trial_id property of the Trial object.
        
            idx = cellfun(@(x)(x.trial_id), obj.trials) == trial_id;
            trial = obj.trials{idx};
        end
        
        
        
        function original_trial = get_original_trial_for_replay(obj, replay_trial_id)
        %%get_original_trial_for_replay For a replay trial, return the original trial
        %
        %   ORIGINAL_TRIAL = get_original_trial_for_replay(REPLAY_TRIAL_ID)
        %   for a trial which is a replay, return the original trial.
        %   REPLAY_TRIAL_ID is an integer refering to the ID of the replay
        %   trial and ORIGINAL_TRIAL is a Trial object referring to the
        %   original trial.
        
            replay_trial = obj.get_trial_with_id(replay_trial_id);
            
            if ~replay_trial.is_replay
                original_trial = [];
                return
            end
            
            original_trial = obj.get_trial_with_id(replay_trial.original_trial_id);
        end
        
        
        
        function offset = match_replay_to_original(obj, trial_id)
        %%match_replay_to_original Match the replay trial to the original trial
        %
        %   OFFSET = match_replay_to_original(TRIAL_ID)
        %   matches a replay trial with trial ID TRIAL_ID to its original,
        %   and returns the number of sample points the replay is offset
        %   from the original. # samples is returned as an integer OFFSET.
        
            
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
        
        
        
        function all_bouts = get_motion_bouts_for_trial_group(obj, trial_group_label, options)
        %%get_motion_bouts_for_trial_group For all trials in a trial
        %%group, get all of the motion bouts 
        %
        %   BOUTS_ARRAY = get_motion_bouts_for_trial_group(GROUP_LABELS, OPTIONS)
        %
        %   Returns all bouts in the trials in GROUP_LABELS. GROUP_LABELS 
        %   can be a single string with a group label, or a cell array of
        %   group labels. See `get_trials_with_trial_group_label`.
        %
        %   OPTIONS can be supplied to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for 
        %                                   a bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low
        %
        %   BOUTS_ARRAY is a cell array of MotionBout objects.
        %
        %   See also: MotionBout, Trial.motion_bouts,
        %   get_trials_with_trial_group_label 
        
            trials = obj.get_trials_with_trial_group_label(trial_group_label); %#ok<*PROPLC>
            
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
        %%apply_offsets Applies the cached table of offsets to the
        %%individual trials
        %
        %   apply_offsets(OFFSETS_TABLE) takes the MATLAB table of offsets
        %   and applies them to all the trials in the session.
        %
        %   See also: FormattedData.apply_offsets
        
        
        %%applies offsets to replay trials from the table of offsets
            for ii = 1 : length(obj.trials)
                idx = obj.trials{ii}.trial_id == offsets_tbl.replay_trial_id;
                if any(idx)
                    obj.trials{ii}.replay_offset = offsets_tbl.sample_offset(idx);
                end
            end
        end
        
        
        
        function out_trial = to_aligned(obj, in_trial)
        %%to_aligned For a trial, return the corresponding AlignedTrial
        %%object
        %
        %   ALIGNED_TRIAL = to_aligned(TRIAL)
        %   for a trial in the session, TRIAL (of class Trial), find the
        %   original trial and create an AlignedTrial object and return as
        %   ALIGNED_TRIAL.
        
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
