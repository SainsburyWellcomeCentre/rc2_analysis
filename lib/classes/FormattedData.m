classdef FormattedData < handle
    
    properties (SetAccess = private)
        
        probe_id
        
        clusters = Cluster.empty()
        anatomy = Anatomy.empty()
        sessions = {}
        
        svm_table
        replay_offsets
    end
    
    properties (Hidden = true)
        
        ctl
        data
    end
    
    
    
    methods
        
        function obj = FormattedData(ctl, data)
        %%methods
        %   selected_clusters
        %   VISp_clusters
        %   VISp_cluster_ids
        %   spiking_class_ids   
        
            obj.ctl = ctl;
            obj.data = data;
            
            obj.probe_id = obj.data.probe_id;
            
            % create cluster object for each cluster
            for ii = 1 : length(obj.data.clusters)
                obj.clusters(ii) = Cluster(obj.data.clusters(ii));
            end
            
            obj.anatomy = Anatomy(obj.data.anatomy);
            
            for ii = 1 : length(obj.data.sessions)
                
                protocol = obj.ctl.get_protocol_from_session_id(obj.data.probe_id, obj.data.sessions(ii).session_id);
                
                switch protocol
                    case 'locovest_loco_vest'
                        obj.sessions{ii} = LocoVestLocoVestSession(obj.data.sessions(ii));
                    case 'locovest_loco_vest_darkness'
                        obj.sessions{ii} = LocoVestLocoVestDarknessSession(obj.data.sessions(ii));
                    case 'four_way_protocol'
                        obj.sessions{ii} = FourWayProtocolSession(obj.data.sessions(ii));
                    case 'head_tilt_mice'
                        obj.sessions{ii} = HeadTiltMiceSession(obj.data.sessions(ii));
                    case 'mismatch_nov2020'
                        obj.sessions{ii} = MismatchNov2020Session(obj.data.sessions(ii));
                    case 'mismatch_jul2021'
                        obj.sessions{ii} = MismatchJul2021Session(obj.data.sessions(ii));
                    case 'mismatch_darkness_oct2021'
                        obj.sessions{ii} = MismatchDarknessOct2021Session(obj.data.sessions(ii));
                    case 'sparse_noise'
                        obj.sessions{ii} = SparseNoiseSession(obj.data.sessions(ii));
                    case 'sf_tf'
                        obj.sessions{ii} = SparseNoiseSession(obj.data.sessions(ii));
                end
            end
            
            obj.svm_table = obj.ctl.load_svm_table(obj.data.probe_id);
            obj.replay_offsets = obj.ctl.load_offsets_table(obj.data.probe_id);
            
            if isempty(obj.replay_offsets)
                return
            end
            
            obj.apply_offsets();
        end
        
        
        
        function tuning_curve = load_tuning_curves(obj, cluster_id, trial_group)
            
            tbl = obj.ctl.load_tuning_curves(obj.probe_id);
            group_idx = cellfun(@(x)(isequal(x, trial_group)), tbl.trial_groups);
            tuning_curves = tbl.tuning_curves{group_idx};
            cluster_idx = [tuning_curves(:).cluster_id] == cluster_id;
            tuning_curve = tuning_curves(cluster_idx);
        end
        
        
        
        function apply_offsets(obj)
        %%applies offset sample for a replay trial
            sessions = obj.motion_sessions();
            for ii = 1 : length(sessions)
                sessions{ii}.apply_offsets(obj.replay_offsets);
            end
        end
        
        
        
        function session = get_session_with_id(obj, session_id)
        %%returns session object with id
            session_ids = cellfun(@(x)(x.session_id), obj.sessions, 'UniformOutput', false);
            idx = strcmp(session_ids, session_id);
            session = obj.sessions{idx};
        end
        
        
        
        function sessions = motion_sessions(obj)
        %%gets the sessions which have linear motion of some kind
        %   TODO: this is suboptimal... would have to change every time
        %   experiment is done and different session order occurs
        
            if strcmp(obj.probe_id, 'CAA-1112872_rec1_rec1b_rec2_rec3')
                sessions = obj.sessions(1:2);
            else
                sessions = obj.sessions(1);
            end
        end
        
        
        
        function trials = motion_trials(obj)
        %%gets the trials which are have linear motion of some kind
            
            sessions = obj.motion_sessions();
            t = cellfun(@(x)(x.trials), sessions, 'UniformOutput', false);
            trials = [t{:}]; %#ok<*PROP>
        end
        
        
        
        function trials = get_trials_with_trial_group_label(obj, trial_group_label)
        %%gets the trials which have trial group label
            
            sessions = obj.motion_sessions(); %#ok<*PROPLC>
            
            trials = {};
            for ii = 1 : length(sessions)
                these_trials = sessions{ii}.get_trials_with_trial_group_label(trial_group_label);
                trials = [trials, these_trials]; %#ok<*AGROW>
            end
        end
        
        
        
        function trials = get_trials_with_trial_ids(obj, trial_ids)
        %%return cell array of trials with trial IDs specified by
        % 'trial_ids'
        % If a single number is given returns the trial itself rather than
        % a cell array of length 1.
        
            trials = obj.motion_trials();
            idx = cellfun(@(x)(ismember(x.trial_id, trial_ids)), trials);
            trials = trials(idx);
            
            if length(trials) == 1
                trials = trials{1};
            end
        end
        
        
        
        function cluster = get_cluster_with_id(obj, cluster_id)
        %%returns cluster object for cluster with id = 'cluster_id'
            idx = [obj.clusters(:).id] == cluster_id;
            cluster = obj.clusters(idx);
        end
        
        
        
        function selected_clusters = selected_clusters(obj, spiking_class)
        %%retun cluster objects for manually selected clusters
        %   Inputs:     spiking_class - 'any' (default), 'narrow', 'wide'
        %   Outputs:    array of cluster objects
        
            VariableDefault('spiking_class', 'any');
            
            search_ids = intersect(obj.data.selected_clusters, obj.spiking_class_ids(spiking_class));
            idx = ismember([obj.clusters(:).id], search_ids);
            selected_clusters = obj.clusters(idx);
        end
        
        
        
        function ids = selected_cluster_ids(obj, spiking_class)
        %%eturn the ids of selected clusters
        %   Inputs:     spiking_class - 'any' (default), 'narrow', 'wide'
        %   Outputs:    list of cluster ids
        
            VariableDefault('spiking_class', 'any');
            
            selected_clusters = obj.selected_clusters(spiking_class);
            ids = [selected_clusters(:).id];
        end
        
        
        
        function visp_clusters = VISp_clusters(obj, is_selected, spiking_class)
        %%return cluster objects for VISp clusters
        %   Inputs:     is_selected - boolean, indicates whether to return
        %                             only "selected" VISp clusters (default) or all
        %               spiking_class - 'any' (default), 'narrow' or 'wide'
        %   Outputs:    array of cluster objects
        
            VariableDefault('is_selected', true);
            VariableDefault('spiking_class', 'any');
            
            if is_selected
                selected_clusters = obj.selected_clusters(spiking_class);
            else
                selected_clusters = obj.clusters;
            end
            
            % look for clusters in VISp
            idx = regexp({selected_clusters(:).region_str}, 'VISp[\dX]');
            idx = ~cellfun(@isempty, idx);
            
            % restrict clusters
            visp_clusters = selected_clusters(idx);
            
            % further restrict based on spiking class
            if ismember(spiking_class, {'wide', 'narrow'})
                idx = strcmp({visp_clusters(:).spiking_class}, spiking_class);
                visp_clusters = visp_clusters(idx);
            end
        end
        
        
        
        function ids = VISp_cluster_ids(obj, is_selected, spiking_class)
        %%return the ids of VISp clusters
        %   Inputs:     is_selected - boolean, indicates whether to return
        %                             only "selected" VISp cluster ids or all
        %               spiking_class - 'any', 'narrow' or 'wide'
        %   Outputs:    uint32 list of cluster ids in VISp
        
            VariableDefault('is_selected', true);
            VariableDefault('spiking_class', 'any');
            
            visp_clusters = obj.VISp_clusters(is_selected, spiking_class);
            ids = [visp_clusters(:).id];
        end
        
        
        
        function ids = spiking_class_ids(obj, spiking_class)
        %%return cluster ids for a particular spiking class
        %   Input:  spiking_class - 'any' (default), 'narrow', or 'wide'
        %   Output: uint32 list of cluster ids for spiking class
            
            VariableDefault('spiking_class', 'any');
        
            ids = [obj.clusters(:).id];
            
            if ~strcmp(spiking_class, 'any')
                all_spiking_classes = {obj.clusters(:).spiking_class};
                idx = strcmp(all_spiking_classes, spiking_class);
                ids = ids(idx);
            end
        end
        
        
        
        function [relative_depth, layer] = get_relative_layer_depth_of_cluster(obj, cluster_id)
        %%returns the relative depth of the cluster in a VISp layer
            
            cluster = obj.get_cluster_with_id(cluster_id);
            if ~cluster.is_VISp; relative_depth = -1; layer = ''; return; end
            [relative_depth, layer] = obj.anatomy.VISp_layer_relative_depth(cluster.distance_from_probe_tip);
        end
        
        
        
        function bouts = get_motion_bouts_for_trial_group(obj, trial_group_label, options)
        %%returns an array of 'motion bouts' for a particular
        %%trial_group_label
        %
        %   'options' can be supplied to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low
        
            sessions = obj.motion_sessions(); %#ok<*PROPLC>
            
            bouts = {};
            for ii = 1 : length(sessions)
                these_bouts = sessions{ii}.get_motion_bouts_for_trial_group(trial_group_label, options);
                bouts = [bouts, these_bouts]; %#ok<*AGROW>
            end
        end
        
        
        
        function times = get_motion_bout_start_times(obj, trial_group_label)
        %%returns the start times of all motion bouts for a trial group
        
            bouts = obj.get_motion_bouts_for_trial_group(trial_group_label);
            times = cellfun(@(x)(x.start_time), bouts);
        end
        
        
        
        function [sig, p, direction, stat_med, mot_med, n] = is_stationary_vs_motion_significant(obj, cluster_id, trial_group_label)
        %%calculates whether cluster_id is signficantly modulation by
        %%motion for the trial group type specified by trial_group_label
        %   returns sig - is there a significant modulation for motion
        %           p - the p value of that result
        %           direction - whether it increases (1), decreases (-1) or
        %           no change (0)
        
            valid = obj.check_trial_group(trial_group_label);
            if ~valid
                warning('no such trial group label'); 
                sig = []; p = []; direction = []; stat_med = []; mot_med = []; n = [];
                return
            end
            
            stationary_fr   = obj.stationary_fr_for_trial_group(cluster_id, trial_group_label);
            motion_fr       = obj.motion_fr_for_trial_group(cluster_id, trial_group_label);
            
            idx = isnan(stationary_fr) | isnan(motion_fr);
            
            stationary_fr(idx) = [];
            motion_fr(idx) = [];
            
            [~, ~, p, direction] = compare_groups_with_signrank(stationary_fr, motion_fr);

            sig = p < 0.05;
            
            stat_med = median(stationary_fr);
            mot_med = median(motion_fr);
            n = length(stationary_fr);
        end
        
        
        
        function [sig, p, direction, x_med, y_med, n] = is_motion_vs_motion_significant(obj, cluster_id, trial_group_label_x, trial_group_label_y, max_n_trials)
        %%calculates whether cluster_id has significantly larger motion
        %%modulation for one trial group than another
        %   returns sig - is there a significant difference between motion
        %   in the two groups
        %           p - the p value of that result
        %           direction - whether it increases (1), decreases (-1) or
        %           no change (0)  (note that increase = larger FR value in
        %           group y than group x)
        
            VariableDefault('max_n_trials', []);
        
            valid = obj.check_trial_group(trial_group_label_x) & obj.check_trial_group(trial_group_label_y);
            if ~valid
                warning('no such trial group label'); 
                sig = []; p = []; direction = []; x_med = []; y_med = []; n = [];
                return
            end
            
            [x_fr, x_trial_id] = obj.motion_fr_for_trial_group(cluster_id, trial_group_label_x);
            [y_fr, y_trial_id] = obj.motion_fr_for_trial_group(cluster_id, trial_group_label_y);
            
            n_trials = min([length(x_fr), length(y_fr)]);
            
            x_fr = x_fr(1:n_trials);
            y_fr = y_fr(1:n_trials);
            
            assert(all(abs(y_trial_id(1:n_trials) - x_trial_id(1:n_trials)) < length(obj.sessions{1}.trial_group_ids)));
            
            idx = isnan(x_fr) | isnan(y_fr);
            
            x_fr(idx) = [];
            y_fr(idx) = [];
            
            assert(length(x_fr) == length(y_fr));
            n_trials = length(x_fr);
            
            if ~isempty(max_n_trials)
                n_trials = min(n_trials, max_n_trials);
            end
            
            x_fr = x_fr(1:n_trials);
            y_fr = y_fr(1:n_trials);
            
            [~, ~, p, direction] = compare_groups_with_signrank(x_fr, y_fr);

            sig = p < 0.05;
            
            x_med = median(x_fr);
            y_med = median(y_fr);
            n = length(x_med);
        end
        
        
        
        function [fr, trial_id] = stationary_fr_for_trial_group(obj, cluster_id, trial_group_label)
        %%returns stationary firing rates for cluster_id for trials part of
        %%trial group label
        %   returns vector of firing rates and trial id from where it came
            
            valid = obj.check_trial_group(trial_group_label);
            if ~valid; fr = []; return; end
        
            idx = ismember(obj.svm_table.trial_group_label, trial_group_label) & ...
                  obj.svm_table.cluster_id == cluster_id;
            
            fr = obj.svm_table.stationary_fr(idx);
            trial_id = obj.svm_table.trial_id(idx);
        end
        
        
        
        function [fr, trial_id] = motion_fr_for_trial_group(obj, cluster_id, trial_group_label)
        %%returns motion firing rates for cluster_id for trials part of
        %%trial group label
        %   returns vector of firing rates and trial id from where it came
        
            valid = obj.check_trial_group(trial_group_label);
            if ~valid; fr = []; return; end
            
            idx = ismember(obj.svm_table.trial_group_label, trial_group_label) & ...
                  obj.svm_table.cluster_id == cluster_id;
            
            fr = obj.svm_table.motion_fr(idx);
            trial_id = obj.svm_table.trial_id(idx);
        end
        
        
        
%         function [fr, sd, n, x, shuff, stat_fr, stat_sd, stat_n] = tuning_curve_info(obj, cluster_id, trial_group_label)
%             
%             tbl = obj.load_tuning_curves();
%             
%             idx = tbl.cluster_id == cluster_id;
%             
%             tuning = tbl.tuning{idx};
%             timing = tbl.timing{idx};
%             
%             bin_edges = tbl.bin_edges{idx};
%             stat_fr = tbl.stationary_fr{idx};
%             
%             fr = nanmean(tuning, 2);
%             sd = nanstd(tuning, [], 2);
%             n = sum(~isnan(tuning), 2);
%             x = (bin_edges(1:end-1) + bin_edges(2:end))/2;
%             
%             shuff = ShuffleTuning(tuning, x);
%             
% %             p_anova = s.p;%anova1(tuning', [], 'off');
% %             beta = s.beta(1);
%             
%             stat_fr = nanmean(stat_fr);
%             stat_sd = nanstd(stat_fr);
%             stat_n = sum(~isnan(stat_fr));
%         end
        
        
        
        function valid = check_trial_group(obj, trial_group_label)
        %%is the trial group in the data
            sessions = obj.motion_sessions();
            t = cellfun(@(x)(x.trial_group_labels), sessions, 'uniformoutput', false);
            valid = ismember(trial_group_label, [t{:}]);
        end
        
        
        
        function [sig, p, direction] = is_mismatch_significant(obj, cluster_id, trial_group_label)
        %%returns whether the mismatch response for cluster_id is
        %%significant in trial group 'trial_group_label'
        
            valid = obj.check_trial_group(trial_group_label);
            if ~valid; sig = []; p = []; direction = []; return; end
            
            mm = MismatchAnalysis();
            
            cluster = obj.get_cluster_with_id(cluster_id);
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            [sig, p, direction] = mm.is_response_significant(cluster, trials);
        end
        
        
        
        function delta_fr = get_mismatch_response(obj, cluster_id, trial_group_label)
        %%returns magnitude of response to mismatch
        
            valid = obj.check_trial_group(trial_group_label);
            if ~valid; delta_fr = []; return; end
            
            mm = MismatchAnalysis();
            
            cluster = obj.get_cluster_with_id(cluster_id);
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            
            baseline = mm.get_avg_baseline_fr(cluster, trials);
            response = mm.get_avg_response_fr(cluster, trials);
            delta_fr = response - baseline;
        end
        
        
        
        function times = get_mismatch_onset_times(obj, trial_group_label)
        %%returns the onset times of the mismatch for trials in the trial
        %%group
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            times = cellfun(@(x)(x.mismatch_onset_t), trials);
        end
        
        
        
        function [times, bouts] = get_motion_onset_times(obj, trial_group_label, options)
        %%returns the onset times of motion bouts for a trial group
        %
        %   'options' can be supplied to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low
        
            bouts = obj.get_motion_bouts_for_trial_group(trial_group_label, options);
            times = cellfun(@(x)(x.start_time), bouts);
        end
        
        
        
        function [traces, t, times] = get_traces_around_mismatch_onset(obj, trial_group_label, limits, common_fs)
        %%for all mismatch trials in a trial group, get traces around the
        %%mismatch, padding specified by 'limits' and sampled at 'common_fs' frequency
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            times = obj.get_mismatch_onset_times(trial_group_label);
            [traces, t] = obj.get_traces_around_times(trials, times, limits, common_fs);
        end
        
        
        
        function [traces, t, times] = get_traces_around_solenoid_up(obj, trial_group_label, limits, common_fs)
        %%for all mismatch trials in a trial group, get traces around the
        %%mismatch, padding specified by 'limits' and sampled at 'common_fs' frequency
        
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            idx = cellfun(@(x)(find(diff(x.solenoid > 2.5) == 1) + 1), trials, 'uniformoutput', false);
            times = cellfun(@(x, y)(x.probe_t(y)), trials, idx);    
            [traces, t] = obj.get_traces_around_times(trials, times, limits, common_fs);
        end
        
        
        
        function [traces, t, times] = get_traces_around_motion_onset(obj, trial_group_label, limits, common_fs)
        %%for all motion trials in a trial group, get traces around the
        %%motion onsets, padding specified by 'limits' and sampled at 'common_fs' frequency
            [times, bouts] = obj.get_motion_onset_times(trial_group_label);
            trials = cellfun(@(x)(x.trial), bouts, 'uniformoutput', false);
            [traces, t] = obj.get_traces_around_times(trials, times, limits, common_fs);
        end
        
        
        
        function [fr, t] = get_fr_responses(obj, cluster_id, trial_group_label, type, limits, fs, options)
        %%returns firing rate convolutions around events indicated by
        %%'type' (motion or mismatch), for trials in trial_group_label,
        %%between limits (1x2 vector) sampled at fs.
        %
        %   'options' can be supplied to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low
        
            if strcmp(type, 'motion')
                times = obj.get_motion_onset_times(trial_group_label, options);
            elseif strcmp(type, 'mismatch')
                times = obj.get_mismatch_onset_times(trial_group_label);
            end
            
            [fr, t] = obj.get_fr_around_times(cluster_id, times, limits, fs);
        end
        
        
        
        function tbl = create_svm_table(obj)
        %%creates table with firing rate for each cluster during the
        %%stationary and motion periods
        
            headers = {'probe_id', ...
                       'session_id', ...
                       'trial_id', ...
                       'trial_group_label', ...
                       'cluster_id', ...
                       'stationary_fr', ...
                       'stationary_time', ...
                       'motion_fr', ...
                       'motion_time'};
                   
            tbl = cell2table(cell(0, length(headers)), 'VariableNames', headers);
            
            all_trials = obj.motion_trials();
            
            for ii = 1 : length(all_trials)
                
                trial               = all_trials{ii}.to_aligned();
                
                stationary_mask     = trial.stationary_mask;
                stationary_time     = trial.stationary_time;
                
                motion_mask         = trial.motion_mask;
                motion_time         = trial.motion_time;
            
                for cluster = obj.selected_clusters
                    
                    fr_conv = cluster.fr.get_convolution(trial.probe_t);
                    
                    new_row = {trial.probe_id, ...
                               trial.session_id, ...
                               trial.trial_id, ...
                               trial.trial_group_label, ...
                               cluster.id, ...
                               mean(fr_conv(stationary_mask)), ...
                               stationary_time, ...
                               mean(fr_conv(motion_mask)), ...
                               motion_time};
                    
                    tbl = [tbl; new_row];
                end
            end
        end
        
        
        
        function tbl = create_replay_offsets_table(obj)
        %%creates table with sample offset between replay trials and the original
        %%trials
        
            headers = {'original_trial_id', ...
                       'replay_trial_id', ...
                       'sample_offset'};
            tbl = cell2table(cell(0, length(headers)), 'VariableNames', headers);
            
            % get trials which are replays
            all_trials = obj.motion_trials();
            replay_trials = all_trials(cellfun(@(x)(x.is_replay), all_trials));
            
            for ii = 1 : length(replay_trials)
                
                trial = replay_trials{ii};
                if isempty(trial.original_trial_id); continue; end
                session = obj.get_session_with_id(trial.session_id);
                offset = session.match_replay_to_original(trial.trial_id);
                
                new_row = {trial.original_trial_id, ...
                           trial.trial_id, ...
                           offset};
                
                tbl = [tbl; new_row];
            end
        end
        
        
        
        function tbl = create_tuning_curves(obj, trial_types)
            
            tbl = cell(1, length(trial_types));
            
            for ii = 1 : length(trial_types)
                
                tt = TuningTable(obj.probe_id);
                
                % get all trials of type specified in the cell array
                % 'trial_types'
                trials = obj.get_trials_with_trial_group_label(trial_types{ii});
                
                % if the trials are replays, align them to the original
                for jj = 1 : length(trials)
                    aligned_trials{jj} = trials{jj}.to_aligned;
                end
                
                % add these trials to the TuningTable object
                tt.add_trials(aligned_trials, trial_types{ii});
                
                % all selected clusters add them to the tuning table
                clusters = obj.selected_clusters();
                for jj = 1 : length(clusters)
                    tbl{ii}(jj) = tt.tuning_curve(clusters(jj));
                end
            end
        end
        
        
        
        function [traces, t] = get_traces_around_times(obj, trials, times, limits, fs)
        %%gets traces centered around times in 'times' using trials in
        %%'trials' (must be same length as times and matched to times).
        %%Returns traces limits(1) to limits(2) around those times sampled
        %%at frequency common_fs
        
            t = obj.timebase(limits, fs);
            traces = nan(length(t), length(trials));
            
            for ii = 1 : length(trials)
                these_limits = limits + times(ii);
                this_t = t + times(ii);
                start_idx = find(trials{ii}.probe_t < these_limits(1), 1, 'last');
                end_idx = find(trials{ii}.probe_t > these_limits(2), 1, 'first');
                idx = trials{ii}.probe_t >= trials{ii}.probe_t(start_idx) & ...
                      trials{ii}.probe_t <= trials{ii}.probe_t(end_idx);
%                 idx = trials{ii}.probe_t >= these_limits(1) & ...
%                       trials{ii}.probe_t <= these_limits(2);
                traces(:, ii) = interp1(trials{ii}.probe_t(idx), ...
                                        trials{ii}.velocity(idx), ...
                                        this_t);
            end
        end
        
        
        
        function [fr, t] = get_fr_around_times(obj, cluster_id, times, limits, fs)
        %%returns firing rate convolutions around the times 'times' (Nx1
        %%vector) for the cluster 'cluster_id', spanning limits at a sample
        %%frequency of common_fs
            
            cluster = obj.get_cluster_with_id(cluster_id);
            
            t = obj.timebase(limits, fs);
            fr = nan(length(times), length(t));
            
            for ii = 1 : length(times)
                fr(ii, :) = cluster.fr.get_convolution(times(ii) + t);
            end
        end
    end
    
    
    
    methods (Static = true)
        
        function val = timebase(limits, fs)
        %%returns a timebase running between limits (1x2 vector) sampled at fs
            val = linspace(limits(1), limits(2), round((limits(2) - limits(1))*fs)+1);
        end
        
        
        
        function idx = time_idx_around_trigger_time(trial, trigger_t, limits)
            
            idx = trial.probe_t >= trigger_t + limits(1) & ...
                    trial.probe_t <= trigger_t + limits(2);
        end
    end
end
