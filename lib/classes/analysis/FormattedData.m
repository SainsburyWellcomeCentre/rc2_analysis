classdef FormattedData < handle
% FormattedData Class for handling the formatted data for an experiment and
% extracting important features of the data
%
%   FormattedData Properties:
%       probe_id            - string, the probe recording ID
%       clusters            - #clusters x 1 object array of Cluster objects
%       anatomy             - #probe shanks x 1 object array of Anatomy objects
%       sessions            - #sessions x 1 cell array of Session objects
%       svm_table           - associated MATLAB table of stationary and motion firing rates for 
%                             each cluster during each trial
%       ctl                 - instance of class RC2Analysis
%       data                - contains all the structures in the formatted data file
%
%   FormattedData Methods:
%       list_trial_group_labels             - returns a list of the trial group labels available for this data
%       load_tuning_curves                  - loads previously created .mat files with the tuning curve data for speed
%       load_tuning_curves_acceleration     - loads previously created .mat files with the tuning curve data for acceleration
%       apply_offsets                       - applies the offsets saved in the offsets .csv files to the trials
%       get_session_with_id                 - returns a session object with a session ID
%       motion_sessions                     - returns a cell array of Session objects which are part of an RVTSession
%       motion_trials                       - returns a cell array of Trial objects which are part of an RVTSession
%       get_trials_with_trial_group_label   - returns a cell array of Trials objects with a specified `trial_group_label`
%       get_trials_with_trial_ids           - returns the trial or trials with the specified trial IDs
%       get_cluster_with_id                 - returns the Cluster object for the cluster with integer ID
%       selected_clusters                   - returns an object array of class Clusters with the clusters selected 
%                                             on the manual sorting step (Phy)
%       selected_cluster_ids                - returns a vector of integer IDs corresponding to the clusters selected 
%                                             on the manual sorting ste (Phy)
%       VISp_clusters                       - returns an object array of class Clusters with the clusters which are 
%                                             allocated to a VISp layer
%       VISp_cluster_ids                    - returns a vector of integer IDs corresponding to the clusters selected 
%                                             on the manual sorting step and which are allocated to a VISp layers
%       spiking_class_ids                   - returns a vector of integer IDs corresponding to the clusters with a 
%                                             particular spiking class
%       get_relative_layer_depth_of_cluster - returns the relative depth of a cluster in a VISp layer
%       get_motion_bouts_for_trial_group    - returns a cell array of MotionBout objects extracted from the trials 
%                                             of a particular type
%       is_stationary_vs_motion_significant - for a cluster and a set of trials return whether the difference in 
%                                             firing rate between the stationary and motion periods is significant
%       is_motion_vs_motion_significant     - for a cluster and two sets of trials, return whether the difference in 
%                                             firing rate between the motion periods of the two trial sets is different
%       stationary_fr_for_trial_group       - returns firing rates during stationary periods for a set of trials
%       motion_fr_for_trial_group           - returns firing rates during motion periods for a set of trials
%       check_trial_group                   - checks whether the trial group label appears in the experiment
%       is_mismatch_significant             - for a cluster and a set of mismatch trials, return whether the firing rate 
%                                             after the mismatch is significant
%       get_mismatch_response               - for a cluster and a set of mismatch trials, return the delta firing rate 
%                                             between baseline and response
%       get_mismatch_onset_times            - for a set of mismatch trials, return the time (in probe time) at which 
%                                             the mismatch onset occurs
%       get_motion_onset_times              - for a set of trials, return the time (in probe time) at which the motion periods begin.
%       get_traces_around_mismatch_onset    - return an array of velocity traces around mismatch onset for a set of mismatch trials
%       get_traces_around_solenoid_up       - return an array of velocity traces around the onset of the solenoid 
%                                             going up (i.e. blocking the treadmill)
%       get_traces_around_motion_onset      - return an array of velocity traces around motion onset for a set of trials
%       get_fr_responses                    - return convolved firing rate traces around a set of events
%       create_svm_table                    - creates MATLAB table with, for each cluster and each trial, the firing rate 
%                                             during the stationary and motion periods during that trial.
%       create_replay_offsets_table         - creates MATLAB table with, for each replay trial the amount of offset 
%                                             (in sample points) to take so that taking from that sample point will 
%                                             align it to the original trial
%       create_tuning_curves                - creates cell array of structures with information about velocity tuning 
%                                             of clusters to a set of trials.
%       create_tuning_curves_acceleration   - creates cell array of structures with information about acceleration tuning 
%                                             of clusters to a set of trials.
%       get_traces_around_times             - return an array of velocity traces around a specified set of events
%       get_fr_around_times                 - return convolved firing rate traces around a set of events
%
%       timebase                            - returns a timebase running between limits sampled at particular frequency
%       time_idx_around_trigger_time        - mask trial timebase around an event
%       
%
%   See also: Cluster, Anatomy, Session, RC2Format
%
%   TODO:       1. improve ordering of the methods

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
        % FormattedData
        %
        %   FormattedData(CTL, DATA) takes an object, CTL, of class
        %   RC2Analysis and a structure DATA loaded from the formatted data
        %   file and creates an object to interact with the data.
        %
        %   See also: RC2Format
        
            obj.ctl = ctl;
            obj.data = data;
            
            obj.probe_id = obj.data.probe_id;
            
            % create cluster object for each cluster
            for ii = 1 : length(obj.data.clusters)
                obj.clusters(ii) = Cluster(obj.data.clusters(ii));
            end
            
            for ii = 1 : length(obj.data.anatomy)
                obj.anatomy{ii} = Anatomy(obj.data.anatomy(ii));
            end
            
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
                    case 'mismatch_darkness_jan2022'
                        obj.sessions{ii} = MismatchDarknessJan2022Session(obj.data.sessions(ii));
                    case 'sparse_noise'
                        obj.sessions{ii} = SparseNoiseSession(obj.data.sessions(ii));
                    case 'sf_tf'
                        obj.sessions{ii} = SparseNoiseSession(obj.data.sessions(ii));
                    case 'passive_protocol_always_vis'
                        obj.sessions{ii} = PassiveProtocolAlwaysVisSession(obj.data.sessions(ii));
                    case 'passive_protocol_motion_clouds'
                        obj.sessions{ii} = PassiveProtocolAlwaysVisSession(obj.data.sessions(ii));
                    case 'training'
                        obj.sessions{ii} = TrainingRunningSession(obj.data.sessions(ii));
                end
            end
            
            obj.svm_table = obj.ctl.load_svm_table(obj.data.probe_id);
            obj.replay_offsets = obj.ctl.load_offsets_table(obj.data.probe_id);
            
            if isempty(obj.replay_offsets)
                return
            end
            
            obj.apply_offsets();
        end
        
        
        
        function val = list_trial_group_labels(obj)
        %%list_trial_group_labels Returns a list of the trial group labels
        %%available for this data
        %
        %   TRIAL_GROUP_LABELS = list_trial_group_labels()
        %   returns a cell array of strings with the available trial group
        %   labels for this experiment
        
            sessions = obj.motion_sessions();
            trial_group_labels = cellfun(@(x)(x.trial_group_labels), sessions, 'uniformoutput', false);
            trial_group_labels = [trial_group_labels{:}];
            val = unique(trial_group_labels);
        end
        
        
        function tuning_curve = load_tuning_curves(obj, cluster_id, trial_group)
        %%load_tuning_curves Loads previously created .mat files with the
        %%tuning curve data on speed
        %
        %   TUNING_CURVE_STRUCT = load_tuning_curve(CLUSTER_ID,
        %   TRIAL_GROUP) loads data from a previously saved .mat file.
        %   Loads the data for a single cluster with integer ID,
        %   CLUSTER_ID, and run on trials with trial group, TRIAL_GROUP, a
        %   cell array of trial group labels.
        %
        %   TRIAL_GROUP should be the same as the TRIAL_GROUP supplied to
        %   `create_tuning_curves` during creation of the tuning curve
        %   files.
        %   
        %   TUNING_CURVE_STRUCT is a structure with the tuning curve data
        %   contained in the tuning curve file.
        %
        %   See also: create_tuning_curves, RC2Analysis.create_tuning_curves
        
            tbl = obj.ctl.load_tuning_curves(obj.probe_id);
            group_idx = cellfun(@(x)(isequal(x, trial_group)), tbl.trial_groups);
            tuning_curves = tbl.tuning_curves{group_idx};
            cluster_idx = [tuning_curves(:).cluster_id] == cluster_id;
            tuning_curve = tuning_curves(cluster_idx);
        end
        
        function tuning_curve = load_tuning_curves_acceleration(obj, cluster_id, trial_group)
        %%load_tuning_curves_acceleration Loads previously created .mat files with the
        %%tuning curve data on acceleration
        %
        %   TUNING_CURVE_STRUCT = load_tuning_curve(CLUSTER_ID,
        %   TRIAL_GROUP) loads data from a previously saved .mat file.
        %   Loads the data for a single cluster with integer ID,
        %   CLUSTER_ID, and run on trials with trial group, TRIAL_GROUP, a
        %   cell array of trial group labels.
        %
        %   TRIAL_GROUP should be the same as the TRIAL_GROUP supplied to
        %   `create_tuning_curves` during creation of the tuning curve
        %   files.
        %   
        %   TUNING_CURVE_STRUCT is a structure with the tuning curve data
        %   contained in the tuning curve file.

            tuning_curve = {};
            for i_table = 1 : 3
                tbl = obj.ctl.load_tuning_curves_acceleration(obj.probe_id, i_table);
                group_idx = cellfun(@(x)(isequal(x, trial_group)), tbl.trial_groups);
                tbl = tbl.tuning_curves{group_idx};
                
                cluster_idx = [tbl(:).cluster_id] == cluster_id;
                tuning_curve{i_table} = tbl(cluster_idx);
            end
        end
        
        
        function apply_offsets(obj)
        %%apply_offsets Applies the offsets saved in the offsets .csv files
        %%to the trials.
        %
        %   apply_offsets() takes the table of sample offsets for the
        %   replay trials (relative to their original trial) and applies
        %   those offsets to the `replay_offset` property of the corresponding
        %   Trial objects.
        %
        %   See also: Trial, RC2Format.create_replay_offsets_table
        
            sessions = obj.motion_sessions();
            for ii = 1 : length(sessions)
                sessions{ii}.apply_offsets(obj.replay_offsets);
            end
        end
        
        
        
        function session = get_session_with_id(obj, session_id)
        %%get_session_with_id Returns a session object with a session ID
        %
        %   SESSION = get_session_with_id(SESSION_ID)
        %   recovers the session with string ID, SESSION_ID, (e.g. of form
        %   CA_176_1_rec1_001) and returns the corresponding Session
        %   object.
        %
        %   See also: Session
        
            session_ids = cellfun(@(x)(x.session_id), obj.sessions, 'UniformOutput', false);
            idx = strcmp(session_ids, session_id);
            session = obj.sessions{idx};
        end
        
        
        
        function sessions = motion_sessions(obj)
        %%motion_sessions Returns a cell array of Session objects which are
        %%part of an RVTSession
        %
        %   SESSIONS = motion_session() recovers the sessions which are of
        %   type "RVTSession" and returns them as a cell array of Session
        %   objects in SESSIONS.
        % 
        %   TODO: this is suboptimal atm... would have to change every time
        %   experiment is done and different session order occurs
        
            if strcmp(obj.probe_id, 'CAA-1112872_rec1_rec1b_rec2_rec3') |  strcmp(obj.probe_id, 'CAA-1121416_rec1_rec2')...
                    |  strcmp(obj.probe_id, 'CAA-1121763_rec1_rec2')
                sessions = obj.sessions(1:2);
            elseif strcmp(obj.probe_id, 'CAA-1123304_rec1_rec2_rec3_rec4')
                sessions = obj.sessions(1:4);
            else
                sessions = obj.sessions(1);
            end
        end
        
        
        
        function trials = motion_trials(obj)
        %%motion_trials Returns a cell array of Trial objects which are
        %%part of an RVTSession
        %
        %   TRIALS = motion_trials() recovers the trials which are part of
        %   an "RVTSession" and returns them as a cell array of Tession
        %   objects in TRIALS.
        % 
        %   See also: motion_sessions
            
            sessions = obj.motion_sessions();
            t = cellfun(@(x)(x.trials), sessions, 'UniformOutput', false);
            trials = [t{:}]; %#ok<*PROP>
        end
        
        
        
        function trials = get_trials_with_trial_group_label(obj, trial_group_label)
        %%get_trials_with_trial_group_label Returns a cell array of Trials
        %%objects with a specified `trial_group_label`
        %
        %   TRIALS = get_trials_with_trial_group_label(TRIAL_GROUP_LABELS)
        %   takes a string, or cell array of strings, in
        %   TRIAL_GROUP_LABELS, specifying the types of trial to return. If
        %   TRIAL_GROUP_LABELS is a string, all trials which have
        %   `trial_group_label` matching the string are returned as a cell
        %   array of Trial objects. If TRIAL_GROUP_LABELS is a cell array
        %   of strings, all trials which  match one of the strings are
        %   returned as a cell array of Trial objects.
        %
        %   See also: RVTSession, get_trials_with_trial_ids
            
            sessions = obj.motion_sessions(); %#ok<*PROPLC>
            
            trials = {};
            for ii = 1 : length(sessions)
                these_trials = sessions{ii}.get_trials_with_trial_group_label(trial_group_label);
                trials = [trials, these_trials]; %#ok<*AGROW>
            end
        end
        
        
        
        function trials = get_trials_with_trial_ids(obj, trial_ids)
        %%get_trials_with_trial_ids Returns the trial or trials with the
        %%specified trial IDs
        %
        %   TRIALS = get_trials_with_trial_ids(TRIAL_GROUP_IDS)
        %   takes an integer, or vector of integers, in
        %   TRIAL_GROUP_IDS, specifying the types of trial to return. If
        %   TRIAL_GROUP_LABELS is a string, all trials which have
        %   `trial_group_label` matching the string are returned as a cell
        %   array of Trial objects. If TRIAL_GROUP_LABELS is a cell array
        %   of strings, all trials which  match one of the strings are
        %   returned as a cell array of Trial objects.
        %
        %   See also: RVTSession, get_trials_with_trial_ids
        
            trials = obj.motion_trials();
            idx = cellfun(@(x)(ismember(x.trial_id, trial_ids)), trials);
            trials = trials(idx);
            
            if length(trials) == 1
                trials = trials{1};
            end
        end
        
        
        
        function cluster = get_cluster_with_id(obj, cluster_id)
        %%get_cluster_with_id Returns the Cluster object for the cluster
        %%with ID
        %
        %   CLUSTER = get_cluster_with_id(CLUSTER_ID)
        %   returns an object of class Cluster, for cluster with integer
        %   ID, CLUSTER_ID.
        %
        %   See also: Cluster
        
            idx = ismember([obj.clusters(:).id], cluster_id);
            cluster = obj.clusters(idx);
        end
        
        
        
        function selected_clusters = selected_clusters(obj, spiking_class)
        %%selected_clusters Returns an object array of class Clusters with
        %%the clusters selected on the manual sorting step
        %
        %   CLUSTERS = selected_clusters(SPIKE_CLASS) returns an object
        %   array of class Clusters in CLUSTERS. SPIKE_CLASS is optional,
        %   and if supplied should be a string 'any' (default), 'narrow',
        %   or 'wide', which restricts the Cluster object array to only
        %   narrow or wide spike clusters ('any' returns all selected
        %   clusters).
        %
        %   See also: selected_cluster_ids, VISp_cluster_ids
        
            VariableDefault('spiking_class', 'any');
            
            search_ids = intersect(obj.data.selected_clusters, obj.spiking_class_ids(spiking_class));
            idx = ismember([obj.clusters(:).id], search_ids);
            selected_clusters = obj.clusters(idx);
        end
        
        
        
        function ids = selected_cluster_ids(obj, spiking_class)
        %%selected_cluster_ids Returns a vector of integer IDs
        %%corresponding to the clusters selected on the manual sorting step
        %
        %   CLUSTER_IDS = selected_cluster_ids(SPIKE_CLASS) returns a
        %   vector of integer IDS in CLUSTER_IDS. SPIKE_CLASS is optional,
        %   and if supplied should be a string 'any' (default), 'narrow',
        %   or 'wide', which restricts the Cluster object array to only
        %   narrow or wide spike clusters ('any' returns all selected
        %   clusters).
        %
        %   See also: selected_clusters
        
            VariableDefault('spiking_class', 'any');
            
            selected_clusters = obj.selected_clusters(spiking_class);
            ids = [selected_clusters(:).id];
        end
        
        
        
        function selected_clusters = selected_mua_clusters(obj, region_acronym)
        %%selected_mua_clusters Returns an object array of class Clusters
        %%with the MUA clusters selected on the manual sorting step
        %
        %   CLUSTERS = selected_mua_clusters(REGION_ACRONYM) returns an
        %   object array of class Clusters in CLUSTERS. REGION_ACRONYM is 
        %   optional, and if supplied should be a string 'any' (default),
        %   or a valid brain region acronym.
        %
        %   Only a subset of acronyms are allowed, and specified in the
        %   `region_mapping.csv` file (see `is_region_subregion_of`).
        %
        %   If supplied only clusters allocated to that region will be
        %   returned in the array.
        %
        %   See also: selected_clusters, is_region_subregion_of
        
            VariableDefault('region_acronym', 'any');
            
            idx = ismember([obj.clusters(:).id], obj.data.selected_mua_clusters);
            
            if strcmp(region_acronym, 'any')
                selected_clusters = obj.clusters(idx);
                return
            end
            
            % regions of all clusters
            all_regions = {obj.clusters(:).region_str};
            
            % index of those clusters which are in the specified region
            idx = idx(:) & is_region_subregion_of(region_acronym, all_regions);
            
            selected_clusters = obj.clusters(idx);
        end
        
        
        
        function visp_clusters = VISp_clusters(obj, is_selected, spiking_class)
        %%VISp_clusters Returns an object array of class Clusters with the
        %%clusters which are allocated to a VISp layer
        %
        %   CLUSTERS = VISp_clusters(MANUALLY_SELECTED, SPIKE_CLASS)
        %   returns an object array of class Clusters in CLUSTERS.
        %   MANUALLY_SELECTED is optional and if supplied should be a
        %   boolean indicating whether to further restrict the VISp
        %   clusters to those which were manually selected in Phys (true,
        %   default) or all VISp clusters sorted by kilosort (false).
        %   SPIKE_CLASS is optional, and if supplied should be a string
        %   'any' (default), 'narrow', or 'wide', which restricts the
        %   Cluster object array to only narrow or wide spike clusters
        %   ('any' returns all selected clusters).
        %
        %   See also: selected_clusters, VISp_cluster_ids
        
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
        %%VISp_cluster_ids Returns a vector of integer IDs corresponding to
        %%the clusters selected on the manual sorting step and which are
        %%allocated to a VISp layers
        %
        %   CLUSTER_IDS = VISp_cluster_ids(MANUALLY_SELECTED, SPIKE_CLASS)
        %   returns a vector of integer IDS in CLUSTER_IDS.
        %   MANUALLY_SELECTED is optional and if supplied should be a 
        %   boolean indicating whether to further restrict the VISp
        %   clusters to those which were manually selected in Phys (true,
        %   default) or all VISp clusters sorted by kilosort (false).
        %   SPIKE_CLASS is optional, and if supplied should be a string
        %   'any' (default), 'narrow', or 'wide', which restricts the
        %   Cluster object array to only narrow or wide spike clusters
        %   ('any' returns all selected clusters). 
        %
        %   See also: selected_clusters, VISp_clusters
        
            VariableDefault('is_selected', true);
            VariableDefault('spiking_class', 'any');
            
            visp_clusters = obj.VISp_clusters(is_selected, spiking_class);
            ids = [visp_clusters(:).id];
        end
        
        
        
        function ids = spiking_class_ids(obj, spiking_class)
        %%spiking_class_ids Returns a vector of integer IDs corresponding to
        %%the clusters with a particular spiking class
        %
        %   CLUSTER_IDS = spiking_class_ids(SPIKE_CLASS)
        %   returns a vector of integer IDS in CLUSTER_IDS.
        %   SPIKE_CLASS is optional, and if supplied should be a string
        %   'any' (default), 'narrow', or 'wide', which restricts the
        %   Cluster object array to only narrow or wide spike clusters
        %   ('any' returns all selected clusters). 
        %
        %   NOTE: All clusters sorted by Kilosort and satisfying
        %   the SPIKE_CLASS option are returned and not just the manually
        %   selected clusters after Phy. 
        %
        %   See also: selected_clusters, VISp_clusters
            
            VariableDefault('spiking_class', 'any');
        
            ids = [obj.clusters(:).id];
            
            if ~strcmp(spiking_class, 'any')
                all_spiking_classes = {obj.clusters(:).spiking_class};
                idx = strcmp(all_spiking_classes, spiking_class);
                ids = ids(idx);
            end
        end
        
        
        
        function [relative_depth, layer] = get_relative_layer_depth_of_cluster(obj, cluster_id)
        %%get_relative_layer_depth_of_cluster Returns the relative depth of
        %%a cluster in a VISp layer
        %
        %   [RELATIVE_DEPTH, LAYER] = get_relative_layer_depth_of_cluster(CLUSTER_ID)
        %   returns the relative depth (from upper boundary) of a cluster
        %   with integer ID CLUSTER_ID within the VISp layer to which it
        %   has been assigned. RELATIVE_DEPTH is a value between 0 and 1
        %   with 0 indicating that the cluster is at the upper boundary of
        %   the layer and 1 indicating that the cluster is at the lower
        %   boundary of the layer. LAYER is a string containing which VISp
        %   layer the cluster was assigned to (e.g. 'VISp4').
        %
        %   If the cluster is not in VISp, RELATIVE_DEPTH is -1 and LAYER
        %   is an empty string.
            
            cluster = obj.get_cluster_with_id(cluster_id);
            if ~cluster.is_VISp; relative_depth = -1; layer = ''; return; end
            
            shank_id = cluster.cluster.shank_id;
            idx = cellfun(@(x)(x.shank_id), obj.anatomy) == shank_id;
            [relative_depth, layer] = obj.anatomy{idx}.VISp_layer_relative_depth(cluster.distance_from_probe_tip);
        end
        
        
        
        function bouts = get_motion_bouts_for_trial_group(obj, trial_group_label, options)
        %%get_motion_bouts_for_trial_group Returns a cell array of
        %%MotionBout objects extracted from the trials of a particular type
        %
        %   BOUTS_ARRAY = get_motion_bouts_for_trial_group(GROUP_LABELS, OPTIONS)
        %
        %   Returns all bouts in the trials in GROUP_LABELS. GROUP_LABELS
        %   can be a single string with a group label, or a cell array of
        %   group labels.
        %
        %   OPTIONS is used to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a motion bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low
        %                               to look for motion onset
        %
        %   BOUTS_ARRAY is a cell array of MotionBout objects.
        %
        %   See also: MotionBout, Trial.motion_bouts, RVTSession.get_motion_bouts_for_trial_group
        %
        %   TODO:   1. OPTIONS should become optional
        
            sessions = obj.motion_sessions(); %#ok<*PROPLC>
            
            bouts = {};
            for ii = 1 : length(sessions)
                these_bouts = sessions{ii}.get_motion_bouts_for_trial_group(trial_group_label, options);
                bouts = [bouts, these_bouts]; %#ok<*AGROW>
            end
        end
        
        
        
        function times = get_motion_bout_start_times(obj, trial_group_label)
        %TODO: Unused, remove?  This wouldn't work anyway
        
            bouts = obj.get_motion_bouts_for_trial_group(trial_group_label);
            times = cellfun(@(x)(x.start_time), bouts);
        end
        
        
        
        function [sig, p, direction, stat_med, mot_med, n] = is_stationary_vs_motion_significant(obj, cluster_id, trial_group_label)
        %%is_stationary_vs_motion_significant For a cluster and a set of
        %%trials return whether the difference in firing rate between the
        %%stationary and motion periods is significant
        %
        %   [SIGNIFICANT, P_VALUE, DIRECTION, STATIONARY_MEDIAN, MOTION_MEDIAN, N] = 
        %       is_stationary_vs_motion_significant(CLUSTER_ID, TRIAL_GROUP_LABEL)
        %
        %   For the cluster with integer ID CLUSTER_ID, and the set of
        %   trials specified with TRIAL_GROUP_LABEL, gather the firing
        %   rates on those trials during the stationary and motion periods
        %   and perform a Wilcoxon sign-rank test stationary vs. motion.
        %
        %   Outputs:
        %    SIGNIFICANT - whether the test is significant 1 or not 0.
        %                  significance is determined by the P_VALUE < 0.05
        %    P_VALUE     - p-value of the Wilcoxon sign-rank test
        %    DIRECTION   - either 1, 0, or -1
        %                   1  - motion firing rates are > stationary firing rates
        %                   0  - no significant difference between firing rates
        %                   -1 - motion firing rates are < stationary firing rates
        %   STATIONARY_MEDIAN - median firing rates across the stationary periods
        %   MOTION_MEDIAN     - median firing rates across the motion periods
        %   N            - the number of trials used for the test
        %
        %   See also: compare_groups_with_signrank
        %
        %   Note:   Trials with NaN firing rate values in either the
        %           stationary or motion period are removed before the test
        %           is performed.
        %           If the TRIAL_GROUP_LABEL is not valid (i.e. no trials
        %           are found), then all outputs are empty      
        
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
        %%is_motion_vs_motion_significant For a cluster and two sets of
        %%trials, return whether the difference in firing rate between the
        %%motion periods of the two trial sets is different
        %
        %   [SIGNIFICANT, P_VALUE, DIRECTION, X_MEDIAN, Y_MEDIAN, N] = 
        %       is_stationary_vs_motion_significant(CLUSTER_ID, TRIAL_GROUP_LABEL_X, TRIAL_GROUP_LABEL_Y, MAX_N_TRIALS)
        %
        %   For the cluster with integer ID CLUSTER_ID, and the two sets of
        %   trials specified with TRIAL_GROUP_LABEL_X and
        %   TRIAL_GROUP_LABEL_Y, gather the firing
        %   during the motion periods for the two sets of trials
        %   and perform a Wilcoxon sign-rank test stationary vs. motion.
        %
        %
        %   Outputs:
        %    SIGNIFICANT - whether the test is significant 1 or not 0.
        %                  significance is determined by the P_VALUE < 0.05
        %    P_VALUE     - p-value of the Wilcoxon sign-rank test
        %    DIRECTION   - either 1, 0, or -1
        %                   1  - motion firing rates in trial set X are >
        %                   motion firing rates in trial set Y
        %                   0  - no significant difference between firing rates
        %                   -1 - motion firing rates in trial set X are <
        %                   motion firing rates in trial set Y
        %   X_MEDIAN     - median firing rates across the motion periods in
        %                  trial set X
        %   Y_MEDIAN     - median firing rates across the motion periods in
        %                  trial set Y
        %   N            - the number of trials used for the test
        %
        %   See also: compare_groups_with_signrank
        %
        %   Note:   Performing the sign-rank requires pairing the trials in
        %           some way. As the trials were performed in batches (i.e.
        %           all trial types were performed once, before another
        %           batch of trials was performed), we pair trials in the
        %           same batch.
        %           
        %           Paired trials with NaN firing rate values in either the X or Y
        %           trial are removed before the test is performed.
        %
        %           MAX_N_TRIALS is optional and if set selects only the
        %           first MAX_N_TRIALS pairs of trials for the test.
        %   
        %           The order of TRIAL_GROUP_X and TRIAL_GROUP_Y determines
        %           the DIRECTION, X_MEDIAN and Y_MEDIAN outputs.
        %
        %           If either TRIAL_GROUP_LABEL_X or TRIAL_GROUP_LABEL_Y is not valid (i.e. no trials
        %           are found), then all outputs are empty
        
        
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
        %%stationary_fr_for_trial_group Returns firing rates during
        %%stationary periods for a set of trials
        %
        %   [FIRING_RATES, TRIAL_IDS] = stationary_fr_for_trial_group(CLUSTER_ID, TRIAL_GROUP_LABEL)
        %   gathers together the firing rates for cluster with integer ID
        %   CLUSTER_ID, on the trials specified with TRIAL_GROUP_LABEL,
        %   during the stationary periods for those trials.
        %
        %   Outputs:
        %    FIRING_RATES   - #trials x 1 vector of firing rates during stationary periods
        %    TRIAL_IDS      - #trials x 1 vector of trial IDs corresponding
        %                     to the firing rates
        %
        %   Note:  firing rates are taken from the pre-saved table of
        %           stationary/motion firing rates (see README)
        %
        %   See also:   create_svm_table, RC2Format.create_svm_table
            
            valid = obj.check_trial_group(trial_group_label);
            if ~valid; fr = []; return; end
        
            idx = ismember(obj.svm_table.trial_group_label, trial_group_label) & ...
                  obj.svm_table.cluster_id == cluster_id;
            
            fr = obj.svm_table.stationary_fr(idx);
            trial_id = obj.svm_table.trial_id(idx);
        end
        
        
        
        function [fr, trial_id] = motion_fr_for_trial_group(obj, cluster_id, trial_group_label)
        %%motion_fr_for_trial_group Returns firing rates during
        %%motion periods for a set of trials
        %
        %   [FIRING_RATES, TRIAL_IDS] = motion_fr_for_trial_group(CLUSTER_ID, TRIAL_GROUP_LABEL)
        %   gathers together the firing rates for cluster with integer ID
        %   CLUSTER_ID, on the trials specified with TRIAL_GROUP_LABEL,
        %   during the motion periods for those trials.
        %
        %   Outputs:
        %    FIRING_RATES   - #trials x 1 vector of firing rates during motion periods
        %    TRIAL_IDS      - #trials x 1 vector of trial IDs corresponding
        %                     to the firing rates
        %
        %   Note:  firing rates are taken from the pre-saved table of
        %           stationary/motion firing rates (see README)
        %
        %   See also:   create_svm_table, RC2Format.create_svm_table
        
            valid = obj.check_trial_group(trial_group_label);
            if ~valid; fr = []; return; end
            
            idx = ismember(obj.svm_table.trial_group_label, trial_group_label) & ...
                  obj.svm_table.cluster_id == cluster_id;
            
            fr = obj.svm_table.motion_fr(idx);
            trial_id = obj.svm_table.trial_id(idx);
        end
        
        
        
        function valid = check_trial_group(obj, trial_group_label)
        %%check_trial_group Checks whether the trial group label appears in
        %%the experiment
        %
        %   EXIST = check_trial_group(TRIAL_GROUP_LABEL)
        %   checks whether the labels in TRIAL_GROUP_LABEL are part of the
        %   current experiment.
        %   TRIAL_GROUP_LABEL can be a char vector or cell array of
        %   char vectors with the label. If it is a char vector EXIST will
        %   be 0 (label doesn't exist in experiment) or 1 (label does
        %   exist). If it is a Nx1 cell array of char vectors, EXIST will
        %   be a Nx1 boolean array indicating whether each label exists in
        %   the experiment.
        
            sessions = obj.motion_sessions();
            t = cellfun(@(x)(x.trial_group_labels), sessions, 'uniformoutput', false);
            valid = ismember(trial_group_label, [t{:}]);
        end
        
        
        
        function [sig, p, direction] = is_mismatch_significant(obj, cluster_id, trial_group_label)
        %%is_mismatch_significant For a cluster and a set of mismatch
        %%trials, return whether the firing rate after the mismatch is
        %%significant
        %
        %   [SIGNIFICANT, P_VALUE, DIRECTION] = is_mismatch_significant(CLUSTER_ID, TRIAL_GROUP_LABEL)
        %
        %   For the cluster with integer ID CLUSTER_ID, and the set of
        %   mismatch trials specified with TRIAL_GROUP_LABEL, test whether
        %   the firing rate change after mismatch onset, relative to
        %   baseline, is significant. 
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
        %   If the specified TRIAL_GROUP_LABEL doesn't exist in the
        %   experiment, all outputs are empty.
        %
        %   See also: MismatchAnalysis, MismatchAnalysis.is_response_significant
        
            valid = obj.check_trial_group(trial_group_label);
            if ~valid; sig = []; p = []; direction = []; return; end
            
            mm = MismatchAnalysis();
            
            cluster = obj.get_cluster_with_id(cluster_id);
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            [sig, p, direction] = mm.is_response_significant(cluster, trials);
        end
        
        
        
        function delta_fr = get_mismatch_response(obj, cluster_id, trial_group_label)
        %%get_mismatch_response For a cluster and a set of mismatch
        %%trials, return the delta firing rate between baseline and
        %%response.
        %
        %   DELTA_FIRING_RATE = get_mismatch_response(CLUSTER_ID, TRIAL_GROUP_LABEL)
        %
        %   For the cluster with integer ID CLUSTER_ID, and the set of
        %   mismatch trials specified with TRIAL_GROUP_LABEL, return the
        %   average delta firing rate in Hz (response FR - baseline FR), in
        %   DELTA_FIRING_RATE. 
        %
        %   If the specified TRIAL_GROUP_LABEL doesn't exist in the
        %   experiment, the output is empty.
        %
        %   See also: MismatchAnalysis, MismatchAnalysis.is_response_significant
        
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
        %%get_mismatch_onset_times For a set of mismatch
        %%trials, return the time (in probe time) at which the mismatch
        %%onset occurs.
        %
        %   TIMES = get_mismatch_onset_times(TRIAL_GROUP_LABEL)
        %
        %   For the set of mismatch trials specified with
        %   TRIAL_GROUP_LABEL, return the times (in probe time) at which
        %   the mismatch onset occurs. Return times as TIMES, a #trials x 1
        %   vector with the times in.
        %
        %   If the specified TRIAL_GROUP_LABEL doesn't exist in the
        %   experiment, the output is empty.
        
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            times = cellfun(@(x)(x.mismatch_onset_t), trials);
        end
        
        
        
        function [times, bouts] = get_motion_onset_times(obj, trial_group_label, options)
        %%get_motion_onset_times For a set of trials, return the time (in
        %%probe time) at which the motion periods begin.
        %
        %   [TIMES, BOUTS_ARRAY] = get_mismatch_onset_times(TRIAL_GROUP_LABEL, OPTIONS)
        %
        %   For the set of trials specified with TRIAL_GROUP_LABEL, return
        %   the times (in probe time) at which the motion onset occurs.
        %   Return times as TIMES, a #motion periods x 1 vector with the
        %   times in.  Also returned is BOUTS_ARRAY, a cell array of MotionBout objects
        %   corresponding to each element of TIME.
        %
        %   If the specified TRIAL_GROUP_LABEL doesn't exist in the
        %   experiment, the output is empty.
        %
        %   OPTIONS is used to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a motion bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low 
        %                               to look for motion onset
        %
        %   See also: MotionBout, get_motion_bouts_for_trial_group
        
            bouts = obj.get_motion_bouts_for_trial_group(trial_group_label, options);
            times = cellfun(@(x)(x.start_time), bouts);
        end
        
        
        
        function [traces, t, times] = get_traces_around_mismatch_onset(obj, trial_group_label, limits, common_fs)
        %%get_traces_around_mismatch_onset Return an array of velocity
        %%traces around mismatch onset for a set of mismatch trials
        %
        %   [TRACES, T, TIMES] =
        %   get_traces_around_mismatch_onset(TRIAL_GROUP_LABEL, LIMITS, FS)
        %   for all (mismatch) trials specified by TRIAL_GROUP_LABEL get
        %   velocity traces around mismatch onset. The amount of data
        %   either side of the mismatch onset is specified by LIMITS, a 1x2
        %   vector of the form [start_time, end_time], giving the time in
        %   seconds relative to the mismatch onset to take (e.g. [-1, 1]
        %   would get the trace 1 second before to 1 second after mismatch
        %   onset. FS specifies the sampling rate (in Hz) at which to
        %   sample the traces.
        %
        %   Outputs:
        %    TRACES     - #samples x #trials matrix with the velocity trace
        %                 around mismatch onset along the columns (#samples
        %                 is dependent on FS and LIMITS).
        %    T          - #samples x 1 vector with the common time base for
        %                 all the traces
        %   TIMES       - #trials x 1 vector with the time (in probe time)
        %                 at which the mismatch onset occurs (same outputs
        %                 as `get_mismatch_onset_times`)
        %   
        %   See also: get_traces_around_times, get_mismatch_onset_times
        
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            times = obj.get_mismatch_onset_times(trial_group_label);
            [traces, t] = obj.get_traces_around_times(trials, times, limits, common_fs);
        end
        
        
        
        function [traces, t, times] = get_traces_around_solenoid_up(obj, trial_group_label, limits, common_fs)
        %%get_traces_around_solenoid_up Return an array of velocity
        %%traces around the onset of the solenoid going up (i.e. blocking
        %%the treadmill)
        %
        %   [TRACES, T, TIMES] =
        %   get_traces_around_mismatch_onset(TRIAL_GROUP_LABEL, LIMITS, FS)
        %   for all trials specified by TRIAL_GROUP_LABEL get
        %   velocity traces around onset of solenoid going up (blocking
        %   treadmill). The amount of data 
        %   either side of the onset is specified by LIMITS, a 1x2
        %   vector of the form [start_time, end_time], giving the time in
        %   seconds relative to the onset to take (e.g. [-1, 1]
        %   would get the trace 1 second before to 1 second after solenoid
        %   up). FS specifies the sampling rate (in Hz) at which to
        %   sample the traces.
        %
        %   Outputs:
        %    TRACES     - #samples x #trials matrix with the velocity trace
        %                 around onset along the columns (#samples
        %                 is dependent on FS and LIMITS).
        %    T          - #samples x 1 vector with the common time base for
        %                 all the traces
        %   TIMES       - #trials x 1 vector with the time (in probe time)
        %                 at which the onset occurs
        %   
        %   See also: get_traces_around_times
        
            trials = obj.get_trials_with_trial_group_label(trial_group_label);
            idx = cellfun(@(x)(find(diff(x.solenoid > 2.5) == 1) + 1), trials, 'uniformoutput', false);
            times = cellfun(@(x, y)(x.probe_t(y)), trials, idx);    
            [traces, t] = obj.get_traces_around_times(trials, times, limits, common_fs);
        end
        
        
        
        function [traces, t, times] = get_traces_around_motion_onset(obj, trial_group_label, limits, common_fs, options)
        %%get_traces_around_motion_onset Return an array of velocity
        %%traces around motion onset for a set of trials
        %
        %   [TRACES, T, TIMES] =
        %   get_traces_around_motion_onset(TRIAL_GROUP_LABEL, LIMITS, FS, OPTIONS)
        %   for all trials specified by TRIAL_GROUP_LABEL get
        %   velocity traces around motion period onset. The amount of data
        %   either side of the motion onset is specified by LIMITS, a 1x2
        %   vector of the form [start_time, end_time], giving the time in
        %   seconds relative to the motion onset to take (e.g. [-1, 1]
        %   would get the trace 1 second before to 1 second after motion
        %   onset. FS specifies the sampling rate (in Hz) at which to
        %   sample the traces.
        %
        %   OPTIONS is used to determine which bouts are selected
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a motion bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low 
        %                               to look for motion onset
        %
        %   Outputs:
        %    TRACES     - #samples x #motion bouts matrix with the velocity trace
        %                 around motion onset along the columns (#samples
        %                 is dependent on FS and LIMITS).
        %    T          - #samples x 1 vector with the common time base for
        %                 all the traces
        %   TIMES       - #motion bouts x 1 vector with the time (in probe time)
        %                 at which the motion onset occurs (same outputs
        %                 as `get_mismatch_onset_times`)
        %
        %   See also:   get_motion_onset_times, get_traces_around_times
        
            [times, bouts] = obj.get_motion_onset_times(trial_group_label, options);
            trials = cellfun(@(x)(x.trial), bouts, 'uniformoutput', false);
            [traces, t] = obj.get_traces_around_times(trials, times, limits, common_fs);
        end
        
        
        
        function [fr, t] = get_fr_responses(obj, cluster_id, trial_group_label, type, limits, fs, options)
        %%get_fr_responses Return convolved firing rate traces around a set
        %%of events
        %
        %   [FIRING_RATE, T] = get_fr_responses(CLUSTER_ID, TRIAL_GROUP_LABEL, TYPE, LIMITS, FS, OPTIONS)
        %   for the cluster with integer ID CLUSTER_ID, and all trials
        %   specified by TRIAL_GROUP_LABEL get convolved firing rate traces
        %   around the onset of a set of events, specified with TYPE. TYPE
        %   can be either 'motion' or 'mismatch'. If 'motion' the event is
        %   all motion period onsets in the trials. If 'mismatch' the even
        %   is all mismatch onsets in the trials. The amount of data 
        %   either side of the onset is specified by LIMITS, a 1x2
        %   vector of the form [start_time, end_time], giving the time in
        %   seconds relative to the onset to take (e.g. [-1, 1]
        %   would get the trace 1 second before to 1 second after solenoid
        %   up). FS specifies the sampling rate (in Hz) at which to
        %   sample the traces.
        %
        %   OPTIONS is used to determine which bouts are selected (if TYPE
        %   = 'motion').
        %       it is a structure with fields
        %               min_bout_duration -  the minimum duration in s for
        %                                   a motion bout to be included
        %               include_200ms - whether or not to include the first
        %                               200ms after the solenoid goes low
        %                               to look for motion onset
        %
        %   See also: get_fr_around_times
        
            if strcmp(type, 'motion')
                times = obj.get_motion_onset_times(trial_group_label, options);
            elseif strcmp(type, 'mismatch')
                times = obj.get_mismatch_onset_times(trial_group_label);
            end
            
            [fr, t] = obj.get_fr_around_times(cluster_id, times, limits, fs);
        end
        
        
        
        function tbl = create_svm_table(obj)
        %%create_svm_table Creates MATLAB table with, for each cluster and
        %%each trial, the firing rate during the stationary and motion
        %%periods during that trial.
        %
        %   TABLE = create_svm_table() returns a TABLE with a row for each
        %   selected cluster (returned with `selected_clusters`) and each
        %   trial (all trials in the experiment) (total of
        %   #clusters*#trials rows). Each row has the following
        %   information:
        %   
        %       probe_id - probe recording ID associated with trial
        %       session_id - session ID associated with trial
        %       trial_id - the integer ID of the trial
        %       trial_group_label - the group label associated with the trial
        %       cluster_id - the integer ID of the cluster
        %       stationary_fr - the firing rate of the cluster duing the
        %                       stationary periods of the trial
        %       stationary_time - the amount of time in seconds the trial
        %                         was stationary
        %       motion_fr - the firing rate of the cluster during the
        %                   motion periods of the trial
        %       motion_time - the amoutn of time in seconds the trial was
        %                     in motion
        %
        %   See also: RC2Format.create_svm_table
        
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
        %%create_replay_offsets_table Creates MATLAB table with, for each
        %%replay trial the amount of offset (in sample points) to take so
        %%that taking from that sample point will align it to the original
        %%trial
        %
        %   TABLE = create_replay_offsets_table() returns a TABLE with a
        %   row for each replay trial. Each row has the following
        %   information:
        %   
        %       original_trial_id - the integer ID of the original trial
        %       replay_trial_id - the integer ID of the replay trial
        %       sample_offset - the sample number  which will align the
        %                       replay to the original
        %
        %   See also: RC2Format.create_replay_offsets_table
        
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
        %%create_tuning_curves  Creates cell array of structures with
        %%information about velocity tuning of clusters to a set of trials.
        %
        %   TUNING_INFO = create_tuning_curves(TRIAL_TYPES)
        %   takes a cell array of trial group labels. So TRIAL_TYPES is of
        %   the form {TRIAL_GROUP_LABEL_1, TRIAL_GROUP_LABEL_2, ...}. Each
        %   TRIAL_GROUP_LABEL_N is itself either a char vector indicating a
        %   trial group label or a cell array of char vectors indicating a
        %   set of trial group labels.
        %
        %   In either case, the code loops over entries of TRIAL_TYPES and
        %   selects the trials specified by the TRIAL_GROUP_LABEL_N labels. 
        %
        %   These trials are then passed to the TuningTable object to
        %   create the velocity bins and extract the information about
        %   tuning for each cluster. 
        %
        %   TUNING_INFO is thus a cell array of the same length as
        %   TRIAL_TYPES, and each entry is a structure array (of length
        %   equal to the number of selected clusters returned by
        %   `selected_clusters`), with the tuning curve info returned by
        %   TuningTable.tuning_curve.
        %
        %   See also: TuningTable, TuningTable.tuning_curve, RC2Analysis.create_tuning_curves
        
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
        
        function tbl = create_tuning_curves_acceleration(obj, trial_types)
        %%create_tuning_curves_acceleration  Creates cell array of structures with
        %%information about acceleration tuning of clusters to a set of trials.
        %
        %   TUNING_INFO = create_tuning_curves_acceleration(TRIAL_TYPES)
        %   takes a cell array of trial group labels. So TRIAL_TYPES is of
        %   the form {TRIAL_GROUP_LABEL_1, TRIAL_GROUP_LABEL_2, ...}. Each
        %   TRIAL_GROUP_LABEL_N is itself either a char vector indicating a
        %   trial group label or a cell array of char vectors indicating a
        %   set of trial group labels.
        %
        %   In either case, the code loops over entries of TRIAL_TYPES and
        %   selects the trials specified by the TRIAL_GROUP_LABEL_N labels. 
        %
        %   These trials are then passed to the TuningTable object to
        %   create the acceleration bins and extract the information about
        %   tuning for each cluster. 
        %
        %   TUNING_INFO is thus a cell array of the same length as
        %   TRIAL_TYPES, and each entry is a structure array (of length
        %   equal to the number of selected clusters returned by
        %   `selected_clusters`), with the tuning curve info returned by
        %   TuningTable.tuning_curve.

            modalities = ["all", "acc", "dec"];
            
            tbl = cell(1, length(modalities), length(trial_types));
            
            for ii = 1 : length(trial_types)
                
                tt = TuningTableAcc(obj.probe_id);
                
                % get all trials of type specified in the cell array
                % 'trial_types'
                trials = obj.get_trials_with_trial_group_label(trial_types{ii});
             
                % if the trials are replays, align them to the original
                for jj = 1 : length(trials)
                    aligned_trials{jj} = trials{jj}.to_aligned;
                end

                for i_mod = 1 : length(modalities)
                    % Use all accelerations, or only positive or negative
                    % values
                    tt.mode = modalities(i_mod);
                    
                    % add these trials to the TuningTable object
                    tt.add_trials(aligned_trials, trial_types{ii});

                    % all selected clusters add them to the tuning table
                    clusters = obj.selected_clusters();
                    for jj = 1 : length(clusters)
                        tbl{i_mod}{ii}(jj) = tt.tuning_curve(clusters(jj));
                    end
                end
            end
        end
       
        
        function [traces, t] = get_traces_around_times(obj, trials, times, limits, fs)
        %%get_traces_around_times Return an array of velocity
        %%traces around a specified set of events
        %
        %   [TRACES, T] = get_traces_around_times(TRIALS, TIMES, LIMITS, FS)
        %   gets velocity traces around a set of events. 
        %   Events are specified in TRIALS and TIMES. TRIALS should be a
        %   #events x 1 cell array with the Trial objects in which each
        %   event occurs (the same Trial object can appear mutliple times
        %   in this cell array if different events occur in the same trial). 
        %   TIMES is a #events x 1 vector specified the time of each event
        %   (in probe time). Each event in TIMES must be matched to the
        %   correct trial in TRIALS (i.e. same position in the cell array
        %   as in the vector).
        %   The amount of data either side of the events is specified by
        %   LIMITS, a 1x2 vector of the form [start_time, end_time], giving
        %   the time in seconds relative to the events to take (e.g. [-1, 1]
        %   would get the trace 1 second before to 1 second after event)
        %   FS specifies the sampling rate (in Hz) at which to
        %   sample the traces.
        %
        %   Outputs:
        %    TRACES     - #samples x #events matrix with the velocity trace
        %                 around events along the columns (#samples
        %                 is dependent on FS and LIMITS).
        %    T          - #samples x 1 vector with the common time base for
        %                 all the traces
        
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
        %%get_fr_around_times Return convolved firing rate traces around a set
        %%of events
        %
        %   [FIRING_RATE, T] = get_fr_around_times(CLUSTER_ID, TIMES, LIMITS, FS)
        %   for the cluster with integer ID CLUSTER_ID, get convolved firing rate traces
        %   around the set of events defined in TIMES. TIMES is a #events x
        %   1 vector specified the time of each event (in probe time). The
        %   amount of data either side of the onset is specified by LIMITS,
        %   a 1x2 vector of the form [start_time, end_time], giving the
        %   time in seconds relative to the onset to take (e.g. [-1, 1]
        %   would get the trace 1 second before to 1 second after solenoid
        %   up). FS specifies the sampling rate (in Hz) at which to
        %   sample the traces.
            
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
        %%timebase Returns a timebase running between limits sampled at
        %%particular frequency
        %
        %   TIMEBASE = timebase(LIMITS, FS)
        %   Creates a timebase running between LIMITS (1x2 vector of form
        %   [time_min, time_max], and sampled at FS Hz.
        
            val = linspace(limits(1), limits(2), round((limits(2) - limits(1))*fs)+1);
        end
        
        
        
        function idx = time_idx_around_trigger_time(trial, trigger_t, limits)
        %%time_idx_around_trigger_time Mask trial timebase around an event
        %
        %   MASK = time_idx_around_trigger_time(TRIAL, TIME, LIMITS)
        %   returns a boolean MASK of same length as the number of sample
        %   points in the TRIAL (object of class Trial), with 1 values
        %   around a TIME (single value in probe time).  The amount of data
        %   either side of the events is specified by 
        %   LIMITS, a 1x2 vector of the form [start_time, end_time], giving
        %   the time in seconds relative to the events to take (e.g. [-1, 1]
        %   would get the trace 1 second before to 1 second after the event)
        
            idx = trial.probe_t >= trigger_t + limits(1) & ...
                    trial.probe_t <= trigger_t + limits(2);
        end
    end
end
