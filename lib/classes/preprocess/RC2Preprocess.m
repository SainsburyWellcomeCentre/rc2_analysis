classdef RC2Preprocess < RC2Format
% RC2Preprocess Class for preprocessing raw data
%
%   Stage 1 can be run in two ways:
%     1. preprocess_step_1(probe_id)            - run all 7 steps from the start (complete pipeline)
%     2. run_from_step(probe_id, start_step)    - run from a chosen step onwards, see 'help RC2Preprocess.run_from_step'
%
%   To run only the sorting step, with finer control over where it resumes
%   from, see 'help RC2Preprocess.run_sorting_from_step'.
%
%   See also the main README ("Preprocessing and Analysis Workflow" and
%   "Tips for debugging").
%
%   RC2Preprocess Properties:
%
%   RC2Preprocess Methods:
%       preprocess_step_1               - first step of preprocessing
%       run_from_step                   - run preprocess_step_1 starting from a chosen step
%       si_sorting                      - run the SpikeInterface sorting pipeline
%       run_sorting_from_step           - run the SpikeInterface pipeline with optional controls
%       create_check_clusters_csv       - write the automated (Bombcell+metric) curation .csv and selected_clusters.txt
%       create_check_mua_clusters_csv   - write the automated MUA curation .csv and selected_mua_clusters.txt
%       create_trigger_file             - separate the trigger channel from the the probe .bin file
%       correct_trigger_file            - manually correct the trigger file if necessary
%       create_driftmap                 - create and save a driftmap
%       process_camera_data             - run processing of the camera data
%       hf_power                        - create HighFrequencyPowerProfile object
%       save_hf_power                   - save an analyzed HighFrequencyPowerProfile object
%       create_mock_track               - create a 'track.csv' with a similar format to the real eventual track file
%       create_selected_clusters_txt    - write selected_clusters.txt from the `keep` column of the curation .csv
%       create_selected_mua_clusters_txt - write selected_mua_clusters.txt from the `keep` column of the MUA curation .csv
%       mua_from_tip_um                 - for MUA get disatnce from probe tip
%       move_raw_to_local               - move the raw probe data to a local location
%       patch_meta_NP2013               - remap a newer-SpikeGLX NP2013 meta file to the older NP24 field layout
%       cluster_info                    - create CheckClusterQuality object

    properties

        leave_window_open_on_error = false
    end

    methods

        function obj = RC2Preprocess()
        %%RC2Preprocess
        %
        %   RC2Preprocess() creates the object.

            obj = obj@RC2Format();
            fprintf('See ''help RC2Preprocess'' for how to run stage 1.\n');
        end
        
        
        
        function preprocess_step_1(obj, probe_id)
        %%preprocess_step_1 First step of preprocessing
        %
        %   preprocess_step_1(PROBE_ID) runs stage 1 of the preprocessing
        %   for probe recording PROBE_ID. This includes:
        %       - moving raw probe data to a local location
        %       - patch the NP2013 meta file
        %       - run the SpikeInterface sorting pipeline
        %       - create a .csv with clusters to check
        %       - create a .mat with the trigger channel
        %       - create a driftmap
        %       - process the motion energy from the camera data
        %
        %   To start from a step other than the first one, use
        %   run_from_step.

            obj.run_from_step(probe_id, 'move_raw_to_local');
        end



        function run_from_step(obj, probe_id, start_step, varargin)
        %%run_from_step Run preprocess_step_1 starting from a chosen step
        %
        %   run_from_step(PROBE_ID, START_STEP) runs stage 1 for probe
        %   recording PROBE_ID from the step named START_STEP to the end of
        %   the pipeline (never stops early -- later steps depend on earlier
        %   ones). START_STEP is one of the 7 steps below, or one of the
        %   sorting sub-steps ('kilosort4', 'postprocess', 'bombcell' -- see
        %   run_sorting_from_step), in which case si_sorting resumes from
        %   there and the pipeline still continues through every later step:
        %       'move_raw_to_local'
        %       'patch_meta_NP2013'
        %       'si_sorting'
        %       'create_check_clusters_csv'
        %       'create_trigger_file'
        %       'create_driftmap'
        %       'process_camera_data'
        %
        %   run_from_step(..., NAME, VALUE, ...) forwards optional pipeline
        %   controls ('run_catgt', 'run_tprime', see run_sorting_from_step) to
        %   the sorting step, if it is within the range being run.

            steps = {'move_raw_to_local', ...
                     'patch_meta_NP2013', ...
                     'si_sorting', ...
                     'create_check_clusters_csv', ...
                     'create_trigger_file', ...
                     'create_driftmap', ...
                     'process_camera_data'};

            % the sorting sub-steps; a sub-step name may be given directly
            % as START_STEP (see above).
            sorting_sub_steps = {'kilosort4', 'postprocess', 'bombcell'};

            start_step   = char(start_step);
            sorting_args = varargin;

            % shorthand: when START_STEP names a sorting sub-step, treat it
            % as "start the sorting step from this sub-step". Translated
            % into the si_sorting step plus a 'start_step' control so the
            % loop below still chains through every later step.
            if ismember(start_step, sorting_sub_steps)
                if any(strcmpi(sorting_args(1:2:end), 'start_step'))
                    error('RC2Preprocess:run_from_step:duplicateStartStep', ...
                          ['START_STEP "%s" is a sorting sub-step, so do not ' ...
                           'also pass a ''start_step'' name-value pair.'], ...
                          start_step);
                end
                sorting_args = [{'start_step', start_step}, sorting_args];
                start_step   = 'si_sorting';
            end

            % exact, case-sensitive match against the valid step names
            [is_known_step, start_idx] = ismember(start_step, steps);

            if ~is_known_step
                error('RC2Preprocess:run_from_step:unknownStep', ...
                      ['Unknown step "%s". START_STEP must be the exact ' ...
                       'name of one of the pipeline steps:\n  %s\n' ...
                       'or one of the sorting sub-steps:\n  %s'], ...
                      start_step, strjoin(steps, '\n  '), ...
                      strjoin(sorting_sub_steps, '\n  '));
            end

            si_idx = find(strcmp('si_sorting', steps), 1);

            % pipeline controls (run_catgt, run_tprime) only apply if the
            % sorting step is within the range being run; warn and drop otherwise.
            if ~isempty(sorting_args) && start_idx > si_idx
                warning('RC2Preprocess:run_from_step:ignoredPipelineArgs', ...
                        ['Pipeline controls were supplied but START_STEP ' ...
                         '"%s" is after the sorting step, so they are ' ...
                         'ignored.'], start_step);
                sorting_args = {};
            end

            for ii = start_idx : length(steps)
                fprintf('Running step %i/%i: %s\n', ii, length(steps), steps{ii});
                if strcmp(steps{ii}, 'si_sorting') && ~isempty(sorting_args)
                    % run sorting with the pipeline controls supplied by the caller
                    obj.run_sorting_from_step(probe_id, sorting_args{:});
                else
                    obj.(steps{ii})(probe_id);
                end
            end
        end



        function si_sorting(obj, probe_id)
        %%si_sorting Run the SpikeInterface sorting pipeline
        %
        %   si_sorting(PROBE_ID) runs the full SpikeInterface pipeline for
        %   probe recording PROBE_ID (CatGT + KS4 + SI post-processing,
        %   no TPrime).
        %
        %   To skip CatGT or enable TPrime, use run_sorting_from_step.

            obj.run_sorting_from_step(probe_id);
        end



        function run_sorting_from_step(obj, probe_id, varargin)
        %%run_sorting_from_step Run the SpikeInterface sorting pipeline with optional controls
        %
        %   run_sorting_from_step(PROBE_ID) runs the full SpikeInterface
        %   pipeline for probe recording PROBE_ID (identical to
        %   si_sorting(PROBE_ID)). Runs the sorting step ONLY -- it does not
        %   continue to create_check_clusters_csv or any later stage-1 step
        %   (use run_from_step for that).
        %
        %   run_sorting_from_step(PROBE_ID, NAME, VALUE, ...) with the
        %   optional controls below:
        %
        %       'run_catgt'  - logical, whether to run the CatGT step
        %                      (default true, or false if 'start_step' is
        %                      given and is not 'catgt'). Set false to
        %                      re-sort data that has already been
        %                      CatGT-processed.
        %       'run_tprime' - logical, whether to run the TPrime step at
        %                      the end of the pipeline (default false).
        %       'start_step' - 'catgt' (default), 'kilosort4', 'postprocess'
        %                      or 'bombcell': where in the sorting pipeline
        %                      to resume from. Requires the earlier steps'
        %                      output to already exist on disk for this
        %                      probe:
        %                        'catgt'       - full run from CatGT onwards
        %                        'kilosort4'   - skip CatGT, read the existing
        %                                        CatGT output, destripe, run
        %                                        Kilosort4 and everything after
        %                        'postprocess' - skip CatGT and Kilosort4,
        %                                        reload the existing sort, then
        %                                        recompute SortingAnalyzer,
        %                                        Bombcell, Phy export and CSV
        %                                        export
        %                        'bombcell'    - skip CatGT, Kilosort4 AND
        %                                        SortingAnalyzer -- reload the
        %                                        existing SortingAnalyzer from
        %                                        disk, then only rerun Bombcell,
        %                                        Phy export and CSV export
        %
        %   Examples:
        %       ctl.run_sorting_from_step(probe_id, 'start_step', 'postprocess')
        %   reruns SortingAnalyzer onwards only, reloading the existing
        %   Kilosort4 output -- useful after fixing a Bombcell/export bug
        %   without redoing Kilosort4.
        %
        %       ctl.run_sorting_from_step(probe_id, 'start_step', 'bombcell')
        %   reruns only Bombcell onwards, reloading the existing SortingAnalyzer
        %   -- useful after changing Bombcell thresholds only, since
        %   SortingAnalyzer's metrics do not depend on them.

            valid_start_steps = {'catgt', 'kilosort4', 'postprocess', 'bombcell'};

            parser = inputParser();
            parser.addParameter('run_catgt', [], ...
                                @(x) isempty(x) || (isscalar(x) && (islogical(x) || isnumeric(x))));
            parser.addParameter('run_tprime', false, ...
                                @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
            parser.addParameter('start_step', 'catgt', ...
                                @(x) ismember(x, valid_start_steps));
            parser.parse(varargin{:});

            % default run_catgt to false whenever resuming past CatGT, unless
            % the caller explicitly overrides it (matches the old janelia
            % run_ecephys_from_step behaviour for 'start_module').
            run_catgt = parser.Results.run_catgt;
            if isempty(run_catgt)
                run_catgt = strcmp(parser.Results.start_step, 'catgt');
            end

            si_helper = SortingHelper(obj, probe_id);
            si_helper.leave_window_open_on_error = obj.leave_window_open_on_error;
            si_helper.run_catgt  = logical(run_catgt);
            si_helper.run_tprime = logical(parser.Results.run_tprime);
            si_helper.start_step = parser.Results.start_step;

            fprintf('Running SpikeInterface pipeline (start_step=%s, run_catgt=%d, run_tprime=%d)\n', ...
                    si_helper.start_step, si_helper.run_catgt, si_helper.run_tprime);

            si_helper.run_from_raw();
        end
        
        
        
        function create_check_clusters_csv(obj, probe_id)
        %%create_check_clusters_csv Create the automated (Bombcell + metric) curation .csv
        %
        %   create_check_clusters_csv(PROBE_ID) writes an EDITABLE .csv listing
        %   every cluster with its Bombcell label and two columns pre-filled
        %   with the same automated decision (Bombcell 'good' AND passing the
        %   quality metrics): `keep_pipeline` and `keep`. It then
        %   auto-generates selected_clusters.txt from those defaults, so
        %   curation is fully automated and NO manual step is required.
        %
        %   Optional manual override: edit the `keep` column of
        %   clusters_to_check.csv (optionally after inspecting units in Phy) to
        %   force a cluster in (keep=1) or out (keep=0), then re-run
        %   create_selected_clusters_txt(PROBE_ID) before formatting. Leave
        %   `keep_pipeline` untouched -- it stays a fixed record of the
        %   automated decision, so any hand-edit can be compared back against
        %   it later (e.g. sum(tbl.keep ~= tbl.keep_pipeline) to count how many
        %   clusters were manually overridden). Only `keep` is ever read by
        %   create_selected_clusters_txt.

            cc = RestrictClusters(obj, probe_id);
            tbl = cc.curation_table();
            obj.save.clusters_to_check_csv(probe_id, tbl);
            obj.create_selected_clusters_txt(probe_id);
        end
           
        
        
        function create_check_mua_clusters_csv(obj, probe_id)
        %%create_check_mua_clusters_csv Create the automated MUA curation .csv
        %
        %   create_check_mua_clusters_csv(PROBE_ID) writes an EDITABLE .csv with
        %   a `keep` column pre-filled from the automated MUA decision (Bombcell
        %   'mua' AND passing the MUA metric thresholds), then auto-generates
        %   selected_mua_clusters.txt. Edit the `keep` column and re-run
        %   create_selected_mua_clusters_txt(PROBE_ID) to override.

            cc = RestrictClusters(obj, probe_id);
            tbl = cc.curation_mua_table();
            obj.save.mua_clusters_to_check_csv(probe_id, tbl);
            obj.create_selected_mua_clusters_txt(probe_id);
        end
        
        
        
        function create_trigger_file(obj, probe_id)
        %%create_trigger_file Separate the trigger channel from the the probe .bin file
        %
        %   create_trigger_file(PROBE_ID) separate the trigger channel from 
        %   the the probe .bin file and save in a .mat, for a probe
        %   recording PROBE_ID.
        
            rec = obj.load.spikeglx_ap_recording(probe_id);
            trigger = rec.data(rec.trigger_channel_idx, :);
            obj.save.trigger_mat(probe_id, trigger);
        end
        
        
        
        function ct = correct_trigger_file(obj, probe_id)
        %%correct_trigger_file Manually correct the trigger file if necessary
        %
        %   GUI_HANDLE = correct_trigger_file(PROBE_ID)
        %   opens a small GUI to correct the trigger trace.
        %   Returns the handle to the GUI opened in GUI_HANDLE.
        
            ct = CorrectTrigger(obj, probe_id);
        end
        
        
        
        function create_driftmap(obj, probe_id)
        %%create_driftmap Create and save a driftmap
        %
        %   create_driftmap(PROBE_ID) creates and save a driftmap for a
        %   probe recording PROBE_ID.

            % ksDriftmap (cortex-lab/spikes) loads pc_features.npy /
            % pc_feature_ind.npy unconditionally (params.loadPCs = true) and
            % uses them directly to compute spikeDepths (weighted centroid
            % of the first PC across channels) -- not an optional/unused
            % load. Those two files are deliberately NOT copied to the
            % imec0_ks4 root by the sorting pipeline (copy_ks4_outputs_to_parent,
            % ~2 GB, unused by rc2_analysis's own FileManager/Loader), so
            % point at sorter_output/ instead, where KS4 always writes them,
            % rather than duplicating the files just for this one caller.
            ks4_dir = fullfile(obj.file.imec0_ks4(probe_id), 'sorter_output');
            [spikeTimes, spikeAmps, spikeDepths] = ksDriftmap(ks4_dir);
            
            % render the driftmap off-screen: it is a very dense scatter
            % that is saved straight to PDF and closed, never inspected
            % live, and drawing it on screen can stall the graphics
            % subsystem ("graphics handshaking" timeout)
            h_fig = figure('Visible', 'off');
            plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
            set(h_fig, 'position', [75, 158, 1041, 778]);
            box off;
            title(probe_id, 'interpreter', 'none');
            obj.save.driftmap(probe_id, h_fig);
            close(h_fig);
        end
        
        
        
        function process_camera_data(obj, probe_id)
        %%process_camera_data Process the motion energy for camera data
        %
        %   process_camera_data(PROBE_ID) computes and save the motion
        %   energy for the camera data for the probe recording PROBE_ID. 
        
            session_ids = obj.get_session_ids_list(probe_id);
            for ii = 1 : length(session_ids)
                ch = CameraProcessingHelper(obj, session_ids{ii});
                ch.run_from_raw();
            end
        end
        
        
        
        function hf_power = hf_power(obj, probe_id, shank_id)
        %%hf_power Create HighFrequencyPowerProfile object
        %
        %   HF_POWER = hf_power(PROBE_ID, SHANK_ID)
        %   return HighFrequencyPowerProfile object for probe recording
        %   PROBE_ID and shank SHANK_ID. Used to examine peaks in the
        %   high-frequency power profile.        
            
            recording               = obj.load.spikeglx_ap_recording(probe_id);
            
            hf_power               = HighFrequencyPowerProfile(recording, probe_id, shank_id);
            
            probe_track            = obj.load_track(probe_id, shank_id);
            hf_power.probe_track   = probe_track;
            
            clusters_from_tip_um   = obj.mua_from_tip_um(probe_id, shank_id);
            hf_power.clusters_from_tip_um = clusters_from_tip_um;
        end
        
        
        
        function save_hf_power(obj, hf_power)
        %%save_hf_power  Save the HighFrequencyPowerProfile object
        %
        %   save_hf_power(HF_POWER). After analyzing the high-frequency
        %   power with the HighFrequencyPowerProfile object, this can be
        %   used to save the details.
        
            probe_id = hf_power.probe_id;
            shank_id = hf_power.shank_id;
            
            h_fig = hf_power.plot_summary();
            
            obj.save.create_tracks_dir(probe_id);
            obj.save.track_offset(probe_id, shank_id, hf_power.delta_l5);
            obj.save.hf_power_figure(probe_id, shank_id, h_fig);
            obj.save.hf_power_parameters(probe_id, shank_id, hf_power.get_parameters());
        end
        
        
        
        function create_mock_track(obj, probe_id, shank_id, n_points, visp_n_points_from_tip)
        %%create_mock_track Create a 'track.csv' with a similar format to the real eventual track file
        %
        %   create_mock_track(PROBE_ID, SHANK_ID, N_POINTS, VISP_N_POINTS_FROM_TIP)
        %   creates a .csv file with a similar structure to the anatomical
        %   track.csv file for probe recording PROBE_ID and shank SHANK_ID.
        %    Creates N_POINTS rows in the .csv, with 'VISpX' in the bottom
        %    VISP_N_POINTS_FROM_TIP rows and 'Unknown' above that.
            
            [fname, exists] = obj.file.track_csv(probe_id, shank_id);
            
            if exists
                user = input('File exists, overwrite (Y)?', 's');
                if ~strcmp(user, 'Y')
                    return
                end
            end
            
            fid = fopen(fname, 'w');
            
            % header
            fwrite(fid, 'Position,Region ID,Region acronym,Region name');
            fwrite(fid, newline);
            
            % for each point
            for i = 1 : n_points
                
                if i <= n_points - visp_n_points_from_tip
                    str = sprintf('%i,%i,%s,"%s"', i-1, -2, 'VISpX', 'Primary visual area, layer X');
                    fwrite(fid, str);
                else
                    str = sprintf('%i,%i,%s,"%s"', i-1, -1, 'Unknown', 'Unknown');
                    fwrite(fid, str);
                end
                fwrite(fid, newline);
                
            end
            
            fclose(fid);
            
            % save a zero to the offset.txt file
            obj.save.track_offset(probe_id, shank_id, 0);
        end
        
        
        
        function create_selected_clusters_txt(obj, probe_id)
        %%create_selected_clusters_txt Write selected_clusters.txt from the `keep` column of the curation .csv
        %
        %   create_selected_clusters_txt(PROBE_ID) reads clusters_to_check.csv
        %   (the automated Bombcell + metric curation, possibly hand-edited) and
        %   writes the IDs whose `keep` column is 1 to selected_clusters.txt.
        %   Called automatically by create_check_clusters_csv; re-run it after
        %   editing the .csv to apply manual keep/discard overrides before
        %   formatting.

            tbl = obj.load.clusters_to_check_csv(probe_id);
            selected_clusters = tbl.cluster_id(tbl.keep == 1);
            obj.save.selected_clusters_txt(probe_id, selected_clusters);
        end
        
        
        
        function create_selected_mua_clusters_txt(obj, probe_id)
        %%create_selected_mua_clusters_txt Write selected_mua_clusters.txt from the `keep` column of the MUA curation .csv
        %
        %   create_selected_mua_clusters_txt(PROBE_ID) reads
        %   mua_clusters_to_check.csv (automated MUA curation, possibly
        %   hand-edited) and writes the IDs whose `keep` column is 1 to
        %   selected_mua_clusters.txt. Called automatically by
        %   create_check_mua_clusters_csv; re-run after editing the .csv.

            tbl = obj.load.mua_clusters_to_check_csv(probe_id);
            selected_clusters = tbl.cluster_id(tbl.keep == 1);
            obj.save.selected_mua_clusters_txt(probe_id, selected_clusters);
        end
        
        
        
        function clusters_from_tip_um = mua_from_tip_um(obj, probe_id, shank_id)
        %%mua_from_tip_um For MUA get disatnce from probe tip
        %
        %   FROM_TIP = mua_from_tip_um(PROBE_ID, SHANK_ID) gets the multiunit activity
        %   units for probe recording PROBE_ID, and which lie on shank
        %   SHANK_ID.
            
            clusters                = obj.format_clusters(probe_id);
            good_clusters           = strcmp({clusters(:).class}, 'good');
            on_shank                = [clusters(:).shank_id] == shank_id;
            cluster_mask            = good_clusters & on_shank;
            clusters_from_tip_um    = [clusters(cluster_mask).distance_from_probe_tip];
        end
        
        
        
        function move_raw_to_local(obj, probe_id)
        %%move_raw_to_local Move the raw probe data to a local location
        %
        %   move_raw_to_local(PROBE_ID) moves the .ap.bin files from the 
        %   remote server to the local fast drive.
        %
        %   TODO: add option to choose the remote location
        
            probe_type = obj.get_probe_type_from_experimentlist(probe_id);
            
            remote_ap_bin = obj.file.glx_ap_bin_raw(probe_id, probe_type);
            remote_ap_meta = obj.file.glx_ap_meta_raw(probe_id, probe_type);
            remote_lf_bin = obj.file.glx_lf_bin_raw(probe_id, probe_type);
            remote_lf_meta = obj.file.glx_lf_meta_raw(probe_id, probe_type);
            
            local_ap_bin = obj.file.glx_ap_bin_processed_fast(probe_id);
            local_ap_meta = obj.file.glx_ap_meta_processed_fast(probe_id);
            local_lf_bin = obj.file.glx_lf_bin_processed_fast(probe_id);
            local_lf_meta = obj.file.glx_lf_meta_processed_fast(probe_id);
            
            % copy the files from remote to local
            obj.xcopy(remote_ap_bin, local_ap_bin);
            obj.xcopy(remote_ap_meta, local_ap_meta);
            obj.xcopy(remote_lf_bin, local_lf_bin);
            obj.xcopy(remote_lf_meta, local_lf_meta);
        end



        function patch_meta_NP2013(obj, probe_id)
        %%patch_meta_NP2013 Patch a newer-SpikeGLX NP2013 meta file to the older NP24 layout
        %
        %   patch_meta_NP2013(PROBE_ID) rewrites the locally-copied AP meta
        %   file for probe recording PROBE_ID so that the (newer SpikeGLX)
        %   NP2013 / probe-type 2013 fields are remapped to the NP2010 /
        %   probe-type 24 form, and the `snsGeomMap` key is renamed to
        %   `snsShankMap`. This normalises the meta to the older field layout
        %   that legacy tooling reads (they look up `snsShankMap` directly).
        %   NOTE: the SpikeInterface reader uses `snsGeomMap` natively, so this
        %   step is likely unnecessary for the current pipeline -- review /
        %   remove if no downstream tool needs `snsShankMap`.
        %
        %   No-op when the meta file does not contain `imDatPrb_pn=NP2013`,
        %   so this is safe to call on already-compatible recordings.
        %
        %   Substitutions applied:
        %       imDatPrb_pn=NP2013         -> imDatPrb_pn=NP2010
        %       imDatPrb_type=2013         -> imDatPrb_type=24
        %       ~imroTbl=(2013,            -> ~imroTbl=(24,
        %       ~snsGeomMap=(NP2013,       -> ~snsShankMap=(

            local_ap_meta = obj.file.glx_ap_meta_processed_fast(probe_id);

            if ~isfile(local_ap_meta)
                warning('patch_meta_NP2013: meta file not found: %s', local_ap_meta);
                return
            end

            fid = fopen(local_ap_meta, 'r');
            text = fread(fid, '*char')';
            fclose(fid);

            if ~contains(text, 'imDatPrb_pn=NP2013')
                return
            end

            fprintf('Patching NP2013 meta to NP24 form: %s\n', local_ap_meta);

            text = strrep(text, 'imDatPrb_pn=NP2013',     'imDatPrb_pn=NP2010');
            text = strrep(text, 'imDatPrb_type=2013',     'imDatPrb_type=24');
            text = strrep(text, '~imroTbl=(2013,',        '~imroTbl=(24,');
            text = strrep(text, '~snsGeomMap=(NP2013,',   '~snsShankMap=(');

            fid = fopen(local_ap_meta, 'w');
            fwrite(fid, text);
            fclose(fid);
        end



        function cinfo = cluster_info(obj, probe_id)
        %%cluster_info Create CheckClusterQuality object
        %
        %   cluster_info(PROBE_ID) creates a CheckClusterQuality object for
        %   probe recording PROBE_ID, to view some quality metrics.
        
            cinfo = CheckClusterQuality(obj, probe_id);
        end
    end
    
    
    
    methods (Static = true)
        
        function xcopy(remote, local)
            fprintf('Moving %s to %s\n', remote, local);
            cmd = sprintf('echo f | xcopy /F /D "%s" "%s"', remote, local);
            system(cmd);
        end
    end
end
