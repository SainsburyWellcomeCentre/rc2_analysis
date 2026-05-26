classdef RC2Preprocess < RC2Format
% RC2Preprocess Class for preprocessing raw data
%
%   RC2Preprocess Properties:
%
%   RC2Preprocess Methods:
%       preprocess_step_1               - first step of preprocessing
%       run_from_step                   - run preprocess_step_1 starting from a chosen step
%       janelia_ecephys_spike_sorting   - run ecephys_spike_sorting
%       run_ecephys_from_step           - run ecephys_spike_sorting from a chosen module
%       create_check_clusters_csv       - create the .csv with the clusters to manually check
%       create_check_mua_clusters_csv   - create the .csv with the clusters to manually check for MUA clusters
%       create_trigger_file             - separate the trigger channel from the the probe .bin file
%       correct_trigger_file            - manually correct the trigger file if necessary
%       create_driftmap                 - create and save a driftmap
%       process_camera_data             - run processing of the camera data
%       hf_power                        - create HighFrequencyPowerProfile object
%       save_hf_power                   - save an analyzed HighFrequencyPowerProfile object
%       create_mock_track               - create a 'track.csv' with a similar format to the real eventual track file
%       create_selected_clusters_txt    - create text file with the clusters selected
%       create_selected_mua_clusters_txt - create text file with the MUA clusters selected
%       mua_from_tip_um                 - for MUA get disatnce from probe tip
%       move_raw_to_local               - move the raw probe data to a local location
%       patch_meta_NP2013               - patch NP2013 meta file for Janelia pipeline compatibility
%       cluster_info                    - create CheckClusterQuality object

    properties

        leave_window_open_on_error = false
    end

    methods

        function obj = RC2Preprocess()
        %%RC2Preprocess
        %
        %   RC2Preprocess() creates the object and prints a short usage
        %   summary describing the two ways to run stage 1 of the
        %   preprocessing.

            obj = obj@RC2Format();

            fprintf('RC2Preprocess ready. Stage 1 can be run in two ways:\n');
            fprintf('  1. preprocess_step_1(probe_id)            - run all 7 steps from the start (complete pipeline)\n');
            fprintf('  2. run_from_step(probe_id, start_step)    - run from a chosen step onwards (useful to debug)\n');
            fprintf('\n');
            fprintf('Valid start_step names for run_from_step (in execution order):\n');
            fprintf('  ''move_raw_to_local''\n');
            fprintf('  ''patch_meta_NP2013''\n');
            fprintf('  ''janelia_ecephys_spike_sorting''\n');
            fprintf('  ''create_check_clusters_csv''\n');
            fprintf('  ''create_trigger_file''\n');
            fprintf('  ''create_driftmap''\n');
            fprintf('  ''process_camera_data''\n');
            fprintf('\n');
            fprintf('To start the janelia step from a chosen sorting module AND still\n');
            fprintf('run every later step to the end, name the module directly as the\n');
            fprintf('start_step of run_from_step, e.g.:\n');
            fprintf('  run_from_step(probe_id, ''noise_templates'')\n');
            fprintf('Valid janelia sorting modules (in execution order):\n');
            fprintf('  ''kilosort_helper''\n');
            fprintf('  ''kilosort_postprocessing''\n');
            fprintf('  ''noise_templates''\n');
            fprintf('  ''mean_waveforms''\n');
            fprintf('  ''quality_metrics''\n');
            fprintf('\n');
            fprintf('To run *only* the ecephys_spike_sorting step (without chaining to the\n');
            fprintf('later steps), use:\n');
            fprintf('  3. run_ecephys_from_step(probe_id, ...)   - see ''help RC2Preprocess.run_ecephys_from_step''\n');
        end
        
        
        
        function preprocess_step_1(obj, probe_id)
        %%preprocess_step_1 First step of preprocessing
        %
        %   preprocess_step_1(PROBE_ID) runs stage 1 of the preprocessing
        %   for probe recording PROBE_ID. This includes:
        %       - moving raw probe data to a local location
        %       - patch the NP2013 meta file
        %       - run the ecephys_spike_sorting pipeline
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
        %   run_from_step(PROBE_ID, START_STEP) runs stage 1 of the
        %   preprocessing for probe recording PROBE_ID, starting from the
        %   step named START_STEP and running every subsequent step in
        %   order, *always to the end* of the pipeline. This matters because
        %   the whole analysis is interlinked: each step depends on the
        %   output of the previous ones, so resuming part-way through must
        %   never stop early.
        %
        %   START_STEP must be the name (char or string) of one of the 7
        %   steps below, listed in execution order:
        %       'move_raw_to_local'
        %       'patch_meta_NP2013'
        %       'janelia_ecephys_spike_sorting'
        %       'create_check_clusters_csv'
        %       'create_trigger_file'
        %       'create_driftmap'
        %       'process_camera_data'
        %
        %   START_STEP may ALSO be the name of one of the janelia sorting
        %   sub-modules, in which case the janelia step is resumed from that
        %   module and, crucially, the pipeline then carries on through
        %   create_check_clusters_csv and every later step (this is what
        %   plain run_ecephys_from_step does NOT do). Valid modules, in
        %   execution order:
        %       'kilosort_helper'
        %       'kilosort_postprocessing'
        %       'noise_templates'
        %       'mean_waveforms'
        %       'quality_metrics'
        %
        %   run_from_step(..., NAME, VALUE, ...) forwards the optional
        %   ecephys controls ('run_catgt', 'start_module', 'run_tprime') to
        %   the janelia step (see run_ecephys_from_step). These only have an
        %   effect when the janelia step is within the range being run.
        %
        %   When the janelia step is resumed from a module other than
        %   'kilosort_helper', CatGT has normally already been produced, so
        %   'run_catgt' defaults to false unless set explicitly.
        %
        %   Examples:
        %       ctl.run_from_step(probe_id, 'create_trigger_file')
        %   runs create_trigger_file, create_driftmap and
        %   process_camera_data, skipping the four earlier steps.
        %
        %       ctl.run_from_step(probe_id, 'noise_templates')
        %   resumes the janelia step from noise_templates (CatGT skipped by
        %   default), then runs create_check_clusters_csv, create_trigger_file,
        %   create_driftmap and process_camera_data.

            steps = {'move_raw_to_local', ...
                     'patch_meta_NP2013', ...
                     'janelia_ecephys_spike_sorting', ...
                     'create_check_clusters_csv', ...
                     'create_trigger_file', ...
                     'create_driftmap', ...
                     'process_camera_data'};

            % the janelia sorting sub-modules, in execution order; a module
            % name may be given directly as START_STEP (see below).
            ecephys_modules = JaneliaEcephysHelper.module_order;

            start_step   = char(start_step);
            ecephys_args = varargin;

            % shorthand: when START_STEP names a janelia sub-module, treat it
            % as "start the janelia step from this module". We translate it
            % into the janelia step plus a 'start_module' control so the loop
            % below still chains through every later step.
            if ismember(start_step, ecephys_modules)
                if any(strcmpi(ecephys_args(1:2:end), 'start_module'))
                    error('RC2Preprocess:run_from_step:duplicateStartModule', ...
                          ['START_STEP "%s" is a janelia module, so do not ' ...
                           'also pass a ''start_module'' name-value pair.'], ...
                          start_step);
                end
                ecephys_args = [{'start_module', start_step}, ecephys_args];
                start_step   = 'janelia_ecephys_spike_sorting';
            end

            % exact, case-sensitive match against the valid step names;
            % ismember returns the index of the match in START_IDX so no
            % approximate / partial matching can ever occur.
            [is_known_step, start_idx] = ismember(start_step, steps);

            if ~is_known_step
                error('RC2Preprocess:run_from_step:unknownStep', ...
                      ['Unknown step "%s". START_STEP must be the exact ' ...
                       'name of one of the pipeline steps:\n  %s\n' ...
                       'or one of the janelia sorting modules:\n  %s'], ...
                      start_step, strjoin(steps, '\n  '), ...
                      strjoin(ecephys_modules, '\n  '));
            end

            janelia_idx = find(strcmp('janelia_ecephys_spike_sorting', steps), 1);

            % the ecephys controls only apply if the janelia step is part of
            % the range about to run; warn and drop them otherwise.
            if ~isempty(ecephys_args) && start_idx > janelia_idx
                warning('RC2Preprocess:run_from_step:ignoredEcephysArgs', ...
                        ['Ecephys controls were supplied but START_STEP ' ...
                         '"%s" is after the janelia step, so they are ' ...
                         'ignored.'], start_step);
                ecephys_args = {};
            end

            % when resuming the janelia step from a module other than the
            % first one, CatGT has normally already run; default run_catgt to
            % false unless the caller set it explicitly.
            if ~isempty(ecephys_args) && mod(numel(ecephys_args), 2) == 0
                keys     = ecephys_args(1:2:end);
                sm_pos   = find(strcmpi(keys, 'start_module'), 1);
                has_catgt = any(strcmpi(keys, 'run_catgt'));
                if ~isempty(sm_pos) && ~has_catgt && ...
                        ~strcmp(ecephys_args{2*sm_pos}, 'kilosort_helper')
                    ecephys_args = [ecephys_args, {'run_catgt', false}];
                end
            end

            for ii = start_idx : length(steps)
                fprintf('Running step %i/%i: %s\n', ii, length(steps), steps{ii});
                if strcmp(steps{ii}, 'janelia_ecephys_spike_sorting') && ~isempty(ecephys_args)
                    % resume janelia from the chosen module / with the chosen
                    % switches, then let the loop continue to the later steps
                    obj.run_ecephys_from_step(probe_id, ecephys_args{:});
                else
                    obj.(steps{ii})(probe_id);
                end
            end
        end



        function janelia_ecephys_spike_sorting(obj, probe_id)
        %%janelia_ecephys_spike_sorting Run ecephys_spike_sorting
        %
        %   janelia_ecephys_spike_sorting(PROBE_ID) runs the full
        %   ecephys_spike_sorting pipeline for probe recording PROBE_ID
        %   (CatGT followed by all sorting modules, no TPrime).
        %
        %   To run the pipeline from a chosen point - skipping CatGT and/or
        %   the early sorting modules, or enabling TPrime - use
        %   run_ecephys_from_step.

            obj.run_ecephys_from_step(probe_id);
        end



        function run_ecephys_from_step(obj, probe_id, varargin)
        %%run_ecephys_from_step Run ecephys_spike_sorting from a chosen point
        %
        %   run_ecephys_from_step(PROBE_ID) runs the full ecephys_spike_sorting
        %   pipeline for probe recording PROBE_ID and is identical to
        %   janelia_ecephys_spike_sorting(PROBE_ID).
        %
        %   run_ecephys_from_step(PROBE_ID, NAME, VALUE, ...) runs the
        %   pipeline with the optional name-value controls below, so it can
        %   be re-run from any point without editing the python script:
        %
        %       'run_catgt'    - logical, whether to run the CatGT step
        %                        (default true). Set false to sort data that
        %                        has already been CatGT-processed.
        %       'start_module' - char, the sorting module to start from. That
        %                        module and every module after it are run.
        %                        One of, in execution order:
        %                            'kilosort_helper'         (default)
        %                            'kilosort_postprocessing'
        %                            'noise_templates'
        %                            'mean_waveforms'
        %                            'quality_metrics'
        %       'run_tprime'   - logical, whether to run the TPrime step at
        %                        the end of the pipeline (default false).
        %
        %   Example:
        %       ctl.run_ecephys_from_step(probe_id, 'run_catgt', false, ...
        %                                 'start_module', 'noise_templates')
        %   skips CatGT and the first two sorting modules, running only
        %   noise_templates, mean_waveforms and quality_metrics.
        %
        %   Note: starting from a later step assumes the output of the
        %   earlier steps already exists on disk.
        %
        %   This method runs the ecephys_spike_sorting step ONLY; it does not
        %   continue to create_check_clusters_csv or any later stage-1 step.
        %   To resume from a sorting module AND carry on to the end of the
        %   pipeline, use run_from_step (e.g.
        %   run_from_step(PROBE_ID, 'noise_templates')).

            parser = inputParser();
            parser.addParameter('run_catgt', true, ...
                                @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
            parser.addParameter('run_tprime', false, ...
                                @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
            parser.addParameter('start_module', 'kilosort_helper', ...
                                @(x) ischar(x) || isstring(x));
            parser.parse(varargin{:});

            je_helper = JaneliaEcephysHelper(obj, probe_id);
            je_helper.leave_window_open_on_error = obj.leave_window_open_on_error;
            je_helper.run_catgt  = logical(parser.Results.run_catgt);
            je_helper.run_tprime = logical(parser.Results.run_tprime);
            % the helper's set.start_module validates this against the known
            % module names and errors on an unknown one
            je_helper.start_module = char(parser.Results.start_module);

            fprintf(['Running ecephys_spike_sorting (run_catgt=%d, ' ...
                     'start_module=''%s'', run_tprime=%d)\n'], ...
                    je_helper.run_catgt, je_helper.start_module, je_helper.run_tprime);

            je_helper.run_from_raw();
        end
        
        
        
        function create_check_clusters_csv(obj, probe_id)
        %%create_check_clusters_csv Create the .csv with the clusters to manually check
        %
        %   create_check_clusters_csv(PROBE_ID) creates the .csv with the 
        %   clusters to manually check for probe recording PROBE_ID.
        
            cc = RestrictClusters(obj, probe_id);
            new_tbl = cc.restrict_metrics_table();
            obj.save.clusters_janelia_csv(probe_id, new_tbl);
        end
           
        
        
        function create_check_mua_clusters_csv(obj, probe_id)
        %%create_check_mua_clusters_csv Create the .csv with the clusters
        %%to manually check for MUA clusters
        %
        %   create_check_mua_clusters_csv(PROBE_ID) creates the .csv with
        %   the clusters to manually check for probe recording PROBE_ID. 
        
            cc = RestrictClusters(obj, probe_id);
            new_tbl = cc.restrict_mua_metrics_table();
            obj.save.mua_clusters_janelia_csv(probe_id, new_tbl);
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
        
            ks2_dir = obj.file.imec0_ks2(probe_id);
            [spikeTimes, spikeAmps, spikeDepths] = ksDriftmap(ks2_dir);
            
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
        %%create_selected_clusters_txt Create text file with the clusters selected
        %
        %   create_selected_clusters_txt(PROBE_ID) takes the .xlsx saved
        %   after the manual Phy step and extracts the selected clusters to
        %   save in a text file in the Kilosort directory.
            
            clusters_xlsx = obj.load.clusters_janelia_xlsx(probe_id);
            idx = ~(strcmp(clusters_xlsx.mateo, 'b') & strcmp(clusters_xlsx.lee, 'b'));
            selected_clusters = clusters_xlsx.cluster_id(idx);
            obj.save.selected_clusters_txt(probe_id, selected_clusters);
        end
        
        
        
        function create_selected_mua_clusters_txt(obj, probe_id)
        %%create_selected_mua_clusters_txt Create text file with a list of
        %%the selected MUA clusters 
        %
        %   create_selected_mua_clusters_txt(PROBE_ID) takes the .xlsx
        %   saved after the manual Phy step and extracts the selected MUA
        %   clusters to save in a text file in the Kilosort directory.
            
            clusters_xlsx = obj.load.mua_clusters_janelia_xlsx(probe_id);
            idx = ~((clusters_xlsx.mateo == 1) | (clusters_xlsx.lee == 1));
            selected_clusters = clusters_xlsx.cluster_id(idx);
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
        %%patch_meta_NP2013 Patch NP2013 meta file for Janelia pipeline compatibility
        %
        %   patch_meta_NP2013(PROBE_ID) rewrites the locally-copied AP meta
        %   file for probe recording PROBE_ID so that the (newer SpikeGLX)
        %   NP2013 / probe-type 2013 fields are remapped to the NP2010 /
        %   probe-type 24 form, and the `snsGeomMap` key is renamed to
        %   `snsShankMap`. This matches the format the Janelia
        %   `ecephys_spike_sorting` pipeline expects (`findDisabled()` in
        %   `SGLXMetaToCoords.py` looks up `snsShankMap` directly).
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
