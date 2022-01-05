classdef RC2Preprocess < RC2Format
% RC2Preprocess Class for preprocessing raw data
%
%   RC2Preprocess Properties:
%
%   RC2Preprocess Methods:
%       preprocess_step_1               - first step of preprocessing
%       janelia_ecephys_spike_sorting   - run ecephys_spike_sorting
%       create_check_clusters_csv       - create the .csv with the clusters to manually check
%       create_trigger_file             - separate the trigger channel from the the probe .bin file
%       correct_trigger_file            - manually correct the trigger file if necessary
%       create_driftmap                 - create and save a driftmap
%       process_camera_data             - run processing of the camera data
%       hf_power                        - create HighFrequencyPowerProfile object
%       save_hf_power                   - save an analyzed HighFrequencyPowerProfile object
%       create_mock_track               - create a 'track.csv' with a similar format to the real eventual track file
%       create_selected_clusters_txt    - create text file with the clusters selected
%       mua_from_tip_um                 - for MUA get disatnce from probe tip
%       move_raw_to_local               - move the raw probe data to a local location
%       cluster_info                    - create CheckClusterQuality object

    properties
    end
    
    methods
        
        function obj = RC2Preprocess()
        %%RC2Preprocess
        %
        %   RC2Preprocess()
        
            obj = obj@RC2Format();
        end
        
        
        
        function preprocess_step_1(obj, probe_id)
        %%preprocess_step_1 First step of preprocessing
        %
        %   preprocess_step_1(PROBE_ID) runs stage 1 of the preprocessing 
        %   for probe recording PROBE_ID. This includes:
        %       - moving raw probe data to a local location
        %       - run the ecephys_spike_sorting pipeline
        %       - create a .csv with clusters to check
        %       - create a .mat with the trigger channel
        %       - create a driftmap
        %       - process the motion energy from the camera data
            
            obj.move_raw_to_local(probe_id);
            obj.janelia_ecephys_spike_sorting(probe_id);
            obj.create_check_clusters_csv(probe_id);
            obj.create_trigger_file(probe_id);
            obj.create_driftmap(probe_id);
            obj.process_camera_data(probe_id);
        end
        
        
        
        function janelia_ecephys_spike_sorting(obj, probe_id)
        %%janelia_ecephys_spike_sorting Run ecephys_spike_sorting
        %
        %   janelia_ecephys_spike_sorting(PROBE_ID) runs the
        %   ecephys_spike_sorting pipeline for probe recording PROBE_ID
        
            je_helper = JaneliaEcephysHelper(obj, probe_id);
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
            
            figure()
            plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
            set(gcf, 'position', [75, 158, 1041, 778]);
            box off;
            title(probe_id, 'interpreter', 'none');
            obj.save.driftmap(probe_id, gcf);
            close(gcf);
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
            idx = ~(strcmp(clusters_xlsx.mateo, 'b') | strcmp(clusters_xlsx.lee, 'b'));
            selected_clusters = clusters_xlsx.cluster_id(idx);
            obj.save.selected_clusters_txt(probe_id, selected_clusters);
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
