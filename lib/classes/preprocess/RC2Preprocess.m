classdef RC2Preprocess < RC2Format
    
    properties
    end
    
    methods
        
        function obj = RC2Preprocess()
        %%class for handling preprocessing of the data
        %   subclass of RC2Analysis
        
            obj = obj@RC2Format();
        end
        
        
        
        function preprocess_step_1(obj, probe_id)
        %%first step in preprocessing of the spike data
            
            obj.move_raw_to_local(probe_id);
            obj.janelia_ecephys_spike_sorting(probe_id);
            obj.create_check_clusters_csv(probe_id);
            obj.create_trigger_file(probe_id);
            obj.create_driftmap(probe_id);
            obj.process_camera_data(probe_id);
        end
        
        
        
        function janelia_ecephys_spike_sorting(obj, probe_id)
        %%runs the Janelia version of ecephys_spike_sorting
            je_helper = JaneliaEcephysHelper(obj, probe_id);
            je_helper.run_from_raw();
        end
        
        
        
        function create_check_clusters_csv(obj, probe_id)
            
            cc = RestrictClusters(obj, probe_id);
            new_tbl = cc.restrict_metrics_table();
            obj.save.clusters_janelia_csv(probe_id, new_tbl);
        end
           
        
        
        function create_trigger_file(obj, probe_id)
            
            rec = obj.load.spikeglx_ap_recording(probe_id);
            trigger = rec.data(rec.trigger_channel_idx, :);
            obj.save.trigger_mat(probe_id, trigger);
        end
        
        
        
        function ct = correct_trigger_file(obj, probe_id)
        %%opens window to correct the trigger file
            ct = CorrectTrigger(obj, probe_id);
        end
        
        
        
        function create_driftmap(obj, probe_id)
        %%create and save driftmap
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
            
            session_ids = obj.get_session_ids_list(probe_id);
            for ii = 1 : length(session_ids)
                ch = CameraProcessingHelper(obj, session_ids{ii});
                ch.run_from_raw();
            end
        end
        
        
        
        function hf_power = hf_power(obj, probe_id, shank_id)
        %%return a HighFrequencyPowerProfile object for probe_id and
        %%shank_id
            
            recording               = obj.load.spikeglx_ap_recording(probe_id);
            
            hf_power               = HighFrequencyPowerProfile(recording, probe_id, shank_id);
            
            probe_track            = obj.load_track(probe_id, shank_id);
            hf_power.probe_track   = probe_track;
            
            clusters_from_tip_um   = obj.mua_from_tip_um(probe_id, shank_id);
            hf_power.clusters_from_tip_um = clusters_from_tip_um;
        end
        
        
        
        function save_hf_power(obj, hf_power)
        %%save information from the hf_power object
        
            probe_id = hf_power.probe_id;
            shank_id = hf_power.shank_id;
            
            h_fig = hf_power.plot_summary();
            
            obj.save.track_offset(probe_id, shank_id, hf_power.delta_l5);
            obj.save.hf_power_figure(probe_id, shank_id, h_fig);
            obj.save.hf_power_parameters(probe_id, shank_id, hf_power.get_parameters());
        end
        
        
        
        function create_mock_track(obj, probe_id, shank_id, n_points, visp_n_points_from_tip)
            %%creates a mock probe 'track' with VISpX at top and 'Unknown' at
            %%bottom
            % open the file
            
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
            
            clusters_xlsx = obj.load.clusters_janelia_xlsx(probe_id);
            idx = ~(strcmp(clusters_xlsx.mateo, 'b') | strcmp(clusters_xlsx.lee, 'b'));
            selected_clusters = clusters_xlsx.cluster_id(idx);
            obj.save.selected_clusters_txt(probe_id, selected_clusters);
        end
        
        
        
        function clusters_from_tip_um = mua_from_tip_um(obj, probe_id, shank_id)
        %%for the high-power frequency plots, return the distances of
        %%'good' clusters from the probe tip
            
            clusters                = obj.format_clusters(probe_id);
            good_clusters           = strcmp({clusters(:).class}, 'good');
            on_shank                = [clusters(:).shank_id] == shank_id;
            cluster_mask            = good_clusters & on_shank;
            clusters_from_tip_um    = [clusters(cluster_mask).distance_from_probe_tip];
        end
        
        
        
%         function get_session_bounds(obj, probe_id)
%         %%from the pattern of triggers, get the boundaries of the sessions
%             
%             min_time_between_sessions = 2;  % seconds
%             
%             
%             
%         end
        
        
        
        function move_raw_to_local(obj, probe_id)
        %%moves the .ap.bin files from the remote server to the local fast
        %%drive
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
