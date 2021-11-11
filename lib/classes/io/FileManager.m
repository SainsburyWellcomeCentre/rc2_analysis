classdef FileManager < handle
    
    properties (SetAccess = private)
        
        ctl
        path_config
    end
    
    
    methods
        
        function obj = FileManager(ctl)
        %%class handling the names of files on the setup
            obj.ctl = ctl;
            obj.path_config = ctl.path_config;
        end
        
        
        
        function animal_list = list_animals(obj)
            
            contents = dir(obj.path_config.raw_probe_dir);
            idx = ~cellfun(@isempty, regexp({contents(:).name}, '^CA'));
            animal_list = {contents(idx).name}';
        end
        
        
        
        function [fname, exists] = experiment_list(obj)
            
            fname = fullfile(obj.path_config.summary_data_dir, 'experiment_list.csv');
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = rc2_dname(obj, session_id)
            
            animal_id = obj.animal_id_from_session_id(session_id);
            
            dname = fullfile(obj.path_config.raw_rc2_dir, animal_id, animal_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = rc2_bin(obj, session_id)
            
            fname = fullfile(obj.rc2_dname(session_id), sprintf('%s.bin', session_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = camera_csv_dir_fast(obj, session_id)
            
            camera_session_id = strrep(session_id, '_001', '');  % TODO: maybe session_id should not include the _001
            
            % first check fast dir
            dname = fullfile(obj.path_config.processed_camera_fast_dir, camera_session_id);
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = camera_csv_dir_slow(obj, session_id)
            
            camera_session_id = strrep(session_id, '_001', '');
            % then check slow dir
            dname = fullfile(obj.path_config.processed_camera_slow_dir, camera_session_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = camera_csv_fast(obj, session_id, camera_id)
            
            dname = obj.camera_csv_dir_fast(session_id);
            fname = fullfile(dname, sprintf('%s.csv', camera_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = camera_csv_slow(obj, session_id, camera_id)
            
            dname = obj.camera_csv_dir_slow(session_id);
            fname = fullfile(dname, sprintf('%s.csv', camera_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = raw_camera_dir(obj, session_id)
            
            camera_session_id = strrep(session_id, '_001', '');  % TODO: maybe session_id should not include the _001
            
            dname = fullfile(obj.path_config.raw_camera_dir, camera_session_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = raw_camera_avi(obj, session_id, camera_id)
            
            dname = obj.raw_camera_dir(session_id);
            fname = fullfile(dname, sprintf('camera%i.avi', camera_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_bin_raw(obj, probe_id, probe_type)
        %%from probe_id get ap.bin file name in raw data
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.ap.bin', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.bin', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_bin_raw(obj, probe_id, probe_type)
        %%from probe_id get lf.bin file name in raw data
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.lf.bin', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.bin', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_meta_raw(obj, probe_id, probe_type)
        %%from probe_id get ap.meta file name in raw data
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.ap.meta', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.meta', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_meta_raw(obj, probe_id, probe_type)
        %%from probe_id get lf.meta file name in raw data
            dname = obj.glx_bin_dir_raw(probe_id, probe_type);
            if strcmp(probe_type, '3A')
                fname = fullfile(dname, sprintf('%s_g0_t0.imec.lf.meta', probe_id));
            else
                fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.meta', probe_id));
            end
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_raw(obj, probe_id, probe_type)
            %%from probe_id get the directory in which the raw .ap.bin (and
            %%.lf.bin) files should be located
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            
            if strcmp(probe_type, '3A')
                dname = fullfile(obj.path_config.raw_probe_dir, animal_id);
            else
                dname = fullfile(obj.path_config.raw_probe_dir, ...
                                 animal_id, ...
                                 sprintf('%s_g0', probe_id), ...
                                 sprintf('%s_g0_imec0', probe_id));
            end
            
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = glx_ap_bin_processed_fast(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_bin_processed_fast(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_meta_processed_fast(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_meta_processed_fast(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_fast(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_processed_fast(obj, probe_id)
            
             % first try fast dir
            top_dir = obj.path_config.processed_probe_fast_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
            
        
        function [fname, exists] = glx_ap_bin_processed_slow(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_bin_processed_slow(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.bin', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_ap_meta_processed_slow(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.ap.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = glx_lf_meta_processed_slow(obj, probe_id)
            
            dname = obj.glx_bin_dir_processed_slow(probe_id);
            fname = fullfile(dname, sprintf('%s_g0_t0.imec0.lf.meta', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_processed_slow(obj, probe_id)
            
             % first try fast dir
            top_dir = obj.path_config.processed_probe_slow_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = glx_bin_dir_processed(obj, probe_id)
            
             % first try fast dir
            top_dir = obj.path_config.processed_probe_fast_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
            
            if exists; return; end
            
            % then try slow dir
            top_dir = obj.path_config.processed_probe_slow_dir;
            dname = obj.generate_glx_bin_dir_processed(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = formatted_data(obj, probe_id)
            
            fname = fullfile(obj.path_config.formatted_data_dir, sprintf('%s.mat', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = processed_output_dir_fast(obj, probe_id)
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_fast_dir, animal_id, 'output');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = processed_output_dir_slow(obj, probe_id)
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_slow_dir, animal_id, 'output');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = json_dir_fast(obj, probe_id)
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_fast_dir, animal_id, 'json_files');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = json_dir_slow(obj, probe_id)
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            dname = fullfile(obj.path_config.processed_probe_slow_dir, animal_id, 'json_files');
            exists = isfolder(dname);
        end
        
        
        
        function [dname, exists] = imec0_ks2(obj, probe_id)
            
            % first try fast dir
            top_dir = obj.path_config.processed_probe_fast_dir;
            dname = obj.generate_imec0_ks2(top_dir, probe_id);
            exists = isfolder(dname);
            
            if exists; return; end
            
            % then try slow dir
            top_dir = obj.path_config.processed_probe_slow_dir;
            dname = obj.generate_imec0_ks2(top_dir, probe_id);
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = trigger_mat(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'trigger.mat');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = original_trigger_mat(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'trigger_ori.mat');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = trigger_points_removed(obj, probe_id)
        
            fname = fullfile(obj.imec0_ks2(probe_id), 'trigger_points_removed.mat');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = selected_clusters_txt(obj, probe_id)
            
            dname = obj.imec0_ks2(probe_id);
            fname = fullfile(dname, 'selected_clusters.txt');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = svm_table(obj, probe_id)
            
            fname = fullfile(obj.path_config.summary_data_dir, 'stationary_vs_motion_fr', sprintf('%s.csv', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = offsets_table(obj, probe_id)
            
            fname = fullfile(obj.path_config.summary_data_dir, 'trial_matched_offsets', sprintf('%s.csv', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = tuning_curves(obj, probe_id)
            
            fname = fullfile(obj.path_config.summary_data_dir, 'tuning_curves', sprintf('%s.mat', probe_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = hf_power_parameters(obj, probe_id, shank_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('hf_power_%i.mat', shank_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = track_csv(obj, probe_id, shank_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('track_%i.csv', shank_id));
            exists = isfile(fname);
        end
           
        
        
        function [fname, exists] = track_offset(obj, probe_id, shank_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('offset_%i.txt', shank_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = ks2_npy(obj, probe_id, var)
            
            fname = fullfile(obj.imec0_ks2(probe_id), sprintf('%s.npy', var));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = cluster_groups(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'cluster_groups.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = ks2_label(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'cluster_KSLabel.tsv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = params(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'params.py');
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = imec0_ks2_csv_dir(obj, probe_id)
            
            dname = fullfile(obj.imec0_ks2(probe_id), 'csv');
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = metrics_csv(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'metrics.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = waveform_metrics_csv(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'waveform_metrics.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = waveform_metrics_fixed_csv(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'waveform_metrics_fix.csv');
            exists = isfile(fname);
            
            % if it doesn't exist attempt to find a later version
            if ~exists
                contents = dir(obj.imec0_ks2_csv_dir(probe_id));
                I = regexp({contents(:).name}, 'waveform_metrics_fix');
                idx = find(cellfun(@(x)(~isempty(x)), I));
                if ~isempty(idx)
                    fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), contents(idx(1)).name);
                    exists = isfile(fname);
                end
            end
        end
        
        
        
        function [fname, exists] = clusters_janelia_csv(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'clusters_janelia.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = clusters_janelia_xlsx(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2_csv_dir(probe_id), 'clusters_janelia.xlsx');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = hf_power_figure(obj, probe_id, shank_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'tracks', sprintf('hf_power_%i.pdf', shank_id));
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = tracks_dir(obj, probe_id)
            
            dname = fullfile(obj.imec0_ks2(probe_id), 'tracks');
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = driftmap(obj, probe_id)
            
            fname = fullfile(obj.imec0_ks2(probe_id), 'driftmap.pdf');
            exists = isfile(fname);
        end
        
        
        
        function probe_ids = probe_id_from_animal_id(obj, animal_id)
            
            animal_dir = fullfile(obj.path_config.raw_probe_dir, animal_id);
            
            % get the contents of the directory
            contents = dir(animal_dir);
            
            % find the directories beginning with the animal name
            idx = find(~cellfun(@isempty, regexp({contents(:).name}, animal_id)));
            
            if ~isempty(regexp(contents(idx(1)).name, '.bin', 'once'))
                probe_ids = regexprep({contents(idx).name}, '_g0.+', '');
                probe_ids = unique(probe_ids);
            else
                probe_ids = {contents(idx).name};
            end
        end
        
        
        
        function dname = generate_imec0_ks2(obj, top_dir, probe_id)
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            
            level_1 = fullfile(top_dir, animal_id, 'output');
            level_2 = fullfile(level_1, sprintf('catgt_%s_g0', probe_id));
            level_3 = fullfile(level_2, sprintf('%s_g0_imec0', probe_id));
            dname = fullfile(level_3, 'imec0_ks2');
        end
        
        
        
        function dname = generate_glx_bin_dir_processed(obj, top_dir, probe_id)
            
            animal_id = obj.animal_id_from_probe_id(probe_id);
            
            level_1 = fullfile(top_dir, animal_id);
            level_2 = fullfile(level_1, sprintf('%s_g0', probe_id));
            dname = fullfile(level_2, sprintf('%s_g0_imec0', probe_id));
%             fname = fullfile(level_3, sprintf('%s_g0_t0.imec0.ap.bin', probe_id));
        end
        
        

        function animal_id = animal_id_from_probe_id(obj, probe_id)
            
            animal_id = obj.ctl.animal_id_from_probe_id(probe_id);
            % for now assumes that rec1 is the non-animal string...
%             idx = regexp(probe_id, '_rec1');
            
%             animal_id = probe_id(1:idx(1)-1);
        end
    end
    
    
    
    methods (Static = true)
        
        
        
        
        
        function animal_id = animal_id_from_session_id(session_id)
            
            % for now assumes that rec is the non-animal string...
            idx = regexp(session_id, '_rec');
            
            animal_id = session_id(1:idx(1)-1);
        end
    end
end
