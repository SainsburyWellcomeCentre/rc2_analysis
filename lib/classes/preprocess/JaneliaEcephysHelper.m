classdef JaneliaEcephysHelper < handle
    
    properties
        
        leave_window_open_on_error = false
    end
    
    properties (Hidden = true)
        
        chose_src_spikeglx_dir = false;
    end
    
    properties (SetAccess = private)
        
        ctl
        probe_id
    end
    
    
    
    methods
        
        function obj = JaneliaEcephysHelper(ctl, probe_id)
            obj.ctl = ctl;
            obj.probe_id = probe_id;
        end
        
        
        
        function run_from_raw(obj)
        %%run from raw data
            obj.create_output_dirs();
            obj.overwrite_ecephys_py_script();
            obj.run_ecephys();
            obj.move_csvs();
            obj.move_chanmap();
            obj.fix_waveforms();
        end
        
        
        
        function create_output_dirs(obj)
            
            output_dir = obj.ctl.file.processed_output_dir_fast(obj.probe_id);
            json_dir = obj.ctl.file.json_dir_fast(obj.probe_id);
            
            if isfolder(output_dir)
                fprintf('%s already exists\n', output_dir);
            else
                fprintf('Making %s\n', output_dir);
                mkdir(output_dir)
            end
            
            if isfolder(json_dir)
                fprintf('%s already exists\n', json_dir);
            else
                fprintf('Making %s\n', json_dir);
                mkdir(json_dir)
            end
        end
        
        
        
        function overwrite_ecephys_py_script(obj)
            
            template_script = obj.ctl.file.path_config.ecephys_template;
            dest_script = fullfile(obj.ctl.file.path_config.ecephys_scripts_dir, 'spikeGLX_pipeline_margrie.py');
            
            animal_id = obj.ctl.animal_id_from_probe_id(obj.probe_id);
            
            % gather variables which we need to replace in the script
            log_name = sprintf('''%s_log.csv''', animal_id);
            npx_directory = sprintf('r''%s''', fullfile(obj.ctl.file.path_config.processed_probe_fast_dir, animal_id));
            run_specs = sprintf('[[''%s'', ''0'', ''0,0'', ''0'']]', obj.probe_id);
            catGT_dest = sprintf('r''%s''', obj.ctl.file.processed_output_dir_fast(obj.probe_id));
            json_directory = sprintf('r''%s''', obj.ctl.file.json_dir_fast(obj.probe_id));
            
            % make sure we print the \
            npx_directory = strrep(npx_directory, '\', '\\');
            catGT_dest = strrep(catGT_dest, '\', '\\');
            json_directory = strrep(json_directory, '\', '\\');
            
            % open the default pipeline script
            fid = fopen(template_script, 'r');
            str = fread(fid, inf, '*char')';
            fclose(fid);
            
            % change some of the variables
            str = regexprep(str, '\<logName =[^\n]*\n',         sprintf('logName = %s\n', log_name));
            str = regexprep(str, '\<npx_directory =[^\n]*\n',   sprintf('npx_directory = %s\n', npx_directory));
            str = regexprep(str, '\nrun_specs =[^#]*',          sprintf('run_specs = %s\n\n', run_specs));
            str = regexprep(str, '\<catGT_dest =[^\n]*\n',      sprintf('catGT_dest = %s\n', catGT_dest));
            str = regexprep(str, '\<json_directory =[^\n]*\n',  sprintf('json_directory = %s\n', json_directory));
            
            % write a new pipeline script
            fid = fopen(dest_script, 'w');
            fprintf(fid, '%s', str);
            fclose(fid);
        end
        
        
        
        function run_ecephys(obj)
            
            fprintf('Running ecephys_spike_sorting...');
            
            % run the python pipeline
            script = fullfile(obj.ctl.file.path_config.ecephys_scripts_dir, 'spikeGLX_pipeline_margrie.py');
            if obj.leave_window_open_on_error
                cmd = sprintf('start /wait cmd /k %s %s', obj.ctl.file.path_config.ecephys_python_exe, script);
            else
                cmd = sprintf('start /wait cmd /c %s %s', obj.ctl.file.path_config.ecephys_python_exe, script);
            end
            system(cmd);
            % cleanup
            delete('C_Waves.log');
            delete('CatGT.log');
            
            fprintf('done\n');
        end
        
        
        
        function move_csvs(obj)
            
            fprintf('Moving csvs\n')
            
            ks2_dir = obj.ctl.file.imec0_ks2(obj.probe_id);
            csv_dir = obj.ctl.file.imec0_ks2_csv_dir(obj.probe_id);
            
            if ~isfolder(csv_dir)
                mkdir(csv_dir);
            end
            
            % files to move
            to_move = {'metrics.csv', 'overlap_summary.csv', 'waveform_metrics.csv', 'waveform_metrics_1.csv'};
            
            % move the directories
            for i = 1 : length(to_move)
                cmd = sprintf('move %s %s', fullfile(ks2_dir, to_move{i}), fullfile(csv_dir, to_move{i}));
                system(cmd);
            end
            
            copyfile(fullfile(ks2_dir, 'cluster_group.tsv'), fullfile(ks2_dir, 'cluster_groups.csv'));
        end
        
        
        
        function move_chanmap(obj)
        %%moves the channel map created by the janelia pipeline to the
        %%original .bin directory (required to fix waveforms) only need to
        %%do this for new NP2 probe
        
            if strcmp(obj.ctl.get_probe_type_from_metadata(obj.probe_id), '24')
                imec0_ks2 = obj.ctl.file.imec0_ks2(obj.probe_id);
                top_dir = fileparts(imec0_ks2);
                old_chanmap_fname = fullfile(top_dir, sprintf('%s_g0_tcat.imec0.ap_chanMap.mat', obj.probe_id));
                
                % presumably doing this on the fast drive
                glx_dir = obj.ctl.file.glx_bin_dir_processed_fast(obj.probe_id);
                new_chanmap_fname = fullfile(glx_dir, sprintf('%s_g0_t0.imec0.ap_chanMap.mat', obj.probe_id));
                
                copyfile(old_chanmap_fname, new_chanmap_fname);
            end
        end
        
        
        
        function fix_waveforms(obj)
        %%fixes for the fact that the janelia version of
        %%ecephys_spike_sorting takes the mean waveform properties from the
        %%processed data, not the raw data. The processed data does seem to
        %%change the shapes of the waveforms somewhat so here we go back
        %%and take waveforms from the raw data instead (running similar
        %%code but on the raw data file)
            
            % we will create new json files
            json_dir  = obj.ctl.file.json_dir_fast(obj.probe_id);
            ori_input_json_fname = fullfile(json_dir, [obj.probe_id, '_imec0-input.json']);
            new_input_json_fname = fullfile(json_dir, [obj.probe_id, '_imec0-fix_mean_waveforms_input.json']);
            new_output_json_fname = fullfile(json_dir, [obj.probe_id, '_imec0-fix_mean_waveforms_output.json']);
            
            % we will not overwrite the mean_waveforms.npy and
            % waveform_metrics.csv, but create new files
            new_mean_waveforms_fname = fullfile(obj.ctl.file.imec0_ks2(obj.probe_id), 'mean_waveforms_fix.npy');
            new_waveform_metrics_fname = fullfile(obj.ctl.file.imec0_ks2_csv_dir(obj.probe_id), 'waveform_metrics_fix.csv');
            
            % the raw data
            ap_band_fname = obj.ctl.file.glx_ap_bin_processed_fast(obj.probe_id);
            
            % replace backslashes with four backslashes (which will
            % become two after MATLAB's sprintf, and eventually 1)
            ap_band_fname = strrep(ap_band_fname, '\', '\\\\');
            new_mean_waveforms_fname = strrep(new_mean_waveforms_fname, '\', '\\\\');
            new_waveform_metrics_fname = strrep(new_waveform_metrics_fname, '\', '\\\\');
            
            % open the original input json
            fid = fopen(ori_input_json_fname, 'r');
            str = fread(fid, inf, '*char')';
            fclose(fid);
            
            % modify a few lines of the json
            str = regexprep(str, '\<"ap_band_file":[^\n]*\n', sprintf('"ap_band_file": "%s",\n', ap_band_fname));
            str = regexprep(str, '\<"mean_waveforms_file":[^\n]*\n', sprintf('"mean_waveforms_file": "%s",\n', new_mean_waveforms_fname));
            str = regexprep(str, '\<"waveform_metrics_file":[^\n]*\n', sprintf('"waveform_metrics_file": "%s"\n', new_waveform_metrics_fname));
%             str = regexprep(str, '\<"use_C_Waves":[^\n]*\n', '"use_C_Waves": false\n');
            
            % write to the new input json
            fid = fopen(new_input_json_fname, 'w');
            fprintf(fid, '%s', str);
            fclose(fid);
            
            % go
            if obj.leave_window_open_on_error
                cmd = sprintf('start /wait cmd /k %s -m ecephys_spike_sorting.modules.mean_waveforms --input_json "%s" --output_json "%s"', ...
                            obj.ctl.file.path_config.ecephys_python_exe, new_input_json_fname, new_output_json_fname);
            else
                cmd = sprintf('start /wait cmd /c %s -m ecephys_spike_sorting.modules.mean_waveforms --input_json "%s" --output_json "%s"', ...
                            obj.ctl.file.path_config.ecephys_python_exe, new_input_json_fname, new_output_json_fname);
            end
            system(cmd);
        end
    end
end
