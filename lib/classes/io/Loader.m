classdef Loader < handle
    
    properties
        
        file_manager
    end
    
    
    methods
        
        function obj = Loader(file_manager)
            
            obj.file_manager = file_manager;
        end
        
        
        
        function [data, dt, chan_names, config] = rc2_bin(obj, session_id)
            
            [fname, exists] = obj.file_manager.rc2_bin(session_id);
            if exists
                [data, dt, chan_names, config] = read_rc2_bin(fname);
            else
                data = []; dt = []; chan_names = []; config = [];
            end
        end
        
        
        
        function camera = camera_csv(obj, session_id, camera_id)
            
            [fname, exist] = obj.file_manager.camera_csv_fast(session_id, camera_id);
            if ~exist
                [fname, exist] = obj.file_manager.camera_csv_slow(session_id, camera_id);
            end
            if ~exist
                camera = [];
                return
            end
            tbl = obj.readtable(fname);
            camera = tbl.Var1;
        end
        
        
        
        function data = formatted_data(obj, probe_id)
            
            fname = obj.file_manager.formatted_data(probe_id);
            data = obj.loadmat(fname);
        end
        
        
        
        function data = load_from_formatted_data(obj, probe_id, var)
        %%load a specific variable from the formatted data
            if ~iscell(var)
                var = {var};
            end
            fname = obj.file_manager.formatted_data(probe_id);
            data = load(fname, var{:});
        end
        
        
        
        function e_list = experiment_list(obj)
            
            fname = obj.file_manager.experiment_list();
            e_list = obj.readtable(fname);
        end
        
        
        
        function tbl = svm_table(obj, probe_id)
            
            fname = obj.file_manager.svm_table(probe_id);
            tbl = obj.readtable(fname);
        end
        
        
        
        function tbl = offsets_table(obj, probe_id)
            
            fname = obj.file_manager.offsets_table(probe_id);
            tbl = obj.readtable(fname);
        end
        
        
        
        function tuning_curves = tuning_curves(obj, probe_id)
            
            fname = obj.file_manager.tuning_curves(probe_id);
            load(fname, 'tuning_curves');
        end
        
        
        
        function track = track_csv(obj, probe_id, shank_id)
            
            fname = obj.file_manager.track_csv(probe_id, shank_id);
            warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
            track = obj.readtable(fname);
            warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
        end
        
        
        
        function offset = track_offset(obj, probe_id, shank_id)
            
            [fname, exists] = obj.file_manager.track_offset(probe_id, shank_id);
            
            offset = [];
            if exists
                fid = fopen(fname, 'r');
                offset = textscan(fid, '%f%[^\n\r]');
                offset = offset{1};
                fclose(fid);
            end
        end
        
        
        
        function val = ks2_npy(obj, probe_id, var)
            
            fname = obj.file_manager.ks2_npy(probe_id, var);
            val = obj.readnpy(fname);
        end
        
        
        
        function tbl = cluster_groups(obj, probe_id)
            
            fname = obj.file_manager.cluster_groups(probe_id);
            tbl = obj.readtable(fname);
        end
        
        
        
        function tbl = ks2_label(obj, probe_id)
            
            fname = obj.file_manager.ks2_label(probe_id);
            tbl = [];
            if isfile(fname)
                tbl = readtable(fname, 'FileType', "delimitedtext");
            end
        end
        
        
        
        function trigger = trigger_mat(obj, probe_id)
            
            fname = obj.file_manager.trigger_mat(probe_id);
            mat = obj.loadmat(fname);
            trigger = mat.trigger;
        end
        
        
        
        function data = params(obj, probe_id)
            
            [fname, exists] = obj.file_manager.params(probe_id);
            
            data = []; 
            
            if ~exists; return; end
            
            format_spec = '%s%s%s%[^\n\r]';
            fid = fopen(fname,'r');
            text_str = textscan(fid, format_spec, inf, ...
                'delimiter', ' ', ...
                'multipledelimsasone', true, ...
                'texttype', 'string', ...
                'headerlines', 0, ...
                'returnonerror', false, ...
                'endofline', '\r\n');
            
            for i = 1 : length(text_str{1})
                data.(text_str{1}{i}) = text_str{3}{i};
            end
            
            data.n_channels_dat     = str2double(data.n_channels_dat);
            data.offset             = str2double(data.offset);
            data.sample_rate        = str2double(data.sample_rate);
            
            fclose(fid);
        end
        
        
        
        function tbl = metrics_csv(obj, probe_id)
            
            fname = obj.file_manager.metrics_csv(probe_id);
            tbl = obj.readtable(fname);
        end
        
        
        
        function tbl = waveform_metrics_csv(obj, probe_id)
            
            fname = obj.file_manager.waveform_metrics_csv(probe_id);
            tbl = obj.readtable(fname);
        end
        
        
        
        function tbl = waveform_metrics_fixed_csv(obj, probe_id)
            
            fname = obj.file_manager.waveform_metrics_fixed_csv(probe_id);
            tbl = obj.readtable(fname);
        end
        
        
        
        function tbl = clusters_janelia_xlsx(obj, probe_id)
            
            fname = obj.file_manager.clusters_janelia_xlsx(probe_id);
            tbl = readtable(fname);
        end
        
        
        
        function cluster_ids = selected_clusters(obj, probe_id)
            
            fname = obj.file_manager.selected_clusters_txt(probe_id);
            
            format_spec = '%f%[^\n\r]';
            fid = fopen(fname, 'r');
            data = textscan(fid, format_spec, ...
                            'texttype', 'string', ...
                            'headerlines', 0, ...
                            'returnonerror', false);
            fclose(fid);
            cluster_ids = data{1};
        end
        
        
        
        function meta = spikeglx_ap_metadata(obj, probe_id)
            
            % search first for fast, then for slow drives
            [dname, exists] = obj.file_manager.glx_bin_dir_processed_fast(probe_id);
            
            if ~exists
                [dname, exists] = obj.file_manager.glx_bin_dir_processed_slow(probe_id);
            end
            
            if ~exists
                meta = [];
                return
            end
            
            meta = SpikeGLXMetaData(dname, 'ap');
        end
        
        
        
        function rec = spikeglx_ap_recording(obj, probe_id)
            
            % search first for fast, then for slow drives
            [dname, exists] = obj.file_manager.glx_bin_dir_processed_fast(probe_id);
            
            if ~exists
                [dname, exists] = obj.file_manager.glx_bin_dir_processed_slow(probe_id);
            end
            
            if ~exists
                rec = [];
                return
            end
            
            rec = SpikeGLXRecording(dname, 'ap');
        end
    end
    
    
    
    methods (Static = true)
        
        function tbl = readtable(fname)
            tbl = [];
            if isfile(fname)
                tbl = readtable(fname);
            end
        end
        
        
        
        function struct = loadmat(fname)
            struct = [];
            if isfile(fname)
                struct = load(fname);
            end
        end
        
        
        
        function val = readnpy(fname)
            val = [];
            if isfile(fname)
                val = readNPY(fname);
            end
        end
    end
end
