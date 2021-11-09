classdef Saver < handle
    
    properties (SetAccess = private)
        
        git
        file_manager
    end
    
    
    
    methods
        
        function obj = Saver(file_manager)
        %%class for saving files
            obj.file_manager = file_manager;
            obj.git = Git(file_manager.path_config.git_dir);
        end
        
        
        
        function svm_table(obj, probe_id, tbl)
        
            fname = obj.file_manager.svm_table(probe_id);
            obj.writetable(fname, tbl);
        end
        
        
        
        function offsets_table(obj, probe_id, tbl)
            
            fname = obj.file_manager.offsets_table(probe_id);
            obj.writetable(fname, tbl);
        end
        
        
        
        function tuning_table(obj, probe_id, tbl_struct)
            
            fname = obj.file_manager.tuning_table(probe_id);
            obj.savemat(fname, tbl_struct);
        end
        
        
        
        function formatted_data(obj, probe_id, formatted_data)
            
            fname = obj.file_manager.formatted_data(probe_id);
            obj.savemat(fname, formatted_data);
        end
        
        
        
        function append_to_formatted_data(obj, probe_id, var)
        %%appends a variable to the formatted data mat file
        %   var should be a structure with the fields to append
            fname = obj.file_manager.formatted_data(probe_id);
            save(fname, '-append', '-struct', 'var');
        end
        
        
        
        function hf_power_figure(obj, probe_id, shank_id, h_fig)
            
            fname = obj.file_manager.hf_power_figure(probe_id, shank_id);
            
            if obj.check_save(fname)
                figure(h_fig);
                print(fname, '-bestfit', '-dpdf');
            end
        end
        
        
        
        function track_offset(obj, probe_id, shank_id, val)
            
            fname = obj.file_manager.track_offset(probe_id, shank_id);
            
            if obj.check_save(fname)
                fid = fopen(fname, 'w');
                fprintf(fid, '%.0f', val);
                fclose(fid);
            end
        end
        
        
        
        function hf_power_parameters(obj, probe_id, shank_id, params)
            
            fname = obj.file_manager.hf_power_parameters(probe_id, shank_id);
            obj.savemat(fname, params);
        end
        
        
        
        function clusters_janelia_csv(obj, probe_id, tbl)
            
            fname = obj.file_manager.clusters_janelia_csv(probe_id);
            obj.writetable(fname, tbl);
        end
        
        
        
        function selected_clusters_txt(obj, probe_id, cluster_ids)
            
            % write to file
            fname = obj.file_manager.selected_clusters_txt(probe_id);
            if obj.check_save(fname)
                fid = fopen(fname, 'w');
                for i = 1 : length(cluster_ids)
                    fprintf(fid, '%i\r\n', cluster_ids(i));
                end
                fclose(fid);
            end
        end
        
        
        
        function trigger_mat(obj, probe_id, trigger)
            
            fname = obj.file_manager.trigger_mat(probe_id);
            if obj.check_save(fname)
                save(fname, 'trigger');
            end
        end
        
        
        
        function original_trigger_mat(obj, probe_id, trigger)
            
            fname = obj.file_manager.original_trigger_mat(probe_id);
            if obj.check_save(fname)
                save(fname, 'trigger');
            end
        end
        
        
        
        function trigger_points_removed(obj, probe_id, points)
        %%save the points we have removed from the trigger
            fname = obj.file_manager.trigger_points_removed(probe_id);
            if obj.check_save(fname)
                fid = fopen(fname, 'w');
                for ii = 1 : size(points, 1)
                    fprintf(fid, '%i,%i\r\n', points(ii, 1), points(ii, 2));
                end
                fclose(fid);
            end
        end
        
        
        
        function create_tracks_dir(obj, probe_id)
            
            [dname, exists] = obj.file_manager.tracks_dir(probe_id);
            
            if exists; return; end
            [parent_dir, sub_dir] = fileparts(dname);
            mkdir(parent_dir, sub_dir);
        end
        
        
        
        function driftmap(obj, probe_id, h_fig)
            
            figure(h_fig);  % make current figure
            fname = obj.file_manager.driftmap(probe_id);
            if obj.check_save(fname)
                print(fname, '-bestfit', '-dpdf');
            end
        end
        
        
        
        function writetable(obj, fname, tbl)
            
            if obj.check_save(fname)
                writetable(tbl, fname);
            end
        end
        
        
        
        function savemat(obj, fname, struct)
        %%takes a filename 'fname' and a MATLAB structure 'struct'
        %   will save the fields of the structure in the filename
            
            if obj.check_save(fname)
                save(fname, '-struct', 'struct', '-v7.3');
                git = obj.git.info; %#ok<*PROPLC>
                save(fname, '-append', 'git');
            end
        end
    end
    
    
    
    methods (Static = true)
        
        function go_ahead = check_save(fname)
            
            go_ahead = true;
            if isfile(fname)
                user = input('File exists, overwrite (Y)?', 's');
                go_ahead = strcmp(user, 'Y');
            end
        end
    end
end
