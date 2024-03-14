classdef Saver < handle
% Saver Class for helping with saving different types of data on the system
%
%  Saver Properties:
%    overwrite                  - whether to force overwriting of files (default = false)
%    git                        - instance of class Git
%    file_manager               - instance of class FileManager
%
%  Saver Methods:
%    svm_table                  - save MATLAB table of stationary and motion information to .csv
%    offsets_table              - save MATLAB table of sample offsets for replay trials to .csv
%    tuning_curves              - save tuning curve information for clusters to .mat
%    formatted_data             - save formatted data structure to .mat
%    append_to_formatted_data   - append specific variable of formatted data to a .mat file
%    hf_power_figure            - save the HF power figure to .pdf
%    track_offset               - save the offset between ephys and anatomy L5 for a shank to .txt
%    hf_power_parameters        - save the HF power parameters to .mat
%    clusters_janelia_csv       - save the clusters which pass quality metric criteria to .csv
%    mua_clusters_janelia_csv   - save the clusters which pass quality metric criteria for MUA to .csv
%    selected_clusters_txt      - save the manually selected clusters to a .txt
%    selected_mua_clusters_txt  - save the manually selected MUA clusters to a .txt
%    trigger_mat                - save the sync trigger channel to a separate .mat
%    original_trigger_mat       - if trigger needs to be updated save original trigger channels to separate .mat
%    trigger_points_removed     - if updating trigger indicate which samples were removed
%    create_tracks_dir          - create a 'tracks' directory in the Kilosort directory
%    driftmap                   - save the driftmap to a .pdf
%    writetable                 - general function for saving table to .csv
%    savemat                    - general function for saving data to .mat
%    append_to_git_cfg          - append git information to a text file
%    check_save                 - check whether the user wants to overwrite

    properties
        
        overwrite = false
    end
    
    properties (SetAccess = private)
        
        git
        file_manager
    end
    
    
    
    methods
        
        function obj = Saver(file_manager)
        %%Saver
        %
        %  Saver(FILE_MANAGER) where FILE_MANAGER is an object of class
        %   FileManager and controls information about path and file
        %   names on the current setup.
        
            obj.file_manager = file_manager;
            obj.git = Git(file_manager.path_config.git_work_tree_dir);
        end
        
        
        
        function svm_table(obj, probe_id, tbl)
        %%svm_table Save MATLAB table of stationary and motion information to .csv
        %
        %  svm_table(PROBE_ID, TABLE) save the data in TABLE for probe
        %  recording with ID, PROBE_ID.
        
            fname = obj.file_manager.svm_table(probe_id);
            obj.writetable(fname, tbl);
        end
        
        
        
        function offsets_table(obj, probe_id, tbl)
        %%offsets_table Save MATLAB table of sample offsets for replay trials to .csv
        %
        %  offsets_table(PROBE_ID, TABLE) save the data in TABLE for probe
        %  recording with ID, PROBE_ID.
        
            fname = obj.file_manager.offsets_table(probe_id);
            obj.writetable(fname, tbl);
        end
        
        
        
        function tuning_curves(obj, probe_id, tbl_struct)
        %%tuning_curves Save tuning curve information for clusters to .mat
        %
        %  tuning_curves(PROBE_ID, STRUCT) save the data in STRUCT for probe
        %  recording with ID, PROBE_ID.
        
            fname = obj.file_manager.tuning_curves(probe_id);
            obj.savemat(fname, tbl_struct);
        end
        
        function tuning_curves_acceleration(obj, probe_id, tbl_struct)
        %%tuning_curves Save tuning curve information for clusters to .mat
        %
        %  tuning_curves(PROBE_ID, STRUCT) save the data in STRUCT for probe
        %  recording with ID, PROBE_ID.
        
            fname = obj.file_manager.tuning_curves_acceleration(probe_id);
            obj.savemat(fname, tbl_struct);
        end
        
        
        
        function formatted_data(obj, probe_id, formatted_data)
        %%formatted_data Save formatted data structure to .mat
        %
        %  formatted_data(PROBE_ID, STRUCT) save the data in STRUCT for probe
        %  recording with ID, PROBE_ID.
        
            fname = obj.file_manager.formatted_data(probe_id);
            obj.savemat(fname, formatted_data);
        end
        
        
        
        function append_to_formatted_data(obj, probe_id, var)
        %%append_to_formatted_data Append specific variable of formatted data to a .mat file
        %
        %  append_to_formatted_data(PROBE_ID, VARIABLES) save the data in
        %  the VARIABLES structure to the formatted data .mat for probe
        %  recording, PROBE ID.
        %
        %  VARIABLES should be a structure with fields matching those seen in the formatted data file.
        %   See XX for more information on the types of information in the
        %   formatted data structure.
        
            fname = obj.file_manager.formatted_data(probe_id);
            save(fname, '-append', '-struct', 'var');
        end
        
        
        
        function hf_power_figure(obj, probe_id, shank_id, h_fig)
        %%hf_power_figure Save the HF power figure to .pdf
        %
        %  hf_power_figure(PROBE_ID, SHANK_ID, FIGURE_HANDLE) save 
        %  the figure referred to by FIGURE_HANDLE for probe recording with ID, PROBE_ID.
        
            fname = obj.file_manager.hf_power_figure(probe_id, shank_id);
            
            if obj.check_save(fname)
                figure(h_fig);
                print(fname, '-bestfit', '-dpdf');
            end
        end
        
        
        
        function track_offset(obj, probe_id, shank_id, val)
        %%track_offset Save the offset between ephys and anatomy L5 for a shank to .txt
        %
        %  track_offset(PROBE_ID, SHANK_ID, VALUE) save 
        %  the offset value VALUE for probe recording with ID, PROBE_ID and
        %  shank SHANK_ID (which is an integer (zero-index) indicating
        %  which shank we are saving).
            
            fname = obj.file_manager.track_offset(probe_id, shank_id);
            
            if obj.check_save(fname)
                fid = fopen(fname, 'w');
                fprintf(fid, '%.0f', val);
                fclose(fid);
            end
        end
        
        
        
        function hf_power_parameters(obj, probe_id, shank_id, params)
        %%hf_power_parameters Save the HF power parameters to .mat
        %
        %  hf_power_parameters(PROBE_ID, SHANK_ID, PARAMS) save 
        %  the parameters used to create the HF power profile for probe recording 
        %  with ID, PROBE_ID and shank SHANK_ID (which is an integer 
        %  (zero-index) indicating which shank we are saving).
        
            fname = obj.file_manager.hf_power_parameters(probe_id, shank_id);
            obj.savemat(fname, params);
        end
        
        
        
        function clusters_janelia_csv(obj, probe_id, tbl)
        %%clusters_janelia_csv Save the clusters which pass quality metric criteria to .csv
        %
        %  clusters_janelia_csv(PROBE_ID, TABLE) save the data in TABLE for probe
        %  recording with ID, PROBE_ID.
        
            fname = obj.file_manager.clusters_janelia_csv(probe_id);
            obj.writetable(fname, tbl);
        end
        
        
        
        function mua_clusters_janelia_csv(obj, probe_id, tbl)
        %%mua_clusters_janelia_csv Save the clusters which pass quality
        %%metric criteria for MUA to .csv
        %
        %  mua_clusters_janelia_csv(PROBE_ID, TABLE) save the data in TABLE
        %  for probe recording with ID, PROBE_ID.
        
            fname = obj.file_manager.mua_clusters_janelia_csv(probe_id);
            obj.writetable(fname, tbl);
        end
        
        
        
        function selected_clusters_txt(obj, probe_id, cluster_ids)
        %%selected_clusters_txt Save the manually selected clusters to a .txt
        %
        %  selected_clusters_txt(PROBE_ID, CLUSTER_IDS) save the list of selected cluster IDs for probe
        %  recording with ID, PROBE_ID.
        
            fname = obj.file_manager.selected_clusters_txt(probe_id);
            if obj.check_save(fname)
                fid = fopen(fname, 'w');
                for i = 1 : length(cluster_ids)
                    fprintf(fid, '%i\r\n', cluster_ids(i));
                end
                fclose(fid);
            end
        end
        
        
        
        function selected_mua_clusters_txt(obj, probe_id, cluster_ids)
        %%selected_mua_clusters_txt Save the manually selected MUA clusters
        %%to a .txt 
        %
        %  selected_mua_clusters_txt(PROBE_ID, CLUSTER_IDS) save the list
        %  of selected MUA cluster IDs for probe recording with ID,
        %  PROBE_ID.
        
            fname = obj.file_manager.selected_mua_clusters_txt(probe_id);
            if obj.check_save(fname)
                fid = fopen(fname, 'w');
                for i = 1 : length(cluster_ids)
                    fprintf(fid, '%i\r\n', cluster_ids(i));
                end
                fclose(fid);
            end
        end
        
        
        
        function trigger_mat(obj, probe_id, trigger)
        %%trigger_mat Save the sync trigger channel to a separate .mat
        %
        %  trigger_mat(PROBE_ID, TRIGGER) save the trigger channel in
        %  TRIGGER for probe recording with ID, PROBE_ID.
        
            fname = obj.file_manager.trigger_mat(probe_id);
            if obj.check_save(fname)
                save(fname, 'trigger');
            end
        end
        
        
        
        function original_trigger_mat(obj, probe_id, trigger)
        %%original_trigger_mat If trigger needs to be updated save original trigger channels to separate .mat
        %
        %  original_trigger_mat(PROBE_ID, TRIGGER) backup the trigger channel in
        %  TRIGGER for probe recording with ID, PROBE_ID.
        
            fname = obj.file_manager.original_trigger_mat(probe_id);
            if obj.check_save(fname)
                save(fname, 'trigger');
            end
        end
        
        
        
        function trigger_points_removed(obj, probe_id, points)
        %%trigger_points_removed If updating trigger indicate which samples were removed
        %
        %  trigger_points_removed(PROBE_ID, POINTS) save the points which were
        %  removed from the trigger channel for probe recording with ID, PROBE_ID.
        %  POINTS is a N x 2 array which specify the samples which were
        %  removed.  [FROM1, TO1; FROM2, TO2; ...; FROMN, TON]
        
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
        %%create_tracks_dir Create a 'tracks' directory in the Kilosort directory
        %
        %  create_tracks_dir(PROBE_ID) 
        
            [dname, exists] = obj.file_manager.tracks_dir(probe_id);
            
            if exists; return; end
            [parent_dir, sub_dir] = fileparts(dname);
            mkdir(parent_dir, sub_dir);
        end
        
        
        
        function driftmap(obj, probe_id, h_fig)
        %%driftmap Save the driftmap to a .pdf
        %
        %  driftmap(PROBE_ID, FIGURE_HANDLE) save the driftmap figure
        %  referenced by FIGURE_HANDLE for probe recording PROBE_ID.
        
            figure(h_fig);  % make current figure
            fname = obj.file_manager.driftmap(probe_id);
            if obj.check_save(fname) 
                print(fname, '-bestfit', '-dpdf');
            end
        end
        
        
        
        function writetable(obj, fname, tbl)
        %%writetable General function for saving table to .csv
        %
        %  writetable(FILENAME, TABLE) save MATLAB table in TABLE to
        %  FILENAME.
        
            if obj.check_save(fname)
                writetable(tbl, fname);
                [pathname, filename] = fileparts(fname);
                cfg_fname = fullfile(pathname, [filename, '.cfg']);
                obj.git.save(cfg_fname);
            end
        end
        
        
        
        function savemat(obj, fname, struct)
        %%savemat General function for saving data to .mat
        %
        %  savemat(FILENAME, STRUCT) save fields in MATLAB structure STRUCT
        %  to FILENAME.
            
            if obj.check_save(fname)
                save(fname, '-struct', 'struct', '-v7.3');
                git = obj.git.info; %#ok<*PROPLC>
                save(fname, '-append', 'git');
            end
        end
        
        
        
        function append_to_git_cfg(obj, fname, force_save, prefix)
        %%append_to_git_cfg Append git information to a text file
        %
        %  append_to_git_cfg(FILENAME, FORCE_SAVE, PREFIX) append git
        %  information in Git object to text file FILENAME.
        %
        %  See also: Git.save
        
            obj.git.save(fname, force_save, prefix, true);
        end
        
        
        
        function go_ahead = check_save(obj, fname)
        %%check_save Check whether the user wants to overwrite
        %
        %  GO = check_save(FILENAME) Checks whether FILENAME exists and if
        %  it does, asks user whether they want to overwrite. Returns GO =
        %  true, if user agrees to overwrite, false otherwise.
        %
        %   If `overwrite` property is true, no checks are made and GO is true.
        
            go_ahead = true;
            
            if obj.overwrite
                return
            end
            
            if isfile(fname)
                user = input(sprintf('%s exists, overwrite (Y)?', fname), 's');
                go_ahead = strcmp(user, 'Y');
            end
        end
    end
end
