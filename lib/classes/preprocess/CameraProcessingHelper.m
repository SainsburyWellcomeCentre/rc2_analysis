classdef CameraProcessingHelper < handle
% CameraProcessingHelper Class for helping with processing of camera data
%
%   CameraProcessingHelper Properties:
%       ctl             - instance of class RC2Preprocess
%       session_id      - ID of the session being processed
%       local_dir       - directory where processed camera data will be saved
%       camera_ids      - IDs of the camera recordings associated with the session
%
% CameraProcessingHelper Methods:
%       run_from_raw        - runs full camera processing pipeline
%       move_raw_to_local   - moves raw .avi files locally to `local_dir`
%       run_runningmouse    - runs the `runningmouse` pipeline for processing of videos
%       transform_avi       - transforms the .avi files to .mp4 for space
%       save_git            - saves associated .cfg with git information


    properties (SetAccess = private)
        
        ctl
        session_id
        
        local_dir
        camera_ids
    end
    
    
    
    methods
        
        function obj = CameraProcessingHelper(ctl, session_id)
        %%CameraProcessingHelper
        %
        %   CameraProcessingHelper(CTL, SESSION_ID) creates a object of
        %   this class. CTL is the object of type RC2Preprocess and
        %   SESSION_ID is a string with the name of a recording session.
        
            obj.ctl = ctl;
            obj.session_id = session_id;
        end
        
        
        
        function val = get.local_dir(obj)
        %%local directory where to save processed data
            val = obj.ctl.file.camera_csv_dir_fast(obj.session_id);
        end
        
        
        
        function val = get.camera_ids(obj)
        %%camera recording IDs associated with session ID
            val = obj.ctl.get_camera_ids(obj.session_id);
        end
        
        
        
        function run_from_raw(obj)
        %%run_from_raw Runs full camera processing pipeline
        %
        %   run_from_raw() runs the full camera processing pipeline incl.
        %       - moving raw .avi files locally
        %       - run `runningmouse` program
        %       - transforms .avi to .mp4 to save space
        %       - save git information to .cfg
        
            obj.move_raw_to_local();
            obj.run_runningmouse();
            obj.transform_avi();
            obj.save_git();
        end
        
        
        
        function move_raw_to_local(obj)
        %%move_raw_to_local Moves raw .avi files locally to `local_dir`
        %
        %   move_raw_to_local() moves the camera .avi files associated with
        %   the current session ID (`session_id`), to a local directory.
        %   See README.
        
            remote_dir = obj.ctl.file.raw_camera_dir(obj.session_id);
            
            % copy the files from remote to local
            cmd = sprintf('xcopy /E /D "%s" "%s"', remote_dir, [obj.local_dir, '\']);
            system(cmd);
        end
        
        
        
        function run_runningmouse(obj)
        %%run_runningmouse Runs the `runningmouse` pipeline for processing of videos
        %
        %   run_runningmouse() runs the running mouse 
        
            for ii = 1 : length(obj.camera_ids)
                camera_fname = fullfile(obj.local_dir, sprintf('%s.avi', obj.camera_ids{ii}));
                
                if ii == length(obj.camera_ids)
                    cmd = sprintf('start /wait cmd /c %s %s -v %s --save-csv -o %s', ...
                        obj.ctl.file.path_config.runningmouse_python_exe, ...
                        obj.ctl.file.path_config.runningmouse_main_script, ...
                        camera_fname, obj.local_dir);
                else
                    cmd = sprintf('start cmd /c %s %s -v %s --save-csv -o %s', ...
                        obj.ctl.file.path_config.runningmouse_python_exe, ...
                        obj.ctl.file.path_config.runningmouse_main_script, ...
                        camera_fname, obj.local_dir);
                end
                system(cmd);
            end
        end
        
        
        
        function transform_avi(obj)
            
            % run the python pipeline
            for ii = 1 : length(obj.camera_ids)
                
                camera_fname = fullfile(obj.local_dir, sprintf('%s.avi', obj.camera_ids{ii}));
                
                cmd = sprintf('start /wait cmd /c ffmpeg -i %s %s', camera_fname, strrep(camera_fname, '.avi', '.mp4'));
                system(cmd);
                
                % remove the large .avi
                delete(camera_fname);
            end
        end
        
        
        
        function save_git(obj)
            
            for ii = 1 : length(obj.camera_ids)
                
                cfg_fname = fullfile(obj.local_dir, sprintf('%s.cfg', obj.camera_ids{ii}));
                
                cmd = sprintf('git -C "%s" rev-parse --show-toplevel', fileparts(obj.ctl.file.path_config.runningmouse_main_script));
                [~, output] = system(cmd);
                
                if ispc
                    output = strrep(output, '/', '\');
                end
                output = output(1:end-1);
                
                rm_git = Git(output);
                rm_git.save(cfg_fname, [], 'runningmouse_git');
                
                obj.ctl.save.append_to_git_cfg(cfg_fname, [], 'lib_git');
            end
        end
    end
end