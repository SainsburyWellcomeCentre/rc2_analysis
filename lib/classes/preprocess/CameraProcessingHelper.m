classdef CameraProcessingHelper < handle
    
    properties (SetAccess = private)
        
        ctl
        session_id
    end
    
    
    
    methods
        
        function obj = CameraProcessingHelper(ctl, session_id)
            
            obj.ctl = ctl;
            obj.session_id = session_id;
        end
        
        
        
        function run_from_raw(obj)
            
            obj.move_raw_to_local();
            obj.run_runningmouse();
            obj.transform_avi();
        end
        
        
        
        function move_raw_to_local(obj)
            
            remote_dir = obj.ctl.file.raw_camera_dir(obj.session_id);
            local_dir = obj.ctl.file.camera_csv_dir_fast(obj.session_id);
            
            % copy the files from remote to local
            cmd = sprintf('xcopy /E /D "%s" "%s"', remote_dir, [local_dir, '\']);
            system(cmd);
        end
        
        
        
        function run_runningmouse(obj)
            
            local_dir = obj.ctl.file.camera_csv_dir_fast(obj.session_id);
            camera_ids = obj.ctl.get_camera_ids(obj.session_id);
            
            % run the python pipeline
            for ii = 1 : length(camera_ids)
                camera_fname = fullfile(local_dir, sprintf('%s.avi', camera_ids{ii}));
                
                if ii == length(camera_ids)
                    cmd = sprintf('start /wait cmd /k %s C:\\Users\\lee\\Documents\\mvelez\\atyson\\runningmouse\\difference_video_2\\main.py -v %s --save-csv -o %s', ...
                        obj.ctl.file.path_config.runningmouse_python_exe, ...
                        camera_fname, local_dir);
                else
                    cmd = sprintf('start cmd /k %s C:\\Users\\lee\\Documents\\mvelez\\atyson\\runningmouse\\difference_video_2\\main.py -v %s --save-csv -o %s', ...
                        obj.ctl.file.path_config.runningmouse_python_exe, ...
                        camera_fname, local_dir);
                end
                system(cmd);
            end
        end
        
        
        
        function transform_avi(obj)
            
            local_dir = obj.ctl.file.camera_csv_dir_fast(obj.session_id);
            camera_ids = obj.ctl.get_camera_ids(obj.session_id);
            
            % run the python pipeline
            for ii = 1 : length(camera_ids)
                
                camera_fname = fullfile(local_dir, sprintf('%s.avi', camera_ids{ii}));
                
                cmd = sprintf('start /wait cmd /c ffmpeg -i %s %s', camera_fname, strrep(camera_fname, '.avi', '.mp4'));
                system(cmd);
                
                % remove the large .avi
                delete(camera_fname);
            end
        end
    end
end