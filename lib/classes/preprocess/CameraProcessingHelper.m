classdef CameraProcessingHelper < handle
    
    properties (SetAccess = private)
        
        ctl
        session_id
        
        local_dir
        camera_ids
    end
    
    
    
    methods
        
        function obj = CameraProcessingHelper(ctl, session_id)
            
            obj.ctl = ctl;
            obj.session_id = session_id;
        end
        
        
        
        function val = get.local_dir(obj)
            val = obj.ctl.file.camera_csv_dir_fast(obj.session_id);
        end
        
        
        
        function val = get.camera_ids(obj)
            val = obj.ctl.get_camera_ids(obj.session_id);
        end
        
        
        
        function run_from_raw(obj)
            
            obj.move_raw_to_local();
            obj.run_runningmouse();
            obj.transform_avi();
            obj.save_git();
        end
        
        
        
        function move_raw_to_local(obj)
            
            remote_dir = obj.ctl.file.raw_camera_dir(obj.session_id);
            
            % copy the files from remote to local
            cmd = sprintf('xcopy /E /D "%s" "%s"', remote_dir, [obj.local_dir, '\']);
            system(cmd);
        end
        
        
        
        function run_runningmouse(obj)
            
            % run the python pipeline
            for ii = 1 : length(obj.camera_ids)
                camera_fname = fullfile(obj.local_dir, sprintf('%s.avi', obj.camera_ids{ii}));
                
                if ii == length(obj.camera_ids)
                    cmd = sprintf('start /wait cmd /k %s %s -v %s --save-csv -o %s', ...
                        obj.ctl.file.path_config.runningmouse_python_exe, ...
                        obj.ctl.file.path_config.runningmouse_main_script, ...
                        camera_fname, obj.local_dir);
                else
                    cmd = sprintf('start cmd /k %s %s -v %s --save-csv -o %s', ...
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