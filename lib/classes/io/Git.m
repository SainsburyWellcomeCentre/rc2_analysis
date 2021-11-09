classdef Git < handle
    
    properties (SetAccess = private)
        
        dir
        available
        current_commit
        is_clean
    end
    
    
    
    methods
        
        function obj = Git(main_dir)
            
            obj.dir = main_dir;
        end
        
        
        
        function val = get.available(obj)
            
            [exit_status, ~] = system('git --version');
            
            if exit_status
                warning('git is not available on system');
            end
            if ~isfolder(obj.dir)
                warning('specified git directory doesn''t exist');
            end
            
            % git exists and directory exists
            val = ~exit_status && isfolder(obj.dir);
        end
        
        
        
        function sha1 = get.current_commit(obj)
            
            if obj.available
                [~, sha1] = system(sprintf('git -C "%s" rev-parse HEAD', obj.dir));
                sha1 = sha1(1:end-1);
            else
                sha1 = [];
            end
        end
        
        
        
        function is_clean = get.is_clean(obj)
            
            if obj.available
                [~, output] = system(sprintf('git -C "%s" status --porcelain', obj.dir));
                is_clean = isempty(output);
            else 
                is_clean = 0;
            end
        end
        
        
        
        function info = info(obj)
            
            info.date = datestr(now);
            info.sha1 = obj.current_commit();
            info.git_dir = obj.dir;
            info.git_clean = obj.is_clean();
        end
    end
end
