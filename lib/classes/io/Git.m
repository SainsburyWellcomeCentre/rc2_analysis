classdef Git < handle
    
    properties (SetAccess = private)
        
        work_tree_dir
        git_dir
        available
        current_commit
        is_clean
    end
    
    
    
    methods
        
        function obj = Git(work_tree_dir)
            
            obj.work_tree_dir = work_tree_dir;
            obj.git_dir = fullfile(work_tree_dir, '.git');
        end
        
        
        
        function val = get.available(obj)
            
            [exit_status, ~] = system('git --version');
            
            if exit_status
                warning('git is not available on system');
            end
            if ~isfolder(obj.git_dir)
                warning('specified git directory doesn''t exist');
            end
            
            % git exists and directory exists
            val = ~exit_status && isfolder(obj.git_dir);
        end
        
        
        
        function sha1 = get.current_commit(obj)
            
            if obj.available
                [~, sha1] = system(sprintf('git -C "%s" rev-parse HEAD', obj.work_tree_dir));
                sha1 = sha1(1:end-1);
            else
                sha1 = [];
            end
        end
        
        
        
        function is_clean = get.is_clean(obj)
            
            if obj.available
                [~, output] = system(sprintf('git -C "%s" status --porcelain', obj.work_tree_dir));
                is_clean = isempty(output);
            else 
                is_clean = 0;
            end
        end
        
        
        
        function info = info(obj)
            
            info.available = obj.available;
            info.date = datestr(now);
            info.sha1 = obj.current_commit();
            info.git_dir = obj.git_dir;
            info.git_clean = obj.is_clean();
        end
    end
end
