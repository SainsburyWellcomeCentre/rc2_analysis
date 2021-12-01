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
            info.sha1 = obj.current_commit;
            info.git_dir = obj.git_dir;
            info.git_clean = obj.is_clean;
        end
        
        
        
        function save(obj, fname, force_save, prefix, append)
            
            VariableDefault('force_save', false);
            VariableDefault('prefix', 'git');
            VariableDefault('append', false);
            
            % check whether to save
            if ~force_save && ~append
                if isfile(fname)
                    print_fname = strrep(fname, '\', '\\');
                    msg = sprintf('%s already exists, overwrite (Y)?', print_fname);
                    user = input(msg, 's');
                    if ~strcmp(user, 'Y')
                        fprintf('Not saving...\n');
                        return
                    end
                end
            end
            
            % append to file
            if append
                fid = fopen(fname, 'a');
            else
                fid = fopen(fname, 'w');
            end
            
            info = obj.info();
            
            fprintf(fid, '%s.git_dir = %s\n', prefix, info.git_dir);
            fprintf(fid, '%s.sha1 = %s\n', prefix, info.sha1);
            fprintf(fid, '%s.date = %s\n', prefix, info.date);
            fprintf(fid, '%s.git_clean = %i\n', prefix, info.git_clean);
            
            fclose(fid);
        end
    end
end
