classdef Git
% Git Class for handling information about git
%
%  Git Properties:
%    work_tree_dir  - full path to the working directory which is being tracked
%    git_dir        - full path to .git directory which is tracking
%    available      - does `git_dir` exist
%    current_commit - which commit is HEAD pointing to
%    is_clean       - is the git directory clean or are there changes to be
%                     staged
%
%  Git Methods:
%    info           - return structure containing import information
%    save           - save the info structure to a text file
%
%   Note: currently this class assumes that the git directory doing the
%   tracking is contained in the working directory as `.git`.

    properties (SetAccess = private)
        
        work_tree_dir
        git_dir
        available
        current_commit
        is_clean
    end
    
    
    
    methods
        
        function obj = Git(work_tree_dir)
        %%Git
        %
        %   Git(WORK_TREE) creates Git object where WORK_TREE is the full
        %   path to the working directory which is being tracked.
        %
        %   Currently this class assumes that the git directory doing the
        %   tracking is contained in the working directory as .git.
        
            obj.work_tree_dir = work_tree_dir;
            obj.git_dir = fullfile(work_tree_dir, '.git');
        end
        
        
        
        function val = get.available(obj)
        %%available Property, whether `git_dir` exists.
        %
        %   VAL = available.  VAL is true if git exists on the system and
        %   the `git_dir` exists.
        
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
        %%current_commit Property, the commit hash to which HEAD points.
        %
        %   SHA1 = current_commit.  SHA1 is a string containing the
        %   40-length sha1 code of the commit.
        
            if obj.available
                [~, sha1] = system(sprintf('git -C "%s" rev-parse HEAD', obj.work_tree_dir));
                sha1 = sha1(1:end-1);
            else
                sha1 = [];
            end
        end
        
        
        
        function is_clean = get.is_clean(obj)
        %%is_clean Property, boolean, whether or not working directory is
        %%clean.
        %
        %   IS_CLEAN = is_clean.  True if working directory is clean, false
        %   otherwise.
        
            if obj.available
                [~, output] = system(sprintf('git -C "%s" status --porcelain', obj.work_tree_dir));
                is_clean = isempty(output);
            else 
                is_clean = 0;
            end
        end
        
        
        
        function info = info(obj)
        %%info Return information about current git setup
        %
        %   INFO = info() returns structure INFO with fields:
        %       available - is `git_dir` available
        %       date - the *current* date and time
        %       sha1 - the git commit hash to which HEAD points
        %       git_dir - the git directory doing the tracking
        %       git_clean - whether or not the working directory is clean
        %
        %  See also: available, current_commit, is_clean
        
            info.available = obj.available;
            info.date = datestr(now);
            info.sha1 = obj.current_commit;
            info.git_dir = obj.git_dir;
            info.git_clean = obj.is_clean;
        end
        
        
        
        function save(obj, fname, force_save, prefix, append)
        %%save Save git information 
        %
        %  save(FILENAME, FORCE, PREFIX, APPEND)
        %    Saves info from returned from `info` methods to a text file.
        %    
        %    Args:
        %       FILENAME - the full path to the text file in which to save
        %                  the information
        %       FORCE    - (optional) if false (default), user will be prompted
        %                   whether to overwrite an existing file. If true,
        %                   the file will be overwritted whether or not 
        %       PREFIX   - (optional) when the information is saved it is
        %                   saved as text in the form <PREFIX>.git_dir =
        %                   ...  By default, PREFIX is set to 'git'.
        %       APPEND   - (optional) whether to write to (overwrite) a
        %                   file with the text information or append the
        %                   text information to the end of a text file.
        %                   Default is false = overwrite a file.
        %
        %   Information about git is saved as text in the following format:
        %
        %       <PREFIX>.git_dir = /path/to/git_dir
        %       <PREFIX>.sha1 = a38b384ffab.....
        %       <PREFIX>.date = 01-Jan-2010 13:53:08
        %       <PREFIX>.git_clean = 1
        %
        %   See also: info
        
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
