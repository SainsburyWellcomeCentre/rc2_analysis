classdef CSVManager < handle
% CSVManager Class for helping with creation and saving of .csvs
%
%  CSVManager Properties:
%    save_on  -  whether to save the .csv when `save` is called (default = false, don't save)
%    curr_dir - base directory in which the .csvs are to be saved (default = currenct directory upon object creation)
%    table    - MATLAB table which is going to be saved as a .csv
%
%  CSVManager Methods:
%    set_csv_fulldir - set `curr_dir` to a full path to a directory
%    set_csv_subdir  - reset `curr_dir` relative to existing `curr_dir`
%    create_table    - create a MATLAB table with data
%    save            - save the `table` as a .csv

    properties
        
        save_on = false
    end
    
    properties (SetAccess = private)
        
        curr_dir
        table = table([])
    end
    
    
    
    methods
        
        function obj = CSVManager()
        %%CSVManager
        %
        %  Create object and set `curr_dir` to present working directory
        
            obj.curr_dir = pwd;
        end
        
        
        
        function set_csv_fulldir(obj, fulldir)
        %%set_csv_fulldir Set `curr_dir` to a full pathname
        %
        %  set_csv_fulldir(FULLDIR) sets `curr_dir` to the pathname given
        %  in FULLDIR. If the directory doesn't exist and `save_on` =
        %  true, the directory will be created.
        
            obj.curr_dir = fulldir;
            if obj.save_on && ~isfolder(obj.curr_dir)
                mkdir(obj.curr_dir)
            end
        end
        
        
        
        function set_csv_subdir(obj, varargin)
        %%set_csv_subdir Set `curr_dir` relative to the existing `curr_dir`
        %
        %  set_csv_subdir(VARARGIN) resets `curr_dir` to a pathname
        %  relative to the current `curr_dir`. VARARGIN should be a
        %  comma seaprated sequence of strings: STRING1, STRING2, ...
        %
        %   e.g. if `curr_dir` = /current/path/, set_csv_subdir('one', 'two', 'three') 
        %   would reset `curr_dir` to be /current/path/one/two/three
        %
        %  if the new `curr_dir` doesn't exist and `save_on` = true, the
        %  directory will be created
        
            assert(all(cellfun(@ischar, varargin)), ...
                'Not all arguments are character strings');
            
            save_dir = obj.curr_dir;
            for i = 1 : length(varargin)
                save_dir = fullfile(save_dir, varargin{i});
            end
            
            obj.curr_dir = save_dir;
            
            if obj.save_on && ~isfolder(obj.curr_dir)
                mkdir(obj.curr_dir)
            end
        end
        
        
        
        function create_table(obj, varargin)
        %%create_table Creates a MATLAB table with data
        %
        %  create_table(VARARGIN) creates a MATLAB table in `table` with
        %  data provided in VARARGIN 
        %
        %  VARARGIN should be of the form STRING1, DATA1, STRING2, DATA2,
        %  ... where each STRING is the heading of the table and DATA are
        %  vectors or string cell arrays all of the same length
        %
        %   e.g. create_table('one', ones(10, 1), 'two', rand(10, 1)
        %   would create a table with two columns with headings 'one' and
        %   'two' and with 10 rows: 1's in the first column and random numbers
        %   in the second column
        
            assert(mod(length(varargin), 2) == 0, 'Wrong number of arguments');
            
            keys = varargin(1:2:end);
            vals = varargin(2:2:end);
            
            assert(all(cellfun(@ischar, keys)), ...
                'Not all keys are character strings');
            
            obj.table = table(vals{:}); %#ok<CPROPLC> % put data in the table
            obj.table.Properties.VariableNames = keys;  % rename the headings
        end
        
        
        function save(obj, fname)
        %%save Save the MATLAB table in `table` as csv
        %
        %   save(FILENAME) saves the data in `table` (using writetable) if
        %   `save_on` = true (otherwise, does nothing). FILENAME is a
        %   string (not a full path, and .csv extension will be added)
        %
        %     e.g if `curr_dir` = /current/path, save('my_file') would save
        %     the table to /current/path/my_file.csv
        
            if obj.save_on
                csv_fname = fullfile(obj.curr_dir, sprintf('%s.csv', fname));
                writetable(obj.table, csv_fname); 
            end
        end
    end
end
 