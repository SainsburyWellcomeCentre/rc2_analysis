classdef CSVManager < handle
    
    
    properties
        ctl
        
        save_on = false
        curr_dir
        
        table = table([]);
    end
    
    methods
        
        function obj = CSVManager(ctl)
            
            obj.ctl = ctl;
            obj.curr_dir = ctl.file.path_config.summary_data_dir;
        end
        
        
        
        function set_csv_fulldir(obj, fulldir)
            obj.curr_dir = fulldir;
            if obj.save_on && ~isfolder(obj.curr_dir)
                mkdir(obj.curr_dir)
            end
        end
        
        
        
        function set_csv_subdir(obj, varargin)
            
            assert(all(cellfun(@ischar, varargin)), ...
                'Not all arguments are character strings');
            
            save_dir = obj.ctl.file.path_config.summary_data_dir;
            for i = 1 : length(varargin)
                save_dir = fullfile(save_dir, varargin{i});
            end
            
            obj.curr_dir = save_dir;
            
            if obj.save_on && ~isfolder(obj.curr_dir)
                mkdir(obj.curr_dir)
            end
        end
        
        
        
        function create_table(obj, varargin)
            
            assert(mod(length(varargin), 2) == 0, 'Wrong number of arguments');
            
            keys = varargin(1:2:end);
            vals = varargin(2:2:end);
            
            assert(all(cellfun(@ischar, keys)), ...
                'Not all keys are character strings');
            
            obj.table = table(vals{:}); %#ok<CPROPLC> % put data in the table
            obj.table.Properties.VariableNames = keys;  % rename the headings
        end
        
        
        function save(obj, fname)
            
            if obj.save_on
                csv_fname = fullfile(obj.curr_dir, sprintf('%s.csv', fname));
                writetable(obj.table, csv_fname); 
            end
        end
    end
end
 