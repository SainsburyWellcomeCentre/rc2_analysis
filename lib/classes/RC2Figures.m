classdef RC2Figures < handle
    
    properties
        
        config
        save_on = false
        curr_dir
        n_figs_to_join = 0
        fnames_to_join = {}
    end
    
    methods
        
        function obj = RC2Figures(config)
            
            obj.config = config;
            obj.curr_dir = config.figure_dir;
        end
        
        
        function h_fig = a4figure(obj, orientation)
            
            VariableDefault('orientation', 'portrait');
            h_fig = a4figure(orientation);
            
        end
        
        
        function h_fig = a3figure(obj, orientation)
            
            VariableDefault('orientation', 'portrait');
            h_fig = a3figure(orientation);
        end
        
        
        function set_figure_subdir(obj, varargin)
            
            assert(all(cellfun(@ischar, varargin)), ...
                'Not all arguments are character strings');
            
            save_dir = obj.config.figure_dir;
            for i = 1 : length(varargin)
                save_dir = fullfile(save_dir, varargin{i});
            end
            
            obj.curr_dir = save_dir;
            
            if obj.save_on && ~isfolder(obj.curr_dir)
                mkdir(obj.curr_dir)
            end
            
        end
        
        function save_fig_to_join(obj, use_opengl, dpi)
            
            VariableDefault('use_opengl', false);
            VariableDefault('dpi', 300);
            
            if obj.save_on
                obj.n_figs_to_join = obj.n_figs_to_join + 1;
                fname = sprintf('fig_%03i.pdf', obj.n_figs_to_join);
                
                obj.fnames_to_join{obj.n_figs_to_join} = fullfile(obj.curr_dir, fname);
                if use_opengl
                    print(obj.fnames_to_join{obj.n_figs_to_join}, '-dpdf', '-opengl', sprintf('-r%i', dpi));
                else
                    print(obj.fnames_to_join{obj.n_figs_to_join}, '-bestfit', '-dpdf', '-painters');
                end
            end
        end
        
        
        function join_figs(obj, fname)
            
            if obj.save_on
                fname = fullfile(obj.curr_dir, fname);
                join_pdfs(obj.fnames_to_join, fname, true);
            end
            
        end
        
        
        function save_fig(obj, fname, use_opengl, dpi)
            
            VariableDefault('use_opengl', false);
            VariableDefault('dpi', 300);
            
            if obj.save_on
                fname = fullfile(obj.curr_dir, fname);
                if use_opengl
                    print(fname, '-dpdf', '-opengl', sprintf('-r%i', dpi));
                else
                    print(fname, '-dpdf', '-painters');
                end
            end
        end
        
        
        function clear_figs(obj)
            
            obj.n_figs_to_join = 0;
            obj.fnames_to_join = {};
            close all
        end
        
    end
    
end