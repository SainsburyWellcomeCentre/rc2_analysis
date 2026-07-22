classdef SortingHelper < handle
% SortingHelper Helper class for running the SpikeInterface-based NP2 sorting pipeline.
%
%  SortingHelper Properties:
%       leave_window_open_on_error - true or false (default) whether to leave the Windows
%                                    command prompt open upon an error or not.
%                                    If true, the prompt stays open when the pipeline
%                                    finishes or errors (user must close manually);
%                                    useful for debugging.
%       run_catgt                  - true (default) or false, whether to run the CatGT
%                                    preprocessing step
%       run_tprime                 - true or false (default), whether to run TPrime
%                                    at the end of the pipeline
%       ctl                        - instance of RC2Preprocess
%       probe_id                   - string with the probe recording ID
%       run_script                 - session script (template filled with session paths)
%       template_script            - template script path (from si_np2_template in path_config)
%       python_exe                 - python executable in the spikeinterface conda env
%
%  SortingHelper Methods:
%       run_from_raw             - run full SpikeInterface sorting pipeline
%       create_output_dirs       - create the CatGT destination directory
%       overwrite_sorting_script - fill session paths into template and save as run_script
%       run_sorting              - run the SpikeInterface pipeline

    properties

        leave_window_open_on_error = true
        run_catgt  = true
        run_tprime = false
    end

    properties (SetAccess = private)

        ctl
        probe_id

        run_script
        template_script
        python_exe
    end



    methods

        function obj = SortingHelper(ctl, probe_id)
        %%SortingHelper
        %
        %   SortingHelper(CTL, PROBE_ID) prepares a helper object where CTL
        %   is an RC2Preprocess instance and PROBE_ID is a string with the
        %   probe recording ID.

            obj.ctl      = ctl;
            obj.probe_id = probe_id;
        end



        function run_from_raw(obj)
        %%run_from_raw Run the full SpikeInterface sorting pipeline
        %
        %   run_from_raw() runs the full pipeline:
        %       - creates the CatGT destination directory
        %       - fills session paths into the template script
        %       - runs the SpikeInterface pipeline
        %
        %   All output files (metrics.csv, waveform_metrics.csv,
        %   waveform_metrics_fix.csv, cluster_groups.csv, .npy files)
        %   are written directly by the Python script into the correct
        %   locations under imec{prb}_ks4/. No post-pipeline file moves
        %   are needed from MATLAB.

            obj.create_output_dirs();
            obj.overwrite_sorting_script();
            obj.run_sorting();
        end



        function fname = get.run_script(obj)
        %%session script path (template filled with this session's paths)
        %
        %   Written to a '_generated' subfolder of si_np2_scripts_dir, kept
        %   separate from the template so it is never mistaken for a second
        %   maintained script. This file is overwritten on every run and is
        %   not tracked in git (see .gitignore).

            fname = fullfile(obj.ctl.file.path_config.si_np2_scripts_dir, ...
                             '_generated', 'spikeGLX_pipeline_session.py');
        end



        function fname = get.template_script(obj)
        %%template script path (from si_np2_template in path_config)

            fname = obj.ctl.file.path_config.si_np2_template;
        end



        function fname = get.python_exe(obj)
        %%python executable path (spikeinterface conda env)

            fname = obj.ctl.file.path_config.si_np2_python_exe;
        end



        function create_output_dirs(obj)
        %%create_output_dirs Create the CatGT destination directory
        %
        %   create_output_dirs() creates the output directory for the
        %   pipeline if it does not already exist.

            output_dir = obj.ctl.file.processed_output_dir_fast(obj.probe_id);

            if isfolder(output_dir)
                fprintf('%s already exists\n', output_dir);
            else
                fprintf('Making %s\n', output_dir);
                mkdir(output_dir)
            end
        end



        function overwrite_sorting_script(obj)
        %%overwrite_sorting_script Fill session paths into template and save as run_script
        %
        %   overwrite_sorting_script() reads template_script, replaces the
        %   session-specific variables (logName, npx_directory, run_specs,
        %   catGT_dest, run_CatGT, runTPrime) and writes the result to
        %   run_script.

            animal_id = obj.ctl.animal_id_from_probe_id(obj.probe_id);

            % values to inject into the template
            log_name      = sprintf('''%s_log.csv''', animal_id);
            npx_directory = sprintf('r''%s''', ...
                fullfile(obj.ctl.file.path_config.processed_probe_fast_dir, animal_id));
            run_specs     = sprintf('[[''%s'', ''0'', ''0,0'', ''0'']]', obj.probe_id);
            catGT_dest    = sprintf('r''%s''', ...
                obj.ctl.file.processed_output_dir_fast(obj.probe_id));

            % escape backslashes for Python raw-string literals
            npx_directory = strrep(npx_directory, '\', '\\');
            catGT_dest    = strrep(catGT_dest,    '\', '\\');

            % read template
            fid = fopen(obj.template_script, 'r');
            str = fread(fid, inf, '*char')';
            fclose(fid);

            % ensure the '_generated' output folder exists
            generated_dir = fileparts(obj.run_script);
            if ~isfolder(generated_dir)
                mkdir(generated_dir)
            end

            % inject session-specific values
            str = regexprep(str, '\<logName =[^\n]*\n',       sprintf('logName = %s\n',       log_name));
            str = regexprep(str, '\<npx_directory =[^\n]*\n', sprintf('npx_directory = %s\n', npx_directory));
            str = regexprep(str, '\nrun_specs =[^#]*',        sprintf('\nrun_specs = %s\n\n', run_specs));
            str = regexprep(str, '\<catGT_dest =[^\n]*\n',    sprintf('catGT_dest = %s\n',    catGT_dest));

            % CatGT / TPrime switches as Python booleans
            bool_str   = {'False', 'True'};
            catgt_str  = bool_str{logical(obj.run_catgt)  + 1};
            tprime_str = bool_str{logical(obj.run_tprime) + 1};

            str = regexprep(str, '\<run_CatGT = \w+', sprintf('run_CatGT = %s', catgt_str));
            str = regexprep(str, '\<runTPrime = \w+',  sprintf('runTPrime = %s', tprime_str));

            % write session script
            fid = fopen(obj.run_script, 'w');
            fprintf(fid, '%s', str);
            fclose(fid);
        end



        function run_sorting(obj)
        %%run_sorting Runs the SpikeInterface sorting pipeline
        %
        %   run_sorting() starts run_script with python_exe in a Windows
        %   command prompt. If leave_window_open_on_error is true the prompt
        %   stays open when the pipeline finishes (or on error).

            fprintf('Running SpikeInterface pipeline...');

            if obj.leave_window_open_on_error
                cmd = sprintf('start /wait cmd /k %s %s', obj.python_exe, obj.run_script);
            else
                cmd = sprintf('start /wait cmd /c %s %s', obj.python_exe, obj.run_script);
            end
            system(cmd);

            % clean up log files left in the working directory by CatGT / C_Waves
            if isfile('C_Waves.log'), delete('C_Waves.log'); end
            if isfile('CatGT.log'),   delete('CatGT.log');   end

            fprintf('done\n');
        end
    end
end
