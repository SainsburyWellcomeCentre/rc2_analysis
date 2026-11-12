classdef SortingHelper < handle
% SortingHelper Helper class for running the SpikeInterface-based NP2 sorting pipeline.
%
%  SortingHelper Properties:
%       leave_window_open_on_error - true or false (default) whether to leave the Windows
%                                    command prompt open upon an error or not.
%                                    If true, the prompt stays open when the pipeline
%                                    finishes or errors (user must close manually);
%                                    useful for debugging.
%       start_step                 - 'preprocess' (default), 'kilosort4', 'postprocess' or
%                                    'bombcell': where in the pipeline to resume from
%                                    (see spikeGLX_pipeline_np2.py's User input section)
%       ctl                        - instance of RC2Preprocess
%       probe_id                   - string with the probe recording ID
%       run_script                 - session script (template filled with session paths)
%       template_script            - template script path (from si_np2_template in path_config)
%       python_exe                 - python executable in the spikeinterface conda env
%
%  SortingHelper Methods:
%       run_from_raw             - run full SpikeInterface sorting pipeline
%       create_output_dirs       - create the pipeline's output directory
%       overwrite_sorting_script - fill session paths into template and save as run_script
%       run_sorting              - run the SpikeInterface pipeline

    properties

        leave_window_open_on_error = true
        start_step = 'preprocess'
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
        %       - creates the output directory
        %       - fills session paths into the template script
        %       - runs the SpikeInterface pipeline
        %
        %   All output files (metrics.csv, waveform_metrics.csv,
        %   cluster_groups.csv, .npy files)
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
        %%create_output_dirs Create the pipeline's output directory
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
        %   output_dest, start_step) and writes the result to run_script.

            valid_start_steps = {'preprocess', 'kilosort4', 'postprocess', 'bombcell'};
            if ~ismember(obj.start_step, valid_start_steps)
                error('SortingHelper:overwrite_sorting_script:invalidStartStep', ...
                    'start_step must be one of: %s', strjoin(valid_start_steps, ', '));
            end

            animal_id = obj.ctl.animal_id_from_probe_id(obj.probe_id);

            % values to inject into the template
            log_name      = sprintf('''%s_log.csv''', animal_id);
            npx_directory = sprintf('r''%s''', ...
                fullfile(obj.ctl.file.path_config.processed_probe_fast_dir, animal_id));
            run_specs     = sprintf('[[''%s'', ''0'', ''0,0'', ''0'']]', obj.probe_id);
            output_dest   = sprintf('r''%s''', ...
                obj.ctl.file.processed_output_dir_fast(obj.probe_id));

            % escape backslashes for Python raw-string literals
            npx_directory = strrep(npx_directory, '\', '\\');
            output_dest   = strrep(output_dest,   '\', '\\');

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
            str = regexprep(str, '\<output_dest =[^\n]*\n',   sprintf('output_dest = %s\n',   output_dest));
            str = regexprep(str, '\<start_step = ''\w+''', sprintf('start_step = ''%s''', obj.start_step));

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

            % Exit-code detection cannot go through "cmd /k ... & exit
            % %errorlevel%": with /k, that whole quoted string is ONE
            % command line, so the exit runs immediately after the Python
            % call regardless of outcome -- it closes the window right away
            % on ANY exit (success or failure), defeating the entire point
            % of leave_window_open_on_error (previously found the hard way:
            % a Python-side failure closed the window before anything could
            % be read from it). Instead, have the Python process itself drop
            % a small marker file on successful completion; its absence
            % after the window closes means it failed or was interrupted,
            % without needing the window's own exit code at all.
            done_marker = [tempname(), '.done'];
            if isfile(done_marker), delete(done_marker); end
            marked_script = sprintf('%s & (if not errorlevel 1 echo. > "%s")', ...
                sprintf('%s %s', obj.python_exe, obj.run_script), done_marker);

            if obj.leave_window_open_on_error
                cmd = sprintf('start /wait cmd /k "%s"', marked_script);
            else
                cmd = sprintf('start /wait cmd /c "%s"', marked_script);
            end
            system(cmd);
            succeeded = isfile(done_marker);
            if isfile(done_marker), delete(done_marker); end

            if ~succeeded
                error('SortingHelper:PythonPipelineFailed', ...
                    ['SpikeInterface pipeline did not complete successfully ', ...
                     '(no completion marker was written). Check the Python ', ...
                     'traceback (rerun with leave_window_open_on_error=true ', ...
                     'to see it) before trusting any downstream step.']);
            end

            fprintf('done\n');
        end
    end
end
