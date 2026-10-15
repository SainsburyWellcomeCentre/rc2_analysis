classdef PupilProcessingHelper < handle
% PupilProcessingHelper Class for computing pupil diameter and position
% from the DeepLabCut tracking of camera0 (contralateral eye)
%
%   PupilProcessingHelper Properties:
%       ctl                 - instance of class RC2Preprocess
%       session_id          - ID of the session being processed
%       likelihood_cutoff   - DLC points with likelihood below this value are set to NaN
%       local_dir           - directory containing the DLC .csv, where the outputs are also saved
%
%   PupilProcessingHelper Methods:
%       run_from_dlc        - compute and save pupil diameter and position
%       find_dlc_csv        - find the DLC .csv for camera0 of the session
%       load_dlc            - load the x, y and likelihood of each DLC body part
%       compute             - compute pupil diameter and position from the DLC points
%       save_cfg            - save .cfg with the DLC file, cutoff and git information
%
%   The DLC .csv is expected in path_config.processed_camera_fast_dir\<SESSION_ID>\
%   with a name of the form sub-*_cam0DLC*.csv, and must contain the body
%   parts pupil_top, pupil_bottom, pupil_nasal and pupil_temporal.
%
%   Pupil diameter (pixels) is the mean of the top-bottom and
%   nasal-temporal distances, using whichever of the two is available on
%   each frame. Pupil position (pixels, image coordinates) is the mean of
%   the top-bottom and nasal-temporal midpoints, again using whichever is
%   available. Frames where neither is available are NaN.


    properties (SetAccess = private)

        ctl
        session_id

        likelihood_cutoff = 0.7

        local_dir
    end

    properties (Constant = true)

        body_parts = {'pupil_top', 'pupil_bottom', 'pupil_nasal', 'pupil_temporal'}
    end



    methods

        function obj = PupilProcessingHelper(ctl, session_id)
        %%PupilProcessingHelper
        %
        %   PupilProcessingHelper(CTL, SESSION_ID) creates a object of
        %   this class. CTL is the object of type RC2Preprocess and
        %   SESSION_ID is a string with the name of a recording session.

            obj.ctl = ctl;
            obj.session_id = session_id;
        end



        function val = get.local_dir(obj)
        %%local directory containing the DLC output
            val = obj.ctl.file.camera_csv_dir_fast(obj.session_id);
        end



        function run_from_dlc(obj)
        %%run_from_dlc Compute and save pupil diameter and position
        %
        %   run_from_dlc() finds the DLC .csv for the session, computes
        %   the pupil diameter and position and saves them to
        %   pupil_diameter.csv and pupil_position.csv in `local_dir`. If
        %   no DLC .csv is found a message is printed and nothing is saved.

            dlc_fname = obj.find_dlc_csv();
            if isempty(dlc_fname)
                return
            end

            pts = obj.load_dlc(dlc_fname);
            [diameter, position] = obj.compute(pts);

            n_frames = length(diameter);
            fprintf('[pupil] %s: %i frames, %.1f%% NaN diameter, %.1f%% NaN position (likelihood cutoff %.2f)\n', ...
                obj.session_id, n_frames, 100*mean(isnan(diameter)), ...
                100*mean(isnan(position(:, 1))), obj.likelihood_cutoff);

            diameter_fname = obj.ctl.file.pupil_diameter_fast(obj.session_id);
            writetable(table(diameter, 'VariableNames', {'pupil_diameter'}), diameter_fname);

            position_fname = obj.ctl.file.pupil_position_fast(obj.session_id);
            writetable(table(position(:, 1), position(:, 2), 'VariableNames', {'x', 'y'}), position_fname);

            obj.save_cfg(dlc_fname);

            fprintf('[pupil] %s: saved %s and %s\n', obj.session_id, diameter_fname, position_fname);
        end



        function fname = find_dlc_csv(obj)
        %%find_dlc_csv Find the DLC .csv for camera0 of the session
        %
        %   FILENAME = find_dlc_csv() returns the full path to the DLC .csv
        %   for camera0 in `local_dir`. If there is none, FILENAME is empty
        %   and a message is printed. If there are several, the most recent
        %   is used and a message is printed.

            fname = '';

            if ~isfolder(obj.local_dir)
                fprintf('[pupil] %s: camera directory %s does not exist, skipping\n', ...
                    obj.session_id, obj.local_dir);
                return
            end

            contents = dir(fullfile(obj.local_dir, 'sub-*_cam0DLC*.csv'));

            if isempty(contents)
                fprintf('[pupil] %s: no DLC .csv (sub-*_cam0DLC*.csv) in %s, skipping\n', ...
                    obj.session_id, obj.local_dir);
                return
            end

            if length(contents) > 1
                [~, idx] = max([contents(:).datenum]);
                fprintf('[pupil] %s: %i DLC .csv files in %s, using the most recent: %s\n', ...
                    obj.session_id, length(contents), obj.local_dir, contents(idx).name);
                contents = contents(idx);
            end

            fname = fullfile(contents.folder, contents.name);
        end



        function pts = load_dlc(obj, fname)
        %%load_dlc Load the x, y and likelihood of each DLC body part
        %
        %   PTS = load_dlc(FILENAME) reads the DLC .csv FILENAME and
        %   returns a structure PTS with one field per body part in
        %   `body_parts`, each an N x 3 matrix of [x, y, likelihood] for N
        %   frames. The columns are found from the 'bodyparts' and 'coords'
        %   header rows of the .csv.

            fid = fopen(fname, 'r');
            fgetl(fid);  % scorer
            bodyparts = strsplit(fgetl(fid), ',');
            coords = strsplit(fgetl(fid), ',');
            fclose(fid);

            data = readmatrix(fname, 'NumHeaderLines', 3);

            for ii = 1 : length(obj.body_parts)

                this_part = strcmp(bodyparts, obj.body_parts{ii});

                x_col = find(this_part & strcmp(coords, 'x'));
                y_col = find(this_part & strcmp(coords, 'y'));
                l_col = find(this_part & strcmp(coords, 'likelihood'));

                assert(isscalar(x_col) && isscalar(y_col) && isscalar(l_col), ...
                    'Body part "%s" not found (or found more than once) in %s', ...
                    obj.body_parts{ii}, fname);

                pts.(obj.body_parts{ii}) = data(:, [x_col, y_col, l_col]);
            end
        end



        function [diameter, position] = compute(obj, pts)
        %%compute Compute pupil diameter and position from the DLC points
        %
        %   [DIAMETER, POSITION] = compute(PTS) takes the structure PTS
        %   returned by `load_dlc` and returns DIAMETER, an N x 1 vector of
        %   pupil diameter in pixels, and POSITION, an N x 2 matrix of
        %   [x, y] pupil centre in pixels. Points with likelihood below
        %   `likelihood_cutoff` are set to NaN.

            for ii = 1 : length(obj.body_parts)
                xy = pts.(obj.body_parts{ii})(:, 1:2);
                xy(pts.(obj.body_parts{ii})(:, 3) < obj.likelihood_cutoff, :) = NaN;
                xy_pts.(obj.body_parts{ii}) = xy;
            end

            % distances are NaN if either point of the pair is NaN
            top_bottom_d = vecnorm(xy_pts.pupil_top - xy_pts.pupil_bottom, 2, 2);
            nasal_temporal_d = vecnorm(xy_pts.pupil_nasal - xy_pts.pupil_temporal, 2, 2);

            diameter = mean([top_bottom_d, nasal_temporal_d], 2, 'omitnan');

            % midpoints are NaN if either point of the pair is NaN
            top_bottom_c = (xy_pts.pupil_top + xy_pts.pupil_bottom) / 2;
            nasal_temporal_c = (xy_pts.pupil_nasal + xy_pts.pupil_temporal) / 2;

            position = [mean([top_bottom_c(:, 1), nasal_temporal_c(:, 1)], 2, 'omitnan'), ...
                        mean([top_bottom_c(:, 2), nasal_temporal_c(:, 2)], 2, 'omitnan')];
        end



        function save_cfg(obj, dlc_fname)
        %%save_cfg Save .cfg with the DLC file, cutoff and git information
        %
        %   save_cfg(DLC_FILENAME) writes pupil.cfg to `local_dir`
        %   recording which DLC .csv and likelihood cutoff were used,
        %   followed by the git information of rc2_analysis.

            cfg_fname = fullfile(obj.local_dir, 'pupil.cfg');

            fid = fopen(cfg_fname, 'w');
            fprintf(fid, 'pupil.dlc_csv = %s\n', dlc_fname);
            fprintf(fid, 'pupil.likelihood_cutoff = %.2f\n', obj.likelihood_cutoff);
            fclose(fid);

            obj.ctl.save.append_to_git_cfg(cfg_fname, [], 'lib_git');
        end
    end
end
