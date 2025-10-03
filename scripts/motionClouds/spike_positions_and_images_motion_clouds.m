% Consolidated Motion Clouds spike export: positions and image links
%
% For experiment group 'passive_same_luminance_mc':
% - For each spike, optionally restrict to motion-only periods
% - Compute current position at spike time (cm) from AlignedTrial.position
% - Compute lagged image indices via lags_ms_multi (no fixed delay)
% - Map trial_id -> motion_cloud_id (via motion_cloud_sequence_250414.mat)
% - Map trial_id -> protocol_id (via passive_protocol_sequence_motion_clouds.mat)
% - Convert past_position to image_index assuming 1 image per 0.1 cm
% - Save one CSV per probe with all fields
%
% Output CSV columns:
%   probe_id, session_id, trial_id, protocol_id, cluster_id, spike_time,
%   position_cm, motion_cloud_id, image_rel_path, image_index_lagXXms

close all;

% User parameters
experiment_groups     = {'passive_same_luminance_mc'};
index_base            = 1;      % 1 for 1-based image indexing; 0 for 0-based
cfg = path_config();
images_root           = fullfile(cfg.motion_clouds_root, 'saved', 'motionClouds');
restrict_to_motion    = true;   % if true, only include spikes during motion masks
lags_ms_multi         = 40:10:120; % extra pre-spike lags for multi-lag indices

save_csv              = true;
save_subdirs          = {'motionClouds', 'passive_same_luminance_mc', 'spike_positions_and_images'};

% Initialize controller
ctl       = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});

% CSV Manager
cm = CSVManager();
cm.save_on = save_csv;
out_dir_base = cfg.figure_dir;
cm.set_csv_fulldir(out_dir_base);
cm.set_csv_subdir(save_subdirs{:});

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));

% Build Motion Cloud directory map (alphabetical order)
dir_info = dir(images_root);
names = {dir_info([dir_info.isdir]).name};
names = names(~ismember(names,{'.','..'}));
names = sort(names);

% Load image_folders.mat and REQUIRE exact match to alphabetical directory list
cloud_names = names; % initialize; will be replaced by saved order after validation
saved_file = fullfile(cfg.motion_clouds_root, 'image_folders.mat');
if ~exist(saved_file, 'file')
    error('image_folders.mat not found at %s', saved_file);
end

Ssaved = load(saved_file);
if isstruct(Ssaved) && isfield(Ssaved, 'image_folders') && isstruct(Ssaved.image_folders) && isfield(Ssaved.image_folders, 'name')
    names_struct = Ssaved.image_folders;
    saved_names = {names_struct.name};
else
    saved_names = extract_string_list_from_struct(Ssaved);
end

if isempty(saved_names)
    error('image_folders.mat present but no folder name list could be extracted.');
end

saved_names = saved_names(:)';
if numel(saved_names) ~= numel(names)
    error('Mismatch in count: %d folders under %s vs %d names in image_folders.mat', numel(names), images_root, numel(saved_names));
end

first_mismatch = find(~strcmp(saved_names, names), 1);
if ~isempty(first_mismatch)
    error('Order/name mismatch at index %d: dir="%s" vs mat="%s"', first_mismatch, names{first_mismatch}, saved_names{first_mismatch});
end

% If we reach here, alphabetical directory order matches saved order; use saved order explicitly
cloud_names = saved_names;
fprintf('Verified: alphabetical folders under %s match image_folders.mat order (field "name").\n', images_root);

% Validate that each cloud_names entry exists under images_root when used

% (Removed verification against image_folders.mat; using alphabetical directory order)

% Load Motion Cloud sequence (trial -> motion_cloud_id)
mc_seq_path = fullfile(cfg.motion_clouds_root, 'motion_cloud_sequence_250414.mat');
mc_sequence = [];
try
    if exist(mc_seq_path, 'file')
        S = load(mc_seq_path);
        % Primary source: table `presentation_sequence` (sequence is a single row across columns)
        if isfield(S, 'presentation_sequence')
            ps = S.presentation_sequence;
            if istable(ps)
                try
                    mc_sequence = ps{1, :};
                    fprintf('Loaded motion cloud ids from presentation_sequence (row-wise) in %s.\n', mc_seq_path);
                catch
                    warning('presentation_sequence found but could not extract row 1.');
                end
            elseif isnumeric(ps)
                mc_sequence = ps;
            end
        end
        % Fallbacks: legacy variable names or any numeric array
        if isempty(mc_sequence)
            if isfield(S, 'motion_cloud_sequence')
                mc_sequence = S.motion_cloud_sequence;
            elseif isfield(S, 'mc_sequence')
                mc_sequence = S.mc_sequence;
            else
                fn = fieldnames(S);
                for i = 1:numel(fn)
                    val = S.(fn{i});
                    if isnumeric(val) && (isvector(val) || ismatrix(val))
                        mc_sequence = val;
                        break;
                    end
                end
            end
            if ~isempty(mc_sequence)
                fprintf('Loaded motion cloud sequence from %s.\n', mc_seq_path);
            end
        end
    else
        warning('Motion cloud sequence file not found: %s. motion_cloud_id will be NaN.', mc_seq_path);
    end
catch ME
    warning(ME.identifier, '%s', sprintf('Failed to load motion cloud sequence: %s. motion_cloud_id will be NaN.', ME.message));
    mc_sequence = [];
end

% Load passive protocol sequence (trial_order -> protocol_id)
proto_seq_path = fullfile(cfg.motion_clouds_root, 'passive_protocol_sequence_motion_clouds.mat');
trial_order = [];
try
    if exist(proto_seq_path, 'file')
        P = load(proto_seq_path);
        if isfield(P, 'trial_order')
            trial_order = P.trial_order;
        else
            fns = fieldnames(P);
            for i = 1:numel(fns)
                val = P.(fns{i});
                if isnumeric(val) && (isvector(val) || ismatrix(val))
                    trial_order = val;
                    break;
                end
            end
        end
        fprintf('Loaded protocol sequence (trial_order) from %s.\n', proto_seq_path);
    else
        warning('Protocol sequence file not found: %s. protocol_id will be NaN.', proto_seq_path);
    end
catch ME
    warning(ME.identifier, '%s', sprintf('Failed to load protocol sequence: %s. protocol_id will be NaN.', ME.message));
    trial_order = [];
end

% No fixed delay used; using multi-lag indices only

for pid = 1:length(probe_ids)
    probe_id = probe_ids{pid};
    fprintf('\nProcessing probe %d/%d: %s\n', pid, length(probe_ids), probe_id);

    data              = ctl.load_formatted_data(probe_id);
    selected_clusters = data.selected_clusters();
    sessions          = data.motion_sessions();

    % (Removed unused id_to_region_map construction)

    % Outputs
    out_probe_id        = {};
    out_session_id      = {};
    out_trial_id        = [];
    out_protocol_id     = [];
    out_cluster_id      = [];
    out_spike_time      = [];
    out_position_cm     = [];
    out_motion_cloud_id = [];
    out_image_rel_path  = {};
    
    out_region_str      = {};
    % Multi-lag outputs (one column per lag) and lags listing
    out_image_index_lag40ms  = [];
    out_image_index_lag50ms  = [];
    out_image_index_lag60ms  = [];
    out_image_index_lag70ms  = [];
    out_image_index_lag80ms  = [];
    out_image_index_lag90ms  = [];
    out_image_index_lag100ms = [];
    out_image_index_lag110ms = [];
    out_image_index_lag120ms = [];

    kept = 0;

    for s = 1:length(sessions)
        session = sessions{s};
        trials_in_session = session.trials;

        for t = 1:length(trials_in_session)
            trial_aligned = trials_in_session{t}.to_aligned;

            t_full   = trial_aligned.probe_t();
            pos_full = trial_aligned.position();

            if isempty(t_full) || isempty(pos_full)
                continue;
            end

            motion_mask_full = trial_aligned.motion_mask();
            cc = [];
            if restrict_to_motion
                if ~any(motion_mask_full)
                    continue;
                end
                cc = bwconncomp(motion_mask_full);
                if cc.NumObjects == 0
                    continue;
                end
            end

            % motion_cloud_id and protocol_id for this trial
            tri = trial_aligned.trial_id;
            mc_id = nan; protocol_id = nan;
            if ~isempty(mc_sequence)
                seq = mc_sequence(:)';
                if tri >= 1 && tri <= numel(seq), mc_id = seq(tri); end
            end
            if ~isempty(trial_order)
                ord = trial_order(:)';
                if tri >= 1 && tri <= numel(ord), protocol_id = ord(tri); end
            end

            for c = 1:length(selected_clusters)
                cl = selected_clusters(c);
                st_all = cl.spike_times;
                if isempty(st_all), continue; end

                % Spikes inside trial span
                in_trial = st_all >= t_full(1) & st_all <= t_full(end);
                st_trial = st_all(in_trial);
                if isempty(st_trial), continue; end

                % Optional restriction to motion-only windows
                if restrict_to_motion
                    comp_starts = cellfun(@(idxs)(idxs(1)), cc.PixelIdxList);
                    comp_ends   = cellfun(@(idxs)(idxs(end)), cc.PixelIdxList);
                    t_starts    = t_full(comp_starts);
                    t_ends      = t_full(comp_ends);
                    in_motion = false(size(st_trial));
                    for k = 1:length(t_starts)
                        in_motion = in_motion | (st_trial >= t_starts(k) & st_trial <= t_ends(k));
                    end
                    if ~any(in_motion), continue; end
                    st_trial = st_trial(in_motion);
                end

                % Positions at spike times (current position)
                pos_at_spike = interp1(t_full, pos_full, st_trial, 'linear', 'extrap');

                % Use all spikes in trial; no single-delay indexing
                st_sel  = st_trial;
                pos_sel = pos_at_spike;
                if isempty(st_sel), continue; end

                % Multi-lag image indices (compute per-lag; NaN if out of bounds)
                lags_s_multi = lags_ms_multi / 1000;
                img_idx_lag40ms  = nan(size(st_sel));
                img_idx_lag50ms  = nan(size(st_sel));
                img_idx_lag60ms  = nan(size(st_sel));
                img_idx_lag70ms  = nan(size(st_sel));
                img_idx_lag80ms  = nan(size(st_sel));
                img_idx_lag90ms  = nan(size(st_sel));
                img_idx_lag100ms = nan(size(st_sel));
                img_idx_lag110ms = nan(size(st_sel));
                img_idx_lag120ms = nan(size(st_sel));

                for kk = 1:numel(lags_s_multi)
                    pt_k = st_sel - lags_s_multi(kk);
                    valid_k = pt_k >= t_full(1) & pt_k <= t_full(end);
                    if any(valid_k)
                        past_pos_k = interp1(t_full, pos_full, pt_k(valid_k), 'linear', 'extrap');
                        idx_k = floor(past_pos_k * 10);
                        if index_base == 0
                            idx_k = idx_k - 1;
                        end
                        switch lags_ms_multi(kk)
                            case 40
                                img_idx_lag40ms(valid_k)   = idx_k;
                            case 50
                                img_idx_lag50ms(valid_k)   = idx_k;
                            case 60
                                img_idx_lag60ms(valid_k)   = idx_k;
                            case 70
                                img_idx_lag70ms(valid_k)   = idx_k;
                            case 80
                                img_idx_lag80ms(valid_k)   = idx_k;
                            case 90
                                img_idx_lag90ms(valid_k)   = idx_k;
                            case 100
                                img_idx_lag100ms(valid_k)  = idx_k;
                            case 110
                                img_idx_lag110ms(valid_k)  = idx_k;
                            case 120
                                img_idx_lag120ms(valid_k)  = idx_k;
                        end
                    end
                end

                n_new = numel(st_sel);
                % Region string for this cluster id (if available)
                % Get region directly from cluster property (works for objects)
                reg_str = char(string(cl.region_str));
                % region_str is taken directly from cluster; no debug aborts

                out_probe_id        = [out_probe_id;        repmat({probe_id}, n_new, 1)]; %#ok<AGROW>
                out_session_id      = [out_session_id;      repmat({session.session_id}, n_new, 1)]; %#ok<AGROW>
                out_trial_id        = [out_trial_id;        repmat(tri, n_new, 1)]; %#ok<AGROW>
                out_protocol_id     = [out_protocol_id;     repmat(protocol_id, n_new, 1)]; %#ok<AGROW>
                out_cluster_id      = [out_cluster_id;      repmat(cl.id, n_new, 1)]; %#ok<AGROW>
                out_spike_time      = [out_spike_time;      st_sel(:)]; %#ok<AGROW>
                out_position_cm     = [out_position_cm;     pos_sel(:)]; %#ok<AGROW>
                out_motion_cloud_id = [out_motion_cloud_id; repmat(mc_id, n_new, 1)]; %#ok<AGROW>
                % image_rel_path: cloud directory name (relative)
                cloud_rel = '';
                if ~isnan(mc_id) && mc_id>=1 && mc_id<=numel(cloud_names)
                    candidate = cloud_names{mc_id};
                    % Ensure the directory exists under images_root
                    full_candidate = fullfile(images_root, candidate);
                    if exist(full_candidate, 'dir')
                        cloud_rel = candidate;
                    else
                        error('Expected folder for motion_cloud_id=%d not found: %s', mc_id, full_candidate);
                    end
                end
                out_image_rel_path  = [out_image_rel_path;  repmat({cloud_rel}, n_new, 1)]; %#ok<AGROW>
                
                out_region_str      = [out_region_str;      repmat({reg_str}, n_new, 1)]; %#ok<AGROW>
                % Multi-lag indices per-lag
                out_image_index_lag40ms  = [out_image_index_lag40ms;  img_idx_lag40ms(:)]; %#ok<AGROW>
                out_image_index_lag50ms  = [out_image_index_lag50ms;  img_idx_lag50ms(:)]; %#ok<AGROW>
                out_image_index_lag60ms  = [out_image_index_lag60ms;  img_idx_lag60ms(:)]; %#ok<AGROW>
                out_image_index_lag70ms  = [out_image_index_lag70ms;  img_idx_lag70ms(:)]; %#ok<AGROW>
                out_image_index_lag80ms  = [out_image_index_lag80ms;  img_idx_lag80ms(:)]; %#ok<AGROW>
                out_image_index_lag90ms  = [out_image_index_lag90ms;  img_idx_lag90ms(:)]; %#ok<AGROW>
                out_image_index_lag100ms = [out_image_index_lag100ms; img_idx_lag100ms(:)]; %#ok<AGROW>
                out_image_index_lag110ms = [out_image_index_lag110ms; img_idx_lag110ms(:)]; %#ok<AGROW>
                out_image_index_lag120ms = [out_image_index_lag120ms; img_idx_lag120ms(:)]; %#ok<AGROW>

                kept = kept + n_new;
            end
        end
    end

    fprintf('  Wrote %d spike rows for probe %s.\n', kept, probe_id);

    if save_csv && kept > 0
        % Reordered columns: motion_cloud_id immediately after session_id
        cm.create_table( ...
            'probe_id',           out_probe_id, ...
            'session_id',         out_session_id, ...
            'motion_cloud_id',    out_motion_cloud_id, ...
            'trial_id',           out_trial_id, ...
            'protocol_id',        out_protocol_id, ...
            'cluster_id',         out_cluster_id, ...
            'spike_time',         out_spike_time, ...
            'position_cm',        out_position_cm, ...
            'image_index_lag40ms',  out_image_index_lag40ms, ...
            'image_index_lag50ms',  out_image_index_lag50ms, ...
            'image_index_lag60ms',  out_image_index_lag60ms, ...
            'image_index_lag70ms',  out_image_index_lag70ms, ...
            'image_index_lag80ms',  out_image_index_lag80ms, ...
            'image_index_lag90ms',  out_image_index_lag90ms, ...
            'image_index_lag100ms', out_image_index_lag100ms, ...
            'image_index_lag110ms', out_image_index_lag110ms, ...
            'image_index_lag120ms', out_image_index_lag120ms, ...
            'image_rel_path',     out_image_rel_path, ...
            'region_str',         out_region_str);

        cm.save(probe_id);
        fprintf('  Saved CSV: %s\\%s.csv\n', cm.curr_dir, probe_id);
    else
        fprintf('  No rows to save or saving disabled.\n');
    end
end

fprintf('\nCompleted consolidated export for %d probe(s).\n', length(probe_ids));

function names_cell = extract_string_list_from_struct(S)
% Return a cell array of char row vector(s) from a struct containing
% either a cellstr, a cell of char, a string array, or a char vector.
% If nothing matches, return [].
    names_cell = [];
    try
        fns = fieldnames(S);
        for ii = 1:numel(fns)
            v = S.(fns{ii});
            if iscellstr(v) || (iscell(v) && all(cellfun(@ischar, v)))
                names_cell = v; break;
            elseif isstring(v)
                names_cell = cellstr(v(:)); break;
            elseif ischar(v)
                names_cell = cellstr(v); break;
            end
        end
        if ~isempty(names_cell)
            names_cell = names_cell(:)';
        end
    catch
        names_cell = [];
    end
end