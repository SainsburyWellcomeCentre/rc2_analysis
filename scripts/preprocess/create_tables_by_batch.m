% Create tuning curve tables for velocity filtered by temporal frequency batches.
% These tables are required by scripts that load tuning curves with TF batch filtering.
%
% This script creates separate tuning curve tables for three temporal frequency batches
% based on motion cloud patterns. Each batch has independent shuffle controls for 
% significance testing.
%
% Batch definitions (each matches ANY of the three patterns - OR logic):
%   Batch 1: {'sf00p003_Bsf0p002_VX1p002', 'sf00p006_Bsf0p002_VX0p501', 'sf00p012_Bsf0p002_VX0p250'}
%   Batch 2: {'sf00p003_Bsf0p002_VX2p003', 'sf00p006_Bsf0p002_VX1p002', 'sf00p012_Bsf0p002_VX0p501'}
%   Batch 3: {'sf00p003_Bsf0p002_VX4p006', 'sf00p006_Bsf0p002_VX2p003', 'sf00p012_Bsf0p002_VX1p002'}

%% Configuration
experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};

% Mode selection: 'all' creates single file with all trials (backward compatible)
%                 'separate' creates three separate files by TF batch
batch_mode              = 'separate';  % Options: 'all' or 'separate'

% Batch type: Determines which batch definitions to use and output folder
% Options: 'TF' (temporal frequency), 'SF' (spatial frequency), 'OR' (orientation)
batch_type              = 'TF';  % Change this to 'TF', 'SF', or 'OR'

%% Batches definitions - automatically selected based on batch_type

% Control stimuli to ALWAYS exclude (regardless of batch type)
exclude_patterns = {'theta0p000_Btheta3p142_sf00p006_Bsf0p004_VX0p000_BV2p000'};

if strcmp(batch_type, 'TF')
    % Temporal frequency batch patterns
    % Batch 1 patterns (low temporal frequencies)
    batch1_patterns = {'sf00p003_Bsf0p002_VX1p002', 'sf00p006_Bsf0p002_VX0p501', 'sf00p012_Bsf0p002_VX0p250'};
    
    % Batch 2 patterns (medium temporal frequencies)  
    batch2_patterns = {'sf00p003_Bsf0p002_VX2p003', 'sf00p006_Bsf0p002_VX1p002', 'sf00p012_Bsf0p002_VX0p501'};
    
    % Batch 3 patterns (high temporal frequencies)
    batch3_patterns = {'sf00p003_Bsf0p002_VX4p006', 'sf00p006_Bsf0p002_VX2p003', 'sf00p012_Bsf0p002_VX1p002'};
    
elseif strcmp(batch_type, 'SF')
    % Spatial frequency batch patterns
    % Batch 1 patterns
    batch1_patterns = {'sf00p003'};
    
    % Batch 2 patterns  
    batch2_patterns = {'sf00p006'};
    
    % Batch 3 patterns
    batch3_patterns = {'sf00p012'};
    
elseif strcmp(batch_type, 'OR')
    % Orientation batch patterns
    % Batch 1 patterns
    batch1_patterns = {'theta0p000'};
    
    % Batch 2 patterns  
    batch2_patterns = {'theta0p785'};
    
    % Batch 3 patterns 
    batch3_patterns = {'theta1p571'};
    
    % Batch 4 patterns 
    batch4_patterns = {'theta-0p785'};
    
else
    error('Invalid batch_type: must be ''TF'', ''SF'', or ''OR''');
end 

%% Helper function: Portable substring matcher (older MATLAB lacks contains/startsWith)
contains_any = @(str, subs) any(cellfun(@(s) ...
    (exist('contains','builtin') && contains(str, s)) || (~exist('contains','builtin') && ~isempty(strfind(str, s))), subs));

%% Load motion cloud sequence and folder names
mc_sequence = [];
cloud_names = {};

% Load presentation sequence
proto_seq_path = fullfile('D:\mvelez\mateoData_mc', 'motion_cloud_sequence_250414.mat');
if exist(proto_seq_path,'file')
    P = load(proto_seq_path);
    if isfield(P,'presentation_sequence')
        mc_sequence = P.presentation_sequence;
    else
        fns = fieldnames(P);
        for i3=1:numel(fns)
            v = P.(fns{i3});
            if isnumeric(v) && (isvector(v) || ismatrix(v))
                mc_sequence = v; break;
            end
        end
    end
else
    warning('Motion cloud sequence file not found: %s', proto_seq_path);
end

% Load cloud names from image_folders.mat
folders_path = fullfile('D:\mvelez\mateoData_mc', 'image_folders.mat');
if exist(folders_path,'file')
    S = load(folders_path);
    fns = fieldnames(S);
    for i3=1:numel(fns)
        v = S.(fns{i3});
        if iscell(v)
            cloud_names = v; break;
        elseif isstring(v)
            cloud_names = cellstr(v(:)); break;
        elseif ischar(v)
            cloud_names = cellstr(v); break;
        elseif isstruct(v)
            if isfield(v, 'name')
                try
                    cloud_names = {v.name}; break;
                catch
                end
            end
        end
    end
else
    warning('Image folders file not found: %s', folders_path);
end

if ~isempty(cloud_names)
    cloud_names = cloud_names(:)';
end

% Validate that we have the necessary data for filtering
if isempty(mc_sequence) || isempty(cloud_names)
    error('Cannot proceed: motion cloud sequence or cloud names not loaded. Check file paths.');
end

fprintf('Loaded motion cloud sequence with %d entries and %d cloud names.\n', ...
    length(mc_sequence), length(cloud_names));

%% Initialize RC2Analysis controller
ctl = RC2Analysis();

%% Main processing loop
if strcmp(batch_mode, 'all')
    % Backward compatible mode: create single file with all trials (original behavior)
    fprintf('\nMode: Creating single tuning curve file with all trials (backward compatible)\n');
    
    for jj = 1 : length(trial_group_labels)
        probe_ids = ctl.get_probe_ids(experiment_groups{jj});

        for ii = 1 : length(probe_ids)
            fprintf('Processing probe %d/%d: %s\n', ii, length(probe_ids), probe_ids{ii});
            ctl.create_tuning_curves(probe_ids{ii}, trial_group_labels);
        end
    end
    fprintf('\nCompleted: All tuning curves created in standard format.\n');
    
elseif strcmp(batch_mode, 'separate')
    % New mode: create separate files, one for each batch
    fprintf('\nMode: Creating separate tuning curve files for each %s batch\n', batch_type);
    
    probe_ids = ctl.get_probe_ids(experiment_groups{:});
    
    for ii = 1 : length(probe_ids)
        fprintf('\nProcessing probe %d/%d: %s\n', ii, length(probe_ids), probe_ids{ii});
        
        % Load formatted data for this probe
        data = ctl.load_formatted_data(probe_ids{ii});
        
        % Process each batch separately
        % Determine number of batches based on batch_type
        if strcmp(batch_type, 'OR')
            batch_names = {'batch1', 'batch2', 'batch3', 'batch4'};
            batch_patterns_list = {batch1_patterns, batch2_patterns, batch3_patterns, batch4_patterns};
        else
            batch_names = {'batch1', 'batch2', 'batch3'};
            batch_patterns_list = {batch1_patterns, batch2_patterns, batch3_patterns};
        end
        
        num_batches = length(batch_names);
        for batch_idx = 1:num_batches
            fprintf('  Creating tuning curves for %s...\n', batch_names{batch_idx});
            
            batch_patterns = batch_patterns_list{batch_idx};
            tuning_curves = cell(1, length(trial_group_labels));
            
            for jj = 1 : length(trial_group_labels)
                
                % Skip if the trial group label is not in the experiment
                if ~data.check_trial_group(trial_group_labels{jj})
                    fprintf('    Skipping %s (not found in experiment)\n', trial_group_labels{jj});
                    continue
                end
                
                % Get all trials for this trial group
                all_trials = data.get_trials_with_trial_group_label(trial_group_labels{jj});
                
                if isempty(all_trials)
                    fprintf('    No trials found for %s\n', trial_group_labels{jj});
                    continue
                end
                
                % Filter trials based on batch patterns
                filtered_trials = {};
                n_total = length(all_trials);
                
                for ti = 1:n_total
                    trial = all_trials{ti};
                    trial_id = trial.trial_id;
                    
                    % Check if this trial matches any pattern in current batch
                    if ~isempty(mc_sequence) && ~isempty(cloud_names) && ...
                       trial_id >= 1 && trial_id <= length(mc_sequence)
                        
                        mc_id = mc_sequence(trial_id);
                        
                        if mc_id >= 1 && mc_id <= length(cloud_names)
                            cname = cloud_names{mc_id};
                            
                            % First check if this is a control stimulus to exclude
                            if contains_any(cname, exclude_patterns)
                                continue;  % Skip this trial
                            end
                            
                            % Check if name contains ANY of the batch patterns (OR logic)
                            if contains_any(cname, batch_patterns)
                                filtered_trials{end+1} = trial;
                            end
                        end
                    end
                end
                
                fprintf('    %s: %d/%d trials matched %s patterns\n', ...
                    trial_group_labels{jj}, length(filtered_trials), n_total, batch_names{batch_idx});
                
                % Only process if we have filtered trials
                if ~isempty(filtered_trials)
                    % Create TuningTable object for this batch
                    tt = TuningTable(probe_ids{ii});
                    
                    % Align trials if they are replays
                    aligned_trials = cell(1, length(filtered_trials));
                    for kk = 1 : length(filtered_trials)
                        aligned_trials{kk} = filtered_trials{kk}.to_aligned;
                    end
                    
                    % Add filtered trials to the TuningTable object
                    % This computes velocity bins and prepares for shuffle controls
                    tt.add_trials(aligned_trials, trial_group_labels{jj});
                    
                    % Get all selected clusters and compute tuning curves
                    clusters = data.selected_clusters();
                    
                    for kk = 1 : length(clusters)
                        % tuning_curve method performs independent shuffling for this batch
                        tuning_curves{jj}(kk) = tt.tuning_curve(clusters(kk));
                    end
                    
                    fprintf('      Computed tuning curves for %d clusters\n', length(clusters));
                else
                    fprintf('      No trials to process after filtering\n');
                    tuning_curves{jj} = [];
                end
            end
            
            % Save tuning curves for this batch to a separate file
            tbl_struct.trial_groups = trial_group_labels;
            tbl_struct.tuning_curves = tuning_curves;
            tbl_struct.batch_name = batch_names{batch_idx};
            tbl_struct.batch_type = batch_type;
            tbl_struct.batch_patterns = batch_patterns;
            
            % Save with batch-specific filename including batch type
            save_tuning_curves_with_suffix(ctl, probe_ids{ii}, tbl_struct, batch_names{batch_idx}, batch_type);
            
            fprintf('  Saved tuning curves for %s\n', batch_names{batch_idx});
        end
        
        fprintf('Completed probe %s: all %d batches processed\n', probe_ids{ii}, num_batches);
    end
    
    fprintf('\nCompleted: All tuning curves created with batch separation.\n');
    if strcmp(batch_type, 'OR')
        fprintf('Files saved with suffixes: _batch1.mat, _batch2.mat, _batch3.mat, _batch4.mat\n');
    else
        fprintf('Files saved with suffixes: _batch1.mat, _batch2.mat, _batch3.mat\n');
    end
    
else
    error('Invalid batch_mode: must be ''all'' or ''separate''');
end

%% Helper function to save with custom suffix and batch type subdirectory
function save_tuning_curves_with_suffix(ctl, probe_id, tbl_struct, suffix, batch_type)
    % Get the standard filename and modify it to include batch suffix
    [fname, ~] = ctl.file.tuning_curves(probe_id);
    
    % Insert suffix before .mat extension and organize in batch_type subdirectory
    [filepath, name, ext] = fileparts(fname);
    
    % Create batch_type subdirectory (e.g., TF/, SF/, OR/)
    batch_type_dir = fullfile(filepath, batch_type);
    if ~exist(batch_type_dir, 'dir')
        mkdir(batch_type_dir);
    end
    
    fname_with_suffix = fullfile(batch_type_dir, sprintf('%s_%s%s', name, suffix, ext));
    
    % Save using the modified filename
    save(fname_with_suffix, '-struct', 'tbl_struct');
end
