% Create temporal frequency tuning curve tables for motion cloud experiments
% These tables pool data across all three TF batches (different VR gains) and transform
% velocity to temporal frequency using batch-specific gain factors.
%
% Gain transformation: TF = velocity * gain, where:
%   Batch 1: gain = 1/30 Hz/(cm/s)  →  30 cm/s = 1 Hz
%   Batch 2: gain = 2/30 Hz/(cm/s)  →  30 cm/s = 2 Hz
%   Batch 3: gain = 4/30 Hz/(cm/s)  →  30 cm/s = 4 Hz
%
% Motion cloud batch patterns:
%   Batch 1: {'sf00p003_Bsf0p002_VX1p002', 'sf00p006_Bsf0p002_VX0p501', 'sf00p012_Bsf0p002_VX0p250'}
%   Batch 2: {'sf00p003_Bsf0p002_VX2p003', 'sf00p006_Bsf0p002_VX1p002', 'sf00p012_Bsf0p002_VX0p501'}
%   Batch 3: {'sf00p003_Bsf0p002_VX4p006', 'sf00p006_Bsf0p002_VX2p003', 'sf00p012_Bsf0p002_VX1p002'}
%
% This script creates TF tuning curves with 20 bins of equal data distribution,
% pooling velocity data from all three batches after TF transformation.
%
% Output: .mat files in formatted_data_dir/csvs/tf_tuning_curves/<probe_id>.mat
%
% See also: create_tuning_curves, create_tables_by_batch, TFTuningTable

%% Configuration
experiment_groups       = {'passive_same_luminance_mc'};
trial_group_labels      = {'VT', 'V', 'T_Vstatic'};

%% Initialize parallel pool with 20 workers
% Start parallel pool if not already running, ensure 20 workers
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool with 20 workers...\n');
    pool = parpool(20);
    fprintf('Parallel pool started with %d workers\n\n', pool.NumWorkers);
elseif pool.NumWorkers ~= 20
    fprintf('Restarting parallel pool to use 20 workers (currently has %d)...\n', pool.NumWorkers);
    delete(pool);
    pool = parpool(20);
    fprintf('Parallel pool restarted with %d workers\n\n', pool.NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n\n', pool.NumWorkers);
end

%% Initialize RC2Analysis controller
ctl = RC2Analysis();

%% Get all probe IDs from the specified experiment groups
probe_ids = ctl.get_probe_ids(experiment_groups{:});

fprintf('Creating TF tuning curves for %d probes\n', length(probe_ids));
fprintf('Trial groups: %s\n', strjoin(trial_group_labels, ', '));
fprintf('Pooling data across all 3 TF batches with gain transformation\n');
fprintf('Parallelizing shuffling computations (1000 shuffles per cluster)\n\n');

%% Main processing loop
for ii = 1 : length(probe_ids)
    
    fprintf('Processing probe %d/%d: %s\n', ii, length(probe_ids), probe_ids{ii});
    
    try
        % Create TF tuning curves for this probe
        % This will:
        %   1. Load all trials for each trial group
        %   2. Determine batch gain for each trial based on motion cloud
        %   3. Transform velocity to TF using batch-specific gains
        %   4. Create 20 TF bins with equal data distribution
        %   5. Compute firing rates in TF bins for each cluster
        %   6. Perform model selection (linear, Gaussian, asymmetric Gaussian, sigmoid)
        %   7. Perform shuffling for significance testing with 20 workers
        %   8. Save results to .mat file
        ctl.create_tf_tuning_curves(probe_ids{ii}, trial_group_labels, 'model_selection', true);
        
        fprintf('  Successfully created TF tuning curves\n');
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        fprintf('  Skipping this probe\n');
    end
    
    fprintf('\n');
end

fprintf('Completed: All TF tuning curves created\n');
fprintf('Files saved to: formatted_data_dir/csvs/tf_tuning_curves/\n');
