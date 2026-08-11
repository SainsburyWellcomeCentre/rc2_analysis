% Script to extract and plot speed/TF distributions across batches
% This analyzes the exact speed distributions and resulting TF distributions
% across the three temporal frequency batches.
%
% Key concept: TF = speed * batch_gain, where:
%   Batch 1: gain = 1/30 (30 cm/s -> 1 Hz)
%   Batch 2: gain = 2/30 (30 cm/s -> 2 Hz)
%   Batch 3: gain = 4/30 (30 cm/s -> 4 Hz)
%
% The speed distribution should be similar across batches (same animal behavior),
% but TF distributions will differ due to different VR gains.

%% Setup paths
% Get the directory of this script and add necessary paths
script_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(script_dir);  % Parent of scripts/
addpath(root_dir);
setup_paths;

fprintf('=== Speed and TF Distribution Analysis Across Batches ===\n\n');

%% Initialize
ctl = RC2Analysis();
experiment_groups = {'passive_same_luminance_mc'};

% Batch gain mappings (TF = velocity * gain)
batch_gains = [1/30, 2/30, 4/30];  % Hz/(cm/s) for batch 1, 2, 3
batch_labels = {'Batch 1 (1 Hz @ 30cm/s)', 'Batch 2 (2 Hz @ 30cm/s)', 'Batch 3 (4 Hz @ 30cm/s)'};

% Batch pattern definitions (from create_tables_by_batch.m)
batch_patterns = { ...
    {'sf00p003_Bsf0p002_VX1p002', 'sf00p006_Bsf0p002_VX0p501', 'sf00p012_Bsf0p002_VX0p250'}, ...
    {'sf00p003_Bsf0p002_VX2p003', 'sf00p006_Bsf0p002_VX1p002', 'sf00p012_Bsf0p002_VX0p501'}, ...
    {'sf00p003_Bsf0p002_VX4p006', 'sf00p006_Bsf0p002_VX2p003', 'sf00p012_Bsf0p002_VX1p002'} };

%% Load motion cloud sequence for trial->stimulus mapping
mc_sequence = [];
cloud_names = {};

% Load presentation sequence
proto_seq_path = fullfile(ctl.path_config.motion_clouds_root, 'motion_cloud_sequence_250414.mat');
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
folders_path = fullfile(ctl.path_config.motion_clouds_root, 'image_folders.mat');
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

fprintf('Loaded motion cloud sequence: %d entries, %d cloud names\n', length(mc_sequence), length(cloud_names));

%% Helper function to check if string contains any pattern
contains_any = @(str, subs) any(cellfun(@(s) contains(str, s), subs));

%% Get probe IDs
probe_ids = ctl.get_probe_ids(experiment_groups{:});
fprintf('Found %d probes\n\n', length(probe_ids));

%% Initialize storage for velocity values per batch
speed_per_batch = cell(1, 3);
tf_per_batch = cell(1, 3);

for batch_i = 1:3
    speed_per_batch{batch_i} = [];
    tf_per_batch{batch_i} = [];
end

%% Build trial_id -> batch mapping
trial_to_batch = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
for tid = 1:length(mc_sequence)
    mc_id = mc_sequence(tid);
    if mc_id < 1 || mc_id > length(cloud_names)
        continue;
    end
    cname = cloud_names{mc_id};
    
    for bi = 1:3
        if contains_any(cname, batch_patterns{bi})
            trial_to_batch(int32(tid)) = int32(bi);
            break;
        end
    end
end
fprintf('Mapped %d trials to batches\n\n', trial_to_batch.Count);

%% Collect speed/TF data from actual trials
fprintf('Collecting speed data from trials...\n');
n_trials_total = 0;

% Only need one probe since trials are the same across all probes
pid = probe_ids{1};
fprintf('  Using probe: %s\n', pid);

data = ctl.load_formatted_data(pid);

% Get VT trials (visual + treadmill)
trials = data.get_trials_with_trial_group_label('VT');
fprintf('    %d VT trials\n', length(trials));

for ti = 1:length(trials)
    trial = trials{ti};
    tid = trial.trial_id;
    
    % Check if this trial has a batch mapping
    if ~trial_to_batch.isKey(int32(tid))
        continue;
    end
    batch_idx = trial_to_batch(int32(tid));
    
    % Get velocity during motion periods
    vel = trial.velocity();
    mmask = trial.motion_mask();
    vel_motion = abs(vel(mmask));  % Speed (absolute velocity)
    
    % Store speed values
    speed_per_batch{batch_idx} = [speed_per_batch{batch_idx}; vel_motion];
    
    % Compute TF = speed * batch_gain
    tf_values = vel_motion * batch_gains(batch_idx);
    tf_per_batch{batch_idx} = [tf_per_batch{batch_idx}; tf_values];
    
    n_trials_total = n_trials_total + 1;
end

fprintf('\nTotal trials processed: %d\n\n', n_trials_total);

%% Report statistics
fprintf('=== Speed Distribution Statistics ===\n');
for bi = 1:3
    fprintf('%s:\n', batch_labels{bi});
    fprintf('  N samples: %d\n', length(speed_per_batch{bi}));
    if ~isempty(speed_per_batch{bi})
        fprintf('  Speed: mean=%.2f, std=%.2f, median=%.2f, [%.2f - %.2f] cm/s\n', ...
            mean(speed_per_batch{bi}), std(speed_per_batch{bi}), median(speed_per_batch{bi}), ...
            prctile(speed_per_batch{bi}, 5), prctile(speed_per_batch{bi}, 95));
        fprintf('  TF: mean=%.3f, std=%.3f, median=%.3f, [%.3f - %.3f] Hz\n\n', ...
            mean(tf_per_batch{bi}), std(tf_per_batch{bi}), median(tf_per_batch{bi}), ...
            prctile(tf_per_batch{bi}, 5), prctile(tf_per_batch{bi}, 95));
    else
        fprintf('  No data collected\n\n');
    end
end

%% Create figure with distributions
figure('Position', [100, 100, 1000, 400], 'Color', 'w');
colors = lines(3);

% Plot 1: Speed distributions (should be similar across batches)
subplot(1, 2, 1);
hold on;
speed_edges = linspace(0, 80, 50);
for bi = 1:3
    if ~isempty(speed_per_batch{bi})
        [counts, edges] = histcounts(speed_per_batch{bi}, speed_edges, 'Normalization', 'probability');
        centers = (edges(1:end-1) + edges(2:end))/2;
        plot(centers, counts * 100, 'LineWidth', 2, 'Color', colors(bi,:), 'DisplayName', batch_labels{bi});
    end
end
xlabel('Speed (cm/s)');
ylabel('Probability (%)');
title('Speed Distributions by Batch');
legend('Location', 'northeast', 'FontSize', 8);
grid on;

% Plot 2: TF distributions (should differ due to different gains)
subplot(1, 2, 2);
hold on;
tf_edges = linspace(0, 12, 50);
for bi = 1:3
    if ~isempty(tf_per_batch{bi})
        [counts, edges] = histcounts(tf_per_batch{bi}, tf_edges, 'Normalization', 'probability');
        centers = (edges(1:end-1) + edges(2:end))/2;
        plot(centers, counts * 100, 'LineWidth', 2, 'Color', colors(bi,:), 'DisplayName', batch_labels{bi});
    end
end
xlabel('Temporal Frequency (Hz)');
ylabel('Probability (%)');
title('TF Distributions by Batch (TF = Speed × Gain)');
legend('Location', 'northeast', 'FontSize', 8);
grid on;

%% Save figure
% Use RC2Analysis to get the figures directory
ctl.setup_figures({'tf_tuning_curves'}, true);
output_dir = ctl.figs.curr_dir;

% Save as PNG and PDF
saveas(gcf, fullfile(output_dir, 'tf_speed_distributions_by_batch.png'));

% Save as PDF with proper formatting
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print(gcf, fullfile(output_dir, 'tf_speed_distributions_by_batch.pdf'), '-dpdf', '-bestfit');

fprintf('\nFigures saved to: %s\n', output_dir);

%% Print percentile bin edges for reference
fprintf('\n=== Percentile Bin Edges (Pooled TF) ===\n');
all_tf = [];
for bi = 1:3
    all_tf = [all_tf; tf_per_batch{bi}];
end
if ~isempty(all_tf)
    prc_edges = prctile(all_tf, 0:5:100);
    fprintf('Bin edges (Hz): ');
    fprintf('%.3f ', prc_edges);
    fprintf('\n');
    
    fprintf('\nBin centers (Hz): ');
    bin_centers = (prc_edges(1:end-1) + prc_edges(2:end))/2;
    fprintf('%.3f ', bin_centers);
    fprintf('\n');
end

fprintf('\n=== Analysis Complete ===\n');
