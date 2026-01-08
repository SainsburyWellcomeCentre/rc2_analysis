% Create tuning curve tables for velocity. 
% These tables are required by any script calling `load_tuning_curves`.

experiment_groups    = {'ambient_light'};
trial_group_labels   = {'RT'};

% Parallel processing configuration
n_workers = 4;  % Number of parallel workers (set to 0 to disable parallelization)

ctl                  = RC2Analysis();

% Enable overwrite mode when using parallel processing to avoid user input prompts
if n_workers > 0
    ctl.save.overwrite = true;
end

% Initialize parallel pool if parallelization is enabled
if n_workers > 0
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= n_workers
        delete(gcp('nocreate'));
        parpool(n_workers);
    end
end

for jj = 1 : length(trial_group_labels)
    probe_ids            = ctl.get_probe_ids(experiment_groups{jj});

    if n_workers > 0
        parfor ii = 1 : length(probe_ids)
            fprintf('Processing probe %d of %d\n', ii, length(probe_ids));
            ctl.create_tuning_curves(probe_ids{ii}, trial_group_labels);
        end
    else
        for ii = 1 : length(probe_ids)
            fprintf('Processing probe %d of %d\n', ii, length(probe_ids));
            ctl.create_tuning_curves(probe_ids{ii}, trial_group_labels);
        end
    end
end