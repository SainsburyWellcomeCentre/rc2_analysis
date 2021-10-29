function [rc_probe_time_base, n_trig, t_trig] = synchronize(ks_dir, rc_fnames)

% we expect there to be at least 'expected_min_interval' between different
% AI recordings...
expected_min_interval = 2; % expect at least 2 seconds between sessions

% read params and trigger file
params              = read_ks_data(ks_dir, 'params');
trigger             = read_ks_data(ks_dir, 'trigger');

% get fs from params.py
fs                      = params.sample_rate;

% trigger times on the probe
probe_trig_times        = get_trigger_times(trigger) / fs;

% time between triggers
probe_trig_interval     = diff(probe_trig_times);

% find on which triggers the AI recordings start.
probe_trig_boundaries   = [1, find(probe_trig_interval > expected_min_interval) + 1];

% number of AI recordings
n_rc_recs = length(probe_trig_boundaries);

% make sure that the number of recordings we see on the probe, equals the
% number of files provided
assert(n_rc_recs == length(rc_fnames), ['Number of AI recordings provided doesn''t ', ...
        'match expected number from trigger channel'])

% for each expected AI recording
probe_n_triggers_per_rec = nan(n_rc_recs, 1);
for i = 1 : n_rc_recs-1
    probe_n_triggers_per_rec(i) = probe_trig_boundaries(i+1) - probe_trig_boundaries(i);
end
probe_n_triggers_per_rec(n_rc_recs) = length(probe_trig_times) - probe_trig_boundaries(n_rc_recs) + 1;

% preallocate
rc_probe_time_base = cell(1, n_rc_recs);
t_trig = cell(1, n_rc_recs);

% for each rc2 recording
for i = 1 : n_rc_recs
    
    % get times of triggers for this "session"
    t_trig{i} = probe_trig_times(probe_trig_boundaries(i) + (0:(probe_n_triggers_per_rec(i) - 1)));
    
    % read rc2 data
    [data, ~, ~, config] = read_rc2_bin(rc_fnames{i});
    
    % number of samples in rc2 recording
    n_samples = size(data, 1);
    
    % 
    rc_trig_interval = config.nidaq.co.low_samps + config.nidaq.co.high_samps;
    
    % the number of triggers we have sent out during the recording of the
    % analog input (a few are sent out after this :()
    n_trig_expected = ceil(n_samples / rc_trig_interval);
    
    
    % make sure the number expected is *roughly* the same as actual
%     assert(probe_n_triggers_per_rec(i) < n_trig_expected + 200 && ...
%         probe_n_triggers_per_rec(i) >= n_trig_expected, 'Number of triggers is not close to number of expected triggers');

    % get the start time of the AI acquisition on the probe
    probe_start_t = probe_trig_times(probe_trig_boundaries(i));
    
    if n_trig_expected > probe_n_triggers_per_rec(i)
        
        probe_end_t = probe_trig_times(probe_trig_boundaries(i) + probe_n_triggers_per_rec(i) - 1);
        
        % the number of samples acquired up to the rise of the very last
        % trigger sent during *acquisition*
        n_samples_accounted_for = (probe_n_triggers_per_rec(i) - 1) * rc_trig_interval;
        
        % this ratio is just slightly above 1
        f = n_samples/n_samples_accounted_for;
        
        % the duration of the AI acquisition is f*(end_t - start_t)
        rc_probe_time_base{i} =  linspace(probe_start_t, probe_start_t + f*(probe_end_t - probe_start_t), n_samples)';
        
    else
        % get the time of the last trigger sent *during acquisition* measured
        % on the probe
        probe_end_t = probe_trig_times(probe_trig_boundaries(i) + n_trig_expected - 1);
        
        % the number of samples acquired up to the rise of the very last
        % trigger sent during *acquisition*
        n_samples_accounted_for = (n_trig_expected - 1) * rc_trig_interval;
        
        % the number of samples after the rise of the very last trigger
        err = n_samples - n_samples_accounted_for;
        
        % make sure that this is less that the trigger interval
        assert(err >= 0 && err < rc_trig_interval, '??');
        
        % this ratio is just slightly above 1
        f = n_samples/n_samples_accounted_for;
        
        % the duration of the AI acquisition is f*(end_t - start_t)
        rc_probe_time_base{i} =  linspace(probe_start_t, probe_start_t + f*(probe_end_t - probe_start_t), n_samples)';
    end
    
    
    n_trig(i).rc = n_trig_expected;
    n_trig(i).probe = probe_n_triggers_per_rec(i);
end
