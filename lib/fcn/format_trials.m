function trials = format_trials(session)

% timebase of solenoid trace
% time_vector = (0:length(solenoid)-1)' * dt;
solenoid    = session.solenoid;
config      = session.config;
fs          = session.fs;

% make sure that the solenoid starts high in the session
assert(solenoid(1) > 2.5, 'Solenoid does not start high');

% get the sample points of solenoid up and down
solenoid_down = find(diff(solenoid > 2.5) == -1) + 1;

% HACK: should be solved at acquisition
if ~isfield(config, 'prot')
    solenoid_down = 1;
    config.prot.type = 'Visual';
end

% make sure that the number of trials from the solenoid is equal to the number of protocols
% in the config file
assert(length(solenoid_down) == length(config.prot) || ...
        length(solenoid_down) == length(config.prot) - 1, ...
    'Number of trials from solenoid does not equal number of trials saved in config file.');

% add an extra bound for the last trial
solenoid_down(end+1) = length(solenoid);


for i = 1 : length(solenoid_down)-1
    
    protocol = config.prot(i).type;
    
    switch protocol
        case {'Coupled', 'EncoderOnly', 'CoupledMismatch', 'EncoderOnlyMismatch'}
            chop_idx = (solenoid_down(i) - 4*fs) : solenoid_down(i+1);
        case {'StageOnly', 'ReplayOnly'}
            chop_idx = solenoid_down(i) : solenoid_down(i+1);
        case 'Visual'
            chop_idx = 1 : length(solenoid);
    end
    
    trials(i).id = i;
    trials(i).session_id = session.id;
    trials(i).fs = fs;
    trials(i).protocol = protocol;
    trials(i).start_idx = chop_idx(1);
    trials(i).end_idx = chop_idx(end);
    trials(i).config = config.prot(i);
%     trials(i).rc2_t = time_vector(chop_idx);
%     trials(i).probe_t = [];
end


%     % unpack the data
%     for chan_i = 1 : length(chan_names)
%         
%         % this channel name
%         str = chan_names{chan_i};
%         
%         % get chopped trace for this channel
%         trace = data_labelled.(str)(chop_idx);
%         
%         % store in trials
%         trials(i).(str) = trace;
%     end

