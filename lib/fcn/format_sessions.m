function sessions = format_sessions(rc_names)

% sessions = cell(1, length(rc_names));

for i = 1 : length(rc_names)
    
    % read the bin and cfg file
    [data, dt, chan_names, config] = read_rc2_bin(rc_names{i});
    
    % get filename without full path or extension
    [~, sessions(i).id] = fileparts(rc_names{i});
    
    sessions(i).fs = 1/dt;
    sessions(i).config = config;
    
    sessions(i).n_samples = size(data, 1);
    sessions(i).n_channels = size(data, 2);
    
    % unpack data matrix
    for j = 1 : length(chan_names)
        sessions(i).(chan_names{j}) = data(:, j);
    end
    
    sessions(i).trials = format_trials(sessions(i));
    sessions(i).n_trials = length(sessions(i).trials);
end
