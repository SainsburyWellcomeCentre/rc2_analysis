function data = load_data(formatted_data_fname)

% load the formatted data
data = load(formatted_data_fname);

% use the sync information to insert the time on the probe
for i = 1 : length(data.sessions)
    
    % HACK: this needs to be solved probably best by checking recording names
    % (rc_session must match camera session name)... need to store ID of
    % camera
    if length(data.cameras) < i || length(data.cameras(i).camera0) ~= length(data.t_trig{i})
        
        s_obj(i) = Session(data.sessions(i), data.t_sync{i});
    
    else
        
        s_obj(i) = Session(data.sessions(i), data.t_sync{i}, ...
            data.cameras(i), data.t_trig{i});
    end
end

data = rmfield(data, 'sessions');
data.sessions = s_obj;
