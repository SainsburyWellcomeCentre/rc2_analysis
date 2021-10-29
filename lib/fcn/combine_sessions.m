function sessions = combine_sessions(sessions, cameras, t_sync, t_trig)

for i = 1 : length(sessions)
    
    if ~isempty(cameras)
        sessions(i).camera0 = cameras(i).camera0;
        sessions(i).camera1 = cameras(i).camera1;
    end
    sessions(i).probe_t = t_sync{i};
    sessions(i).rc2_t = (0:sessions(i).n_samples-1)' * (1/sessions(i).fs);
    sessions(i).camera_t = t_trig;
end
