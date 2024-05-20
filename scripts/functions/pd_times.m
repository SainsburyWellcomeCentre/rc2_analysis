function starts = pd_times(animal_id, session_obj, session_type)

pd = session_obj.photodiode;
fs = session_obj.fs;

[th, start_window, period, n_events, fc] = get_parameters_for_photodiode(animal_id, session_type);

% create a butterworth filter
% fc = 5; 
[b, a] = butter(3, fc/(fs/2));
pd_filt = filtfilt(b, a, pd);

idx = find(20*pd_filt(start_window) < th, 1);
starts = [];
starts(1) = start_window(1) + idx - 1;

direction = 1;
next_search = starts(1) + period;
% % 
% figure(1)
% hold on
% plot(pd)
% plot(20*pd_filt)
% yline(th)

while ~isempty(next_search)
%     scatter(starts(end), -0.5, 'black', '*')
    
    max_window = min(length(pd_filt)-next_search, 0.3*10e5);
    if direction == 1
        idx = find(20*pd_filt(next_search:next_search+max_window) > th, 1);
    else
        idx = find(20*pd_filt(next_search:next_search+max_window) < th, 1);
    end
    if isempty(idx)
        break
    end
    starts(end+1) = next_search + idx - 1;
    next_search = starts(end) + period;
    direction = mod(direction + 1, 2); 
end

starts(end) = [];

length(starts)

if strcmp(animal_id, 'CA_176_1') || strcmp(animal_id, 'CA_176_3')
    assert(length(starts) == n_events * 2);
else
    assert(length(starts) == n_events);
end

