function starts = pd_times(animal_id, session_obj)

pd = session_obj.photodiode;
fs = session_obj.fs;

% create a butterworth filter
fc = 8; % this is the cutoff frequency
[b, a] = butter(3, fc/(fs/2));
pd_filt = filtfilt(b, a, pd);

if strcmp(animal_id, 'CAA-1110262')
    th = 3.5;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;
elseif strcmp(animal_id, 'CAA-1110264')
    th = 2.5;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;
elseif strcmp(animal_id, 'CAA-1110265')
    th = 3.5;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;
elseif strcmp(animal_id, 'CAA-1112224')
    th = 1.75;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;
elseif strcmp(animal_id, 'CAA-1112416')
    th = 1.5;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;
elseif strcmp(animal_id, 'CAA-1112417') | strcmp(animal_id, 'CAA-1112531') | strcmp(animal_id, 'CAA-1112532')
    th = 1.5;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;
elseif strcmp(animal_id, 'CAA-1112529') | strcmp(animal_id, 'CAA-1112530')
    th = 2.5;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;
elseif strcmp(animal_id, 'CA_176_1') | strcmp(animal_id, 'CA_176_3')
    th = 3;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 6e4:6.2e4;
elseif strcmp(animal_id, 'CAA-1112872') | strcmp(animal_id, 'CAA-1112874')
    th = 2;  % CAA-1110262 = 3.5; CAA-1110264 = 2.5;
    start_window = 1.2e5:1.22e5;   
elseif strcmp(animal_id, 'CAA-1113220') | strcmp(animal_id, 'CAA-1113222')
    th = 2;  
    start_window = 1.2e5:1.22e5;  
elseif strcmp(animal_id, 'CAA-1114977') | strcmp(animal_id, 'CAA-1114978') | ...
        strcmp(animal_id, 'CAA-1114979') | strcmp(animal_id, 'CAA-1114980') | ...
        strcmp(animal_id, 'CAA-1115689') | strcmp(animal_id, 'CAA-1115691') 
    th = 10; 
    start_window = 1.2e5:1.22e5;    
end
% 
% th = 4.6; 
% start_window = 1.2e5:1.22e5;

idx = find(20*pd_filt(start_window) < th, 1);
starts = [];
starts(1) = start_window(1) + idx - 1;
%period = 0.25 * 10e4; 
period = 0.23 * 10e3;

direction = 1;
next_search = starts(1) + period;

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

if strcmp(animal_id, 'CA_176_1') | strcmp(animal_id, 'CA_176_3')
    assert(length(starts) == 5000);
else
    assert(length(starts) == 2500);
end

