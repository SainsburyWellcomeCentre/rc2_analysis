function starts = photodiode_times(animal_id, session, session_type)
%%PHOTODIODE_TIMES Get starting stimulus times from photodiode signal
%
%   STARTS = photodiode_times(ANIMAL_ID, SESSION, SESSION_TYPE)
%   Uses a Butterworth filter to flatten the signal of the photodiode
%   and simplify the identification of starting times. The signal of the
%   photodiode comes from a flickering monitor and the reading is noisy.
%   Any change from white to black and vice versa marks the starting of
%   a new stimulus.
%   From get_parameters_for_photodiode receives the cutoff frequency (fc)
%   and the pre-selected threshold depending on ANIMAL_ID and SESSION_TYPE.
%   The SESSION object contains photodiode signal and sampling rate information.
%   SESSION is a RC2Session objcet.
%   Returns the starting times in STARTS.
%   
%   See also: get_parameters_for_photodiode

pd = session.photodiode;
fs = session.fs;

% get pre-selected parameters depending on animal id and condition
[th, start_window, period, n_events, fc] = get_parameters_for_photodiode(animal_id, session_type);

% create a butterworth filter
[b, a]  = butter(3, fc/(fs/2));
pd_filt = filtfilt(b, a, pd);

% find first starting point
idx       = find(20*pd_filt(start_window) < th, 1);
starts    = [];
starts(1) = start_window(1) + idx - 1;

% define next search
direction   = 1;
next_search = starts(1) + period;

% search for startinf times
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

% check number of starting points found
if strcmp(animal_id, 'CA_176_1') || strcmp(animal_id, 'CA_176_3')
    assert(length(starts) == n_events * 2);
else
    assert(length(starts) == n_events);
end

