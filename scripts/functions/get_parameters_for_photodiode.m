function [th, start_window, period, n_events, fc] = get_parameters_for_photodiode(animal_id, session_type)
%%GET_PARAMETERS_FOR_PHOTODIODE Get parameters to extract starting times from photodiode signal
%
%   [TH, START_WINDOW, PERIOD, N_EVENTS, FC] = get_parameters_for_photodiode(ANIMAL_ID, SESSION_TYPE)
%   A function to get pre-selected parameters useful to find strarting times 
%   of the photodiode signal. The optimal parameters depend on the experiment, 
%   which is identified with ANIMAL_ID. 
%   There are two experimental protocols that use the photodiode, sparse_noise 
%   and sf_tf, specified in the SESSION_TYPE argument. The first is used to 
%   measure receptive fields, the second to quantify drifting gratings preference.
%
%   The optimal parameters change depending on which protocol is run.
%   This function returns the following objects:
%   - TH: threshold to be applied to the band-passed signal
%   - START_WINDOW: initial part of the experiment before the sequence of stimuli starts
%   - PERIOD: minimum duration of a stimulus
%   - N_EVENTS: number of stimuli in session
%   - FC: cutoff frequency


if session_type == "sparse_noise"
    start_window = 1.2e5:1.22e5; 
    if strcmp(animal_id, 'CAA-1110262')
        th = 3.5;  
    elseif strcmp(animal_id, 'CAA-1110264')
        th = 2.5; 
    elseif strcmp(animal_id, 'CAA-1110265')
        th = 3.5;  
    elseif strcmp(animal_id, 'CAA-1112224')
        th = 1.75;  
    elseif strcmp(animal_id, 'CAA-1112416')
        th = 1.5;  
    elseif strcmp(animal_id, 'CAA-1112417') || strcmp(animal_id, 'CAA-1112531') || strcmp(animal_id, 'CAA-1112532')
        th = 1.5; 
    elseif strcmp(animal_id, 'CAA-1112529') || strcmp(animal_id, 'CAA-1112530')
        th = 2.5; 
    elseif strcmp(animal_id, 'CA_176_1') || strcmp(animal_id, 'CA_176_3')
        th = 3; 
        start_window = 6e4:6.2e4;
    elseif strcmp(animal_id, 'CAA-1112872') || strcmp(animal_id, 'CAA-1112874')
        th = 2;  
    elseif strcmp(animal_id, 'CAA-1113220') || strcmp(animal_id, 'CAA-1113222')
        th = 2;  
    elseif strcmp(animal_id, 'CAA-1114977') || strcmp(animal_id, 'CAA-1114978') || ...
            strcmp(animal_id, 'CAA-1114979') || strcmp(animal_id, 'CAA-1114980') || ...
            strcmp(animal_id, 'CAA-1115689') || strcmp(animal_id, 'CAA-1115691') 
        th = 10; 
    end
    period = 0.23 * 10e3;
    n_events = 2500;
    fc = 8;
     
elseif session_type == "sf_tf"
    if strcmp(animal_id, 'CAA-1110264')
        th = 2.2; 
        reps = 3;
    elseif strcmp(animal_id, 'CAA-1110262') || strcmp(animal_id, 'CAA-1110265')
        th = 2.5;
        reps = 5;
    else
        th = 1.5;
        reps = 5;
    end
    start_window = 1.2e5:1.22e5;
    period = 0.1 * 10e4;
    n_events = 128 * reps;
    fc = 5;
end