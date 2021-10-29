function sessions = data_cleansing(sessions)


if strcmp(sessions(1).session_id, 'CA_176_1_rec1_001')
    % insert log names to make consistent with later data
    basename = 'C:\Users\Mateo\Desktop\mateoData\CA_176_1\CA_176_1\CA_176_1_rec1_001_single_trial';
elseif strcmp(sessions(1).session_id, 'CA_176_3_rec1_001')
    basename = 'C:\Users\Mateo\Desktop\mateoData\CA_176_3\CA_176_3\CA_176_3_rec1_001_single_trial';
end


% these data were run on an earlier version
% make consistent with later data.
if strcmp(sessions(1).session_id, 'CA_176_1_rec1_001') || ...
        strcmp(sessions(1).session_id, 'CA_176_3_rec1_001')
    
    c = 0;
    
    assert(length(sessions(1).trials) == length(sessions(1).config.prot));
    
    for t_i = 1 : length(sessions(1).trials)
        if strcmp(sessions(1).trials(t_i).config.type, 'Coupled') || ...
                strcmp(sessions(1).trials(t_i).config.type, 'EncoderOnly')
            c = c + 1;
            sessions(1).trials(t_i).config.log_fname = ...
                sprintf('%s_%03i.bin', basename, c);
            sessions(1).config.prot(t_i).log_fname = ...
                sprintf('%s_%03i.bin', basename, c);
        else
            sessions(1).trials(t_i).config.log_fname = '';
            sessions(1).config.prot(t_i).log_fname = '';
        end
    end
    
    % check that names match
    t = [sessions(1).trials(:).config];
    idx = find([t(:).follow_previous_protocol] == 1);
    
    assert(isequal({t(idx).wave_fname}, {t(idx-1).log_fname}));
    
    for t_i = 1 : length(sessions(1).trials)
        sessions(1).trials(t_i).config = rmfield(sessions(1).trials(t_i).config, 'follow_previous_protocol');
    end
    sessions(1).config.prot = rmfield(sessions(1).config.prot, 'follow_previous_protocol');
end


% some of these mice were restricted on the time they were run
if strcmp(sessions(1).session_id, 'CAA-1112414_rec1_001')
    sessions(1).trials(29:end) = [];
    sessions(1).config.prot(29:end) = [];
    sessions(1).n_trials = length(sessions(1).trials);
end
if strcmp(sessions(1).session_id, 'CAA-1110263_rec1_001')
    sessions(1).trials(1:36) = [];
    sessions(1).config.prot(1:36) = [];
    sessions(1).n_trials = length(sessions(1).trials);
end
if strcmp(sessions(1).session_id, 'CAA-1110265_rec1_001')
    sessions(1).trials(1:12) = [];
    sessions(1).config.prot(1:12) = [];
    sessions(1).n_trials = length(sessions(1).trials);
end


% remove last trial of these mice
if any(strcmp(sessions(1).session_id, {'CAA-1110262_rec1_001', 'CAA-1110263_rec1_001', ...
        'CAA-1110264_rec1_001', 'CAA-1110265_rec1_001', 'CAA-1112221_rec1_001', ...
        'CAA-1112222_rec1_001', 'CAA-1112223_rec1_001', ...
        'CAA-1112872_rec1_001'}))
    sessions(1).trials(end) = [];
    sessions(1).n_trials = length(sessions(1).trials);
    sessions(1).config.prot(end) = [];
end


% incorrect scaling of filtered_teensy_2
for sess_i = 1 : length(sessions)
    if any(strcmp(sessions(sess_i).session_id, {'CAA-1112872_rec1_001', 'CAA-1112872_rec1b_001', 'CAA-1112872_rec2_001', ...
            'CAA-1112874_rec1_001', 'CAA-1112874_rec2_001', ...
            'CAA-1113219_rec1_001', 'CAA-1113219_rec2_001', 'CAA-1113219_rec3_001', ...
            'CAA-1113220_rec1_001', 'CAA-1113220_rec2_001', 'CAA-1113220_rec3_001', ...
            'CAA-1113221_rec1_001', 'CAA-1113221_rec2_001', 'CAA-1113221_rec3_001', ...
            'CAA-1113222_rec1_001', 'CAA-1113222_rec2_001', 'CAA-1113222_rec3_001', ...
            'CAA-1114977_rec1_001', 'CAA-1114977_rec2_001', ...
            'CAA-1114978_rec1_001', 'CAA-1114978_rec2_001', ...
            'CAA-1114979_rec1_001', 'CAA-1114979_rec2_001', ...
            'CAA-1114980_rec1_001', 'CAA-1114980_rec2_001', ...
            'CAA-1115688_rec1_001', 'CAA-1115688_rec2_001', ...
            'CAA-1115689_rec1_001', 'CAA-1115689_rec2_001', ...
            'CAA-1115690_rec1_001', 'CAA-1115690_rec2_001', ...
            'CAA-1115691_rec1_001', 'CAA-1115691_rec2_001'}))
        
        t = sessions(sess_i).filtered_teensy_2;
        t = t / sessions(sess_i).config.nidaq.ai.scale(2);
        t = t + sessions(sess_i).config.nidaq.ai.offset(2);
        t = t - sessions(sess_i).config.nidaq.ai.offset(1);
        t = t * sessions(sess_i).config.nidaq.ai.scale(1);
        sessions(sess_i).filtered_teensy_2 = t;
    end
end



% correct trial ids for second session of CAA-1112872_rec1_rec1b_rec2_rec3
if strcmp(sessions(1).session_id, 'CAA-1112872_rec1_001')
    for ii = 1 : sessions(2).n_trials
        sessions(2).trials(ii).trial_id = sessions(2).trials(ii).trial_id + sessions(1).n_trials;
    end
end

