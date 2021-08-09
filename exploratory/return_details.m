function return_details(data, recording_id, config)

this_data = get_data_for_recording_id(data, recording_id);

if ismember(recording_id, experiment_details('mismatch_nov20'))
    exp_obj = MismatchExperiment(this_data, config);
elseif ismember(recording_id, experiment_details('visual_flow'))
    exp_obj = VisualFlowExperiment(this_data, config);
end

fprintf('\n\n\n%s\n', recording_id);
fprintf('# trials: %i\n', length(exp_obj.trials));


for i = 1 : length(exp_obj.protocol_ids)
    
    these_trials = exp_obj.trials_of_type(exp_obj.protocol_ids(i));
    fprintf('# trials, %s: %i\n', exp_obj.protocol_label{i}, length(these_trials));
    
    if ismember(recording_id, experiment_details('mismatch_nov20'))
        
        n_failed = 0;
        for j = 1 : length(these_trials)
            
            mm_onset_t = these_trials(j).mismatch_onset_t();
            mm_offset_t = these_trials(j).mismatch_offset_t();
            
            if mm_offset_t - mm_onset_t < 0.05
                n_failed = n_failed + 1;
            end
        end
        
        fprintf(' # failed: %i\n', n_failed);
    end
end
