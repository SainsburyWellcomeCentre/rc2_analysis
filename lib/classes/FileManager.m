classdef FileManager < handle
    
    properties
        
        path_config
    end
    
    
    methods
        
        function obj = FileManager(path_config)
            
            obj.path_config = path_config;
        end
        
        
        
        function animal_list = list_animals(obj)
            
            contents = dir(obj.path_config.raw_probe_dir);
            idx = ~cellfun(@isempty, regexp({contents(:).name}, '^CA'));
            animal_list = {contents(idx).name}';
        end
        
        
        
        function [fname, exists] = experiment_list(obj)
            
            fname = fullfile(obj.path_config.summary_data_dir, 'experiment_list.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = session_list(obj)
            
            fname = fullfile(obj.path_config.summary_data_dir, 'session_list.csv');
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = processed_data_raw_ap_probe(obj, recording_id)
            
            animal_id = obj.animal_id_from_recording_id(recording_id);
            level_1 = fullfile(obj.path_config.processed_probe_slow_dir, animal_id);
            level_2 = fullfile(level_1, sprintf('%s_g0', recording_id));
            level_3 = fullfile(level_2, sprintf('%s_g0_imec0', recording_id));
            fname = fullfile(level_3, sprintf('%s_g0_t0.imec.ap.bin', recording_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = formatted_data(obj, recording_id)
            
            fname = fullfile(obj.path_config.formatted_data_dir, sprintf('%s.mat', recording_id));
            exists = isfile(fname);
        end
        
        
        
        function [fname, exists] = waveform_metrics_fix(obj, recording_id)
            
            dname = obj.imec0_ks2(recording_id);
            fname = fullfile(dname, 'csv', 'waveform_metrics_fix.csv');
            exists = isfile(fname);
        end
        
        
        
        function [dname, exists] = imec0_ks2(obj, recording_id)
            
            animal_id = obj.animal_id_from_recording_id(recording_id);
            level_1 = fullfile(obj.path_config.processed_probe_slow_dir, animal_id, 'output');
            level_2 = fullfile(level_1, sprintf('catgt_%s_g0', recording_id));
            level_3 = fullfile(level_2, sprintf('%s_g0_imec0', recording_id));
            dname = fullfile(level_3, 'imec0_ks2');
            exists = isfolder(dname);
        end
        
        
        
        function [fname, exists] = trigger_mat(obj, recording_id)
            
            dname = obj.imec0_ks2(recording_id);
            fname = fullfile(dname, 'trigger.mat');
            exists = isfile(fname);
        end
        
        
        
        function recording_ids = recording_id_from_animal_id(obj, animal_id)
            
            animal_dir = fullfile(obj.path_config.raw_probe_dir, animal_id);
            
            % get the contents of the directory
            contents = dir(animal_dir);
            
            % find the directories beginning with the animal name
            idx = find(~cellfun(@isempty, regexp({contents(:).name}, animal_id)));
            
            if ~isempty(regexp(contents(idx(1)).name, '.bin', 'once'))
                recording_ids = regexprep({contents(idx).name}, '_g0.+', '');
                recording_ids = unique(recording_ids);
            else
                recording_ids = {contents(idx).name};
            end
        end
    end
    
    
    
    methods (Static = true)
        
        function animal_id = animal_id_from_recording_id(recording_id)
            
            % for now assumes that rec1 is the non-animal string...
            idx = regexp(recording_id, '_rec1');
            
            animal_id = recording_id(1:idx(1)-1);
        end
    end
end
