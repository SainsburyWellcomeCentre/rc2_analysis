classdef SpikeGLXMetaData < handle
    
    properties (SetAccess = private)
        
        type
        
        bin_fname
        meta_fname
        
        shank_ids_used
        n_shanks_used
        
        n_saved_channels
        n_samples
        fs
        n_bit
        uV_per_bit
        
        all_channel_ids
        
        trigger_channel_idx
        reference_channel_idx
        
        reference_channel_ids
        
        config
        
        probe
    end
    
    properties (SetAccess = private, Hidden = true)
        
        ai_volts
    end
    
    properties (Dependent = true)
        
        gain
        n_electrode_columns_per_shank
    end
    
    
    
    methods
        
        function obj = SpikeGLXMetaData(glx_dir, type)
       %% SpikeGLXMetaData()
        %   Class containing config information about the SpikeGLX recording
        %   contained in the directory 'glx_dir'.
        %   There can only be one *.ap.bin or *.lf.bin file in glx_dir
        %   Metadata file must also exist for this to work
        %
        %       glx_dir = directory with SpikeGLX ap.bin and ap.meta data
        %       type = 'ap' or 'lf'
        
            % 'ap' or 'lf'
            obj.type = type;
            
            % find the *.bin data
            dir_contents = dir(glx_dir);
            
            % find entry with the right *.bin file
            file_idx = contains({dir_contents(:).name}, sprintf('%s.bin', type));
            
            % make sure there is one and only one file
            assert(sum(file_idx) ~= 0, 'no .bin file found in %s.', glx_dir);
            assert(sum(file_idx) == 1, 'too many .bin files found in %s.', glx_dir);
            
            % store the bin fname
            obj.bin_fname = fullfile(dir_contents(file_idx).folder, ...
                                     dir_contents(file_idx).name);
            
            % path to associated meta file
            obj.meta_fname = strrep(obj.bin_fname, ...
                                    sprintf('.%s.bin', obj.type), ...
                                    sprintf('.%s.meta', obj.type));
            
            % make sure it exists
            assert(isfile(obj.meta_fname), 'associated .meta file doesn''t exist in %s', glx_dir)
            
            % load the meta data
            obj.config              = obj.read_spikeglx_config(obj.meta_fname);
            
            % set probe type
            obj.probe               = Probe(obj.get_probe_type());
            
            % store some variables for ease
            obj.n_saved_channels    = str2double(obj.config.nSavedChans);
            obj.fs                  = str2double(obj.config.imSampRate);
            obj.n_samples           = dir_contents(file_idx).bytes / (2 * obj.n_saved_channels);
            
            obj.uV_per_bit          = 1e6 * (str2double(obj.config.imAiRangeMax)-str2double(obj.config.imAiRangeMin)) / (2^obj.n_bit) / obj.gain;
        end
        
        
        
        function val = get.n_bit(obj)
            
            if strcmp(obj.get_probe_type(), '3A')
                val = 10;
            elseif strcmp(obj.get_probe_type(), '24')
                val = 14;
            end
        end
        
        
        
        function val = get.gain(obj)
        %%gain, dependent on whether 'ap' or 'lf'
            if strcmp(obj.get_probe_type(), '3A')
                if strcmp(obj.type, 'ap')
                    val = 500;
                elseif strcmp(obj.type, 'lf')
                    val = 250;
                end
            elseif strcmp(obj.get_probe_type(), '24')
                val = 80;
            end
        end
        
        
        
        function val = get.trigger_channel_idx(obj)
        %%returns the index in the saved data corresponding to the trigger
        %%channel
            entries     = textscan(obj.config.snsChanMap, '(%2s%d;%d:%d', 'EndOfLine', ')', 'HeaderLines', 1);
            val         = find(strcmp(entries{1}, 'SY'));
        end
        
        
        
        function val = get.all_channel_ids(obj)
        %%for each saved channel return the channel ID
        %   trigger channel comes back as -1
            entries     = textscan(obj.config.snsChanMap, '(%2s%d;%d:%d', 'EndOfLine', ')', 'HeaderLines', 1);
            val = entries{2};
            idx = strcmp(entries{1}, 'SY');
            val(idx) = -1;
        end
        
        
        
        function val = get.reference_channel_idx(obj)
        %%get the index of the reference channels in the data
            electrode_ids = obj.electrode_id_from_channel_id(obj.all_channel_ids);
            val = find(ismember(electrode_ids, obj.probe.reference_electrode_ids));
        end
        
        
        
        function val = get.reference_channel_ids(obj)
        %%get the channel id on which the reference channels occur
            index = obj.reference_channel_idx;
            val = obj.all_channel_ids(index);
        end
        
        
        
        function val = get.n_shanks_used(obj)
        %%get the number of shanks used
            val = length(obj.shank_ids_used);
        end
        
        
        
        function val = get.shank_ids_used(obj)
        %%get which shank IDs were used (zero indexed)
            [~, val] = obj.electrode_id_from_channel_id(obj.all_channel_ids);
            val(isnan(val)) = [];
            val = unique(val);
        end
        
        
        
        function val = get.n_electrode_columns_per_shank(obj)
        %%return number of columns of electrodes per shank
        %   staggering of electrode sites is accounted for
            val = obj.probe.n_electrode_columns_per_shank();
        end
        
        
        
        function val = channel_from_tip_um(obj, channel_id)
        %%return microns from tip of channel with ID 'channel_id' (zero
        %   indexed)
            
            assert(all(ismember(channel_id, obj.all_channel_ids)) & all(channel_id >= 0), 'not all channel ids exist')
            electrode_id = obj.electrode_id_from_channel_id(channel_id);
            val = obj.probe.electrode_from_tip_um(electrode_id);
        end
        
        
        
        function val = channel_from_left_um(obj, channel_ids)
        %%return microns from leftmost electrode of channel with ID 'channel_id' (zero
        %   indexed)
        
            assert(all(ismember(channel_ids, obj.all_channel_ids)) & all(channel_ids >= 0), 'not all channel ids exist')
            [electrode_ids, shank_ids] = obj.electrode_id_from_channel_id(channel_ids);
            val = obj.probe.electrode_from_left_um(electrode_ids, shank_ids);
        end
        
        
        
        function type = get_probe_type(obj)
        %%get type of probe from metadata file
        %   e.g. 3A for Phase 3A,  24 for Neuropixels 2.0, 4 shank
        
            if isfield(obj.config, 'imDatPrb_type')
                type = obj.config.imDatPrb_type;
            else
                type = '3A'; % 3A probe
            end
        end
        
        
        
        function [electrode_id, shank_id] = electrode_id_from_channel_id(obj, channel_ids)
        %%for a channel ID (or vector of channel IDs), get the shank and
        %   electrode which it occurred
        
            if strcmp(obj.get_probe_type, '3A')
                entries         = textscan(obj.config.imroTbl, '(%d %d %*s %*s %*s', 'EndOfLine', ')', 'HeaderLines', 1);
                channel_map_id  = entries{1};
                bank_id         = entries{2};
                shank_map_id        = zeros(length(entries{1}), 1);
                electrode_map_id    = 384*bank_id + channel_map_id;
            elseif strcmp(obj.get_probe_type, '24')
                entries         = textscan(obj.config.imroTbl, '(%d %d %*s %*s %d', 'EndOfLine', ')', 'HeaderLines', 1);
                channel_map_id  = entries{1};
                shank_map_id        = entries{2};
                electrode_map_id    = entries{3};
            end
            
            electrode_id    = nan(length(channel_ids), 1);
            shank_id        = nan(length(channel_ids), 1);
            
            for ii = 1 : length(channel_ids)
                idx = find(channel_map_id == channel_ids(ii));
                if ~isempty(idx)
                    electrode_id(ii) = electrode_map_id(idx);
                    shank_id(ii) = shank_map_id(idx);
                end
            end
        end
        
        
        
        function idx = get_idx_of_channel_id(obj, channel_ids)
        %%returns the index in the data of the channel IDs
        %   channel ID must be >= 0 and exist in the data or error is
        %   thrown
        
            assert(all(ismember(channel_ids, obj.all_channel_ids)) & all(channel_ids >= 0), 'not all channel ids exist')
            idx = arrayfun(@(x)(find(obj.all_channel_ids == x)), channel_ids);
        end
        
        
        
        function channel_ids = get_channel_ids_along_shank(obj, shank_id)
        %%given a shank ID, returns a list of channel IDs working from left
        %   to right, bottom to top
            
            % remove trigger channel
            channel_ids = obj.all_channel_ids;
            channel_ids(channel_ids == -1) = [];
            
            % get shank of each channel
            [~, shank_ids] = obj.electrode_id_from_channel_id(channel_ids);
            
            % restrict channels on this shank
            idx = shank_ids == shank_id;
            channel_ids = channel_ids(idx);
            
            % sort the channels from left
            from_left = obj.channel_from_left_um(channel_ids);
            [~, from_left_sort] = sort(from_left, 'ascend');
            channel_ids = channel_ids(from_left_sort);
            
            % sort the channels from tip
            from_tip = obj.channel_from_tip_um(channel_ids);
            [~, from_tip_sort] = sort(from_tip, 'ascend');
            channel_ids = channel_ids(from_tip_sort);
        end
        
        
        
        function channel_ids = get_channel_ids_along_column(obj, shank_id, column_id)
        %%given a shank ID and column ID (zero indexed) on that shank, 
        %   returns a list of channel ids from bottom to top of the probe
            
            % channels on the shank
            channel_ids = obj.get_channel_ids_along_shank(shank_id);
            
            % distance of channel from leftmost column
            channel_from_left_um = obj.channel_from_left_um(channel_ids);
            
            % columns from left
            column_from_left_um = unique(channel_from_left_um, 'stable');
            
            % get the index of the channels in specified column
            idx = ismember(channel_from_left_um, column_from_left_um(column_id+1));
            channel_ids = channel_ids(idx);
        end
        
        
        
        function channel_ids = order_channel_ids(obj)
        %%get a "natural" ordering of all the channels on the probe from 
        %   left to right, bottom to top, across all shanks
        
            n_shanks = obj.probe.n_shanks;
            channel_ids = [];
            for ii = 1 : n_shanks
                channel_ids = [channel_ids; obj.get_channel_ids_along_shank(ii-1)];
            end
        end
    end
    
    
    
    methods (Static = true)
        
        function config = read_spikeglx_config(bin_fname)
        % get the configuration data (.meta file) for an associated .bin file
        % you supply the .bin filename
            
            % generated the .meta filename from supplied bin
            meta_fname = strrep(bin_fname, '.bin', '.meta');
            
            % read in the config file as a single string
            fid = fopen(meta_fname, 'r');
            str = fread(fid, 'char');
            str = char(str(:))';
            fclose(fid);
            
            % split the config into list
            str_lines = strsplit(str, '\r\n');
            
            % remove last entry
            str_lines(end) = [];
            
            for i = 1 : length(str_lines)
                
                % split on '='
                keyval = strsplit(str_lines{i}, '=');
                
                % some start with ~
                if keyval{1}(1) == '~'
                    keyval{1}(1) = [];
                end
                
                % store in config structure
                config.(keyval{1}) = keyval{2};
            end
        end
    end
end
