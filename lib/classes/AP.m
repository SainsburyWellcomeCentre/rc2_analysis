classdef AP < handle
    
    properties (SetAccess = private)
        
        bin_fname
        meta_fname
        
        data
        t
        
        n_channels
        n_samples
        fs
        uV_per_bit
        
        trigger_channel
        reference_channels = [37, 76, 113, 152, 189, 228, 265, 304, 341, 380];
        channels
        
        config
        
        probe = Probe()
    end
    
    properties (SetAccess = private, Hidden = true)
        
        mmf_t
        ai_volts
    end

    methods
        
        function obj = AP(glx_dir)
        %% AP()
        %   Class controlling access to raw AP data and computation of
        %   properties.
        %       glx_dir = directory with SpikeGLX ap.bin and ap.meta data
        
            % find the lf.bin data
            dir_contents = dir(glx_dir);
            
            % find entry with the lf.bin file
            lf_idx = contains({dir_contents(:).name}, 'ap.bin');
            
            % make sure there is one and only one file
            assert(sum(lf_idx) == 1, 'more or less than 1 .ap.bin file found.');
            
            % store the bin fname
            obj.bin_fname = fullfile(dir_contents(lf_idx).folder, ...
                                     dir_contents(lf_idx).name);
            
            % path to associated meta file
            obj.meta_fname = strrep(obj.bin_fname, '.ap.bin', '.ap.meta');
            
            % make sure it exists
            assert(isfile(obj.meta_fname), 'associated ap.meta file doesn''t exist')
            
            % load the meta data
            obj.config              = read_spikeglx_config(obj.meta_fname);
            
            % store some variables for ease
            obj.n_channels          = str2double(obj.config.nSavedChans);
            obj.fs                  = str2double(obj.config.imSampRate);
            obj.n_samples           = dir_contents(lf_idx).bytes / (2 * obj.n_channels);
            obj.uV_per_bit          = 1e6 * str2double(obj.config.imAiRangeMax) / 512 / 500;
            
            % remove reference channels above the number of recorded
            % channels
            obj.reference_channels(obj.reference_channels > obj.n_channels) = [];
            I = regexp(obj.config.snsChanMap, '(');
            obj.trigger_channel     = find(arrayfun(@(x)(strcmp(obj.config.snsChanMap(x+(0:3)), '(SY0')), I))-1;
            
            % genuine channels
            obj.channels = 1 : obj.n_channels;
            obj.channels([obj.reference_channels, obj.trigger_channel]) = [];
            
            % create a memory map to the bin file
            obj.mmf_t = memmapfile(obj.bin_fname, ...
                'format', {'int16', [obj.n_channels, obj.n_samples], 'data'});
            
            % point directly to the data
            obj.data = obj.mmf_t.Data.data;
            
            % time base for the LFP recording
            obj.t = (0:obj.n_samples-1)*(1/obj.fs);
        end
        
        
        
        function [val, t_out] = get_data_between_t(obj, chans, t_in, type)
            
            % by default just return the int16 data
            VariableDefault('type', []);
            
            % select the time points between the two values given
            idx = obj.t > t_in(1) & obj.t < t_in(2);
            
            % change the data output depending on type
            if strcmp(type, 's')
                % single value of data
                val = single(obj.data(chans, idx));
            elseif strcmp(type, 'uv')
                % converted to uV
                val = obj.uV_per_bit * single(obj.data(chans, idx));
            else
                % just the data
                val = obj.data(chans, idx);
            end
            
            % timebase for outputted data
            t_out = obj.t(idx);
        end
        
        
        
        function val = sdata(obj, chans, samps)
            VariableDefault('chans', ':');
            VariableDefault('samps', ':');
            val = single(obj.data(chans, samps));
        end
        
        
        
        function val = vdata(obj, chans, samps)
            VariableDefault('chans', ':');
            VariableDefault('samps', ':');
            val = obj.uV_per_bit * single(obj.data(chans, samps));
        end
        
        
        function val = bank(obj, n)
            val = intersect(obj.channels, n:4:obj.n_channels);
        end
        
        
%         function get_power(
%             val = bandpower(single([mmf.Data.x(i, :)'; mmf.Data.x(i, end)]), fs, [500, min(fs/2, 2500)]);
%         end
    end
end
