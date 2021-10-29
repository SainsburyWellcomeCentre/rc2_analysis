classdef SpikeGLXRecording < SpikeGLXMetaData
    
    properties (SetAccess = private)
        
        data
        t
        
    end
    
    properties (SetAccess = private, Hidden = true)
        
        mmf_t
    end
    
    
    
    methods
        
        function obj = SpikeGLXRecording(glx_dir, type)
       %% SpikeGLXRecording()
        %   Class controlling access to raw AP or LFP data.
        %
        %       glx_dir = directory with SpikeGLX ap.bin and ap.meta data
        %       type = 'ap' or 'lf'
            
            obj = obj@SpikeGLXMetaData(glx_dir, type);
            
            % create a memory map to the bin file
            obj.mmf_t = memmapfile(obj.bin_fname, ...
                'format', {'int16', [obj.n_saved_channels, obj.n_samples], 'data'});
            
            % point directly to the data
            obj.data = obj.mmf_t.Data.data;
            
            % time base for the recording
            obj.t = (0:obj.n_samples-1)*(1/obj.fs);
        end
        
        
        
        function [val, t_out] = get_data_between_t(obj, channel_ids, t_in, type, legacy_on)
            
            % by default just return the int16 data
            VariableDefault('type', []);
            VariableDefault('legacy_on', false);
            
            % get index of channels in the data
            chan_idx = obj.get_idx_of_channel_id(channel_ids);
            
            % select the time points between the two values given
            if legacy_on
                % this option is left available to replicate values in old
                % work, but it's incorrectly indexed
                start_idx = obj.fs * t_in(1);
                n_samps = obj.fs * (t_in(2) - t_in(1));
                idx = start_idx + (0:n_samps-1);
            else
                idx = find(obj.t >= t_in(1) & obj.t < t_in(2));
            end
            
            
            % change the data output depending on type
            if strcmp(type, 's')
                % single value of data
                val = single(obj.data(chan_idx, idx));
            elseif strcmp(type, 'uv')
                % converted to uV
                val = obj.uV_per_bit * single(obj.data(chan_idx, idx));
            else
                % just the data
                val = obj.data(chan_idx, idx);
            end
            
            % timebase for outputted data
            t_out = obj.t(idx);
        end
        
        
        
        function val = trigger_data(obj)
        %%return int16 value of trigger data
            val = obj.data(obj.trigger_channel_idx, :);
        end
    end
end
