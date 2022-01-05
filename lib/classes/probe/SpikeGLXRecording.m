classdef SpikeGLXRecording < SpikeGLXMetaData
% SpikeGLXRecording Class for reading and handling data in SpikeGLX
% .bin and .meta files
%
%   SpikeGLXRecording Properties:
%       data    - #samples x #channels memory map to data file
%       t       - #samples x 1, timebase of probe recording
%       mmf_t   - 
%
%   SpikeGLXRecording Methods:
%       get_data_between_t
%       trigger_data
%
%   This is a subclass of SpikeGLXMetaData.
%
%   See also: SpikeGLXMetaData


    properties (SetAccess = private)
        
        data
        t
        
    end
    
    properties (SetAccess = private, Hidden = true)
        
        mmf_t
    end
    
    
    
    methods
        
        function obj = SpikeGLXRecording(glx_dir, type)
        % SpikeGLXRecording
        %
        %   SpikeGLXRecording(DIRECTORY, RECORDING_TYPE)
        %
        %   Class controlling access to the SpikeGLX data in .bin
        %   contained in the directory DIRECTORY.
        %
        %   There can only be one *.ap.bin or *.lf.bin file in DIRECTORY
        %   The metadata file must also exist for this to work.
        %
        %   RECORDING_TYPE is either 'ap' or 'lf' and specifies which file
        %   to work with.  By default RECORDING_TYPE is 'ap'
            
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
        %%get_data_between_t Gets data between two timepoints
        %
        %   [DATA, T_OUT] = get_data_between_t(CHANNEL_IDS, T_IN, TYPE)
        %   gets the data between two time points.
        %
        %   Inputs:
        %       CHANNEL_IDS     - list of channel IDs for which to return data
        %       T_IN            - 1 x 2 array with time limits between which to return data
        %       TYPE            - 's', uv', or empty (default). If TYPE is
        %                         empty, the data is return as is without transformation. If
        %                         TYPE is 's', the data is converted to single precesion. If
        %                         TYPE is 'uv', the data returned is in microvolts.
        
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
        %%trigger_data Return the trigger channel
        %
        %   trigger_data() returns the int16 value of trigger channel
        
            val = obj.data(obj.trigger_channel_idx, :);
        end
    end
end
