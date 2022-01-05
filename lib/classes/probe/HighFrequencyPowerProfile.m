classdef HighFrequencyPowerProfile < handle
% HighFrequencyPowerProfile Class for analyzing the high-frequency power
% profile across the probe
%
%   HighFrequencyPowerProfile Properties:
%       first_batch_s           - start time of first batch (default = 60s)
%       batch_interval_s        - interval between batches (default = 300s)
%       batch_duration_s        - duration of each batch (default = 10s)
%       n_batches_to_process    - number of batches to process (default = inf i.e. all data)
%       lower_hz                - lower frequency of bandpower to calculate (default = 500Hz)
%       upper_hz                - upper frequency of bandpower to calculate (default = 5000Hz)
%       smoothing_distance_um   - size of the sliding window filter (default = 100um)
%       interp_distance_um      - resolution of interpolation of power profile (default = 1um)
%       channel_ids_to_remove   - IDs of the channels to remove from processing
%       search_above            - threshold in um to search above for the peak power (default = 1000)
%       search_below            - threshold in um to search below for the peak power (default = inf)
%       batches_to_use          - which batches to use (default 1 to 10)
%       power_raw               - raw power on each channel for each batch
%       probe_track             - instance of ProbeTrack
%       clusters_from_tip_um    - all MUA units and their distance from probe tip
%       ephys_l5                - location of the electrophysiological L5
%       anatomy_l5              - location of the anatomical L5 (from ProbeTrack)
%       delta_l5                - difference between ephys_l5 and anatomy_l5
%       probe_id                - string with probe recording ID
%
%   HighFrequencyPowerProfile Methods:
%       channel_ids_to_process              - ID of the channels to process for the shank
%       raw_channels_from_tip               - distance of each `channel_ids_to_process` from the probe tip
%       interp_points_from_tip              - space to interpolate over (`interp_distance_um` resolution)
%       run                                 - main method for running full analysis
%       compute_power_on_all_channels       - computes the bandpower on all channels
%       interpolate_across_channels         - interpoates the bandpower across the channels at `interp_distance_um` resolution
%       smooth_across_channels              - smooth the interpolated profile from `interpolate_across_channels`
%       bandpower                           - computes the bandpower for a channel
%       search_for_l5                       - searchs the smoothed powerfor the peak
%       separate_into_electrode_columns     - separates the full data into electrode columns
%       remove_unuseful_channels            - removes reference and bad channels from the data
%       set_unuseful_channels_to_nan        - sets reference and bad channels to nan (rather than removing)
%       get_parameters                      - return a structure with the analysis parameters
%       load_saved_data                     - loads previously saved data
%       props_to_save                       - which properties of this class to save
%       apply_parameters                    - applies parameters from a structure
%       interactive_plot                    - creates interactive plot to view the averaged power profile with anatomical boundaries overlaid
%       plot_raw_batches                    - plots the power across all channels separately for each batch in a separate subplot
%       overlay_raw_batches                 - plots power across all channels for each batch on the same plot
%       plot_interpolated_batches           - similar to `plot_raw_batches` but plots the power interpolated across 
%                                             banks of electrodes instead of every channel
%       plot_summary                        - create figure with 3 axes showing various features of the power profile

    properties
        
        first_batch_s = 60   % seconds
        batch_interval_s = 300
        batch_duration_s = 10
        n_batches_to_process = inf
        
        lower_hz = 500
        upper_hz = 5000
        
        smoothing_distance_um = 50
        interp_distance_um = 1
        
        channel_ids_to_remove
        
        search_above = 1000  % um from probe tip
        search_below = inf
        
        batches_to_use = 1:10;
        
        power_raw
        
        probe_track
        clusters_from_tip_um
        
        ephys_l5
        anatomy_l5
        delta_l5
        
        probe_id
    end
    
    properties (SetAccess = private)
        
        shank_id
        
        power_interp
        power_smooth
        
        rec
    end
    
    properties (Dependent = true, Hidden = true)
        
        bad_electrode_ids
        bad_channel_ids
        sample_start_t
    end
    
    
    
    methods
        
        function obj = HighFrequencyPowerProfile(recording, probe_id, shank_id)
        % HighFrequencyPowerProfile
        %
        %   HighFrequencyPowerProfile(RECORDING, PROBE_ID, SHANK_ID)
        %   class for computing the high frequency ("MUA spectral") power
        %   profile along a probe shank (c.f. Senzai et al., 2018). 
        %
        %   Inputs:
        %       RECORDING - an object of class SpikeGLXRecording
        %       PROBE_ID -  string with the probe recording ID
        %       SHANK_ID -  integer with the shank to compute the power on (zero
        %                   indexed) (default = 0)
        %
        %   See also: SpikeGLXRecording
        
            VariableDefault('shank_id', 0);
            
            obj.rec = recording;
            
            obj.upper_hz = min(obj.upper_hz, obj.rec.fs/2);
            
            obj.probe_id = probe_id;
            obj.shank_id = shank_id;
            
            obj.batches_to_use = 1:min(length(obj.sample_start_t), 10);
        end
        
        
        
        function channel_ids = channel_ids_to_process(obj)
        %%channel_ids_to_process ID of the channels to process for the shank
        %
        %   CHANNEL_IDS = channel_ids_to_process()
        %   get channel IDs for which we will compute the raw power, by default
        %   we process all channels along a shank (excluding sync
        %   channels). Returned as a vector CHANNEL_IDS.
            
            % get all channel IDs for the current shank (this is ordered
            % from bottom of the shank to top
            channel_ids = obj.rec.get_channel_ids_along_shank(obj.shank_id);
        end
           
        
        
        function val = get.bad_electrode_ids(obj)
        %%list of bad electrodes known for some probes
            
            val = [];
            % probe serial number
            if isfield(obj.rec.config, 'imProbeSN')
                if strcmp(obj.rec.config.imProbeSN, '641250910')
                    val = [0, 0;
                           1, 0;
                           2, 0;
                           4, 0;
                           26, 0;
                           71, 0;
                           97, 0;
                           104, 0;
                           123, 0;
                           151, 0;
                           155, 0;
                           156, 0;
                           167 0];
                elseif strcmp(obj.rec.config.imProbeSN, '641250913')
                    val = [0, 0;
                           1, 0;
                           2, 0;
                           3, 0;
                           4, 0];
%                 elseif strcmp(obj.rec.config.imProbeSN
                end
            elseif isfield(obj.rec.config, 'imDatPrb_sn')
                if strcmp(obj.rec.config.imDatPrb_sn, '19011119353')
                    val = [69, 1;
                           98, 1];  
                end
            end
        end
        
        
        
        function val = get.bad_channel_ids(obj)
        %%from the bad electrode IDs, identifies the channel IDs which will be bad
        
            channel_ids = obj.rec.all_channel_ids();
            [electrode_ids, shank_ids] = obj.rec.electrode_id_from_channel_id(channel_ids);
            idx = ismember(electrode_ids, obj.bad_electrode_ids(:, 1)) & ismember(shank_ids, obj.bad_electrode_ids(:, 2));
            val = channel_ids(idx);
        end
        
        
        
        function val = get.ephys_l5(obj)
        %property just calling `search_for_l5` getting distance of peak of
        %HF power from the probe tip
        
            val = obj.search_for_l5();
        end
        
        
        
        function val = get.anatomy_l5(obj)
        %%returns the distance of the middle of VISp5 from the probe tip in um
        
            val = nan;
            if ~isempty(obj.probe_track)
                val = obj.probe_track.mid_l5_visp();
            end
        end
        
        
        
        function val = get.delta_l5(obj)
        %%returns the difference between `ephys_l5` and `anatomy_l5` in um
        
            val = 0;
            if ~isnan(obj.anatomy_l5)
                val = obj.ephys_l5 - obj.anatomy_l5;
            end
        end
        
        
        
        function val = raw_channels_from_tip(obj)
        %%raw_channels_from_tip Distance of each `channel_ids_to_process` from the probe tip
        %
        %   DISTANCE = raw_channels_from_tip()
        %   returns the distance in microns of each channel from the probe
        %   tip. The distance is computed for each channel returned by
        %   `channel_ids_to_process`
        %
        %   See also: `channel_ids_to_process`
        
            channel_ids = obj.channel_ids_to_process();
            val = obj.rec.channel_from_tip_um(channel_ids);
        end
        
        
        
        function val = interp_points_from_tip(obj)
        %%interp_points_from_tip Space to interpolate over (`interp_distance_um` resolution)
        %
        %   DISTANCE = interp_points_from_tip()
        %   returns a spatial map at `interp_distance_um` resolution from
        %   the smallest to the largest distance returned by
        %   `raw_channels_from_tip`. This map is used to interpolate the HF
        %   power profile at a higher spatial resolution.
        
            % we will interpolate at resolution
            from_tip = obj.raw_channels_from_tip();
            val = from_tip(1) : obj.interp_distance_um : from_tip(end);
        end
        
        
        
        function val = get.sample_start_t(obj)
        %%start times of each batch
        
            max_t = (obj.rec.n_samples/obj.rec.fs) - obj.batch_duration_s;
            val = obj.first_batch_s : obj.batch_interval_s : max_t;
            val = val(1:min(obj.n_batches_to_process, length(val)));
        end
        
        
        
        function run(obj)
        %%run Main method for running full analysis
        %
        %   run()
        %   computes all steps in the HF power pipeline:
        %       1. compute raw power on all channels
        %       2. interpolation over channels at obj.interp_distance_um
        %       resolution
        %       3. smooth across channels with obj.smoothing_distance_um
        %       sliding window
            
            obj.compute_power_on_all_channels();
            obj.interpolate_across_channels();
            obj.smooth_across_channels();
        end
        
        
        
        function compute_legacy(obj)
        %%computes all steps in the HF power pipeline (using legacy
        %   approach) This is just here to replicate old values
        
            obj.compute_power_on_all_channels_legacy();
            obj.interpolate_across_channels_legacy();
            obj.smooth_across_channels_legacy();
        end
        
        
        
        function compute_power_on_all_channels(obj, channel_ids)
        %%compute_power_on_all_channels Computes the bandpower on all channels
        %
        %   compute_power_on_all_channels(CHANNEL_IDS)
        %   computes HF power on channels with ID CHANNEL_IDS. 
        %   If CHANNEL_IDS is not supplied or is empty, all channels
        %   returned by `channel_ids_to_process` will be processed.
        %   The rest of the module will not be useable if CHANNEL_IDS is
        %   supplied, so this optional argument is only useful if you really just want
        %   the raw power on a set of channels.
        
            VariableDefault('channel_ids', []);
            
            if isempty(channel_ids)
                channel_ids = obj.channel_ids_to_process();
            end
            
            obj.power_raw = nan(length(channel_ids), length(obj.sample_start_t));
            
            fprintf('Computing high frequency power for %s\n', obj.rec.bin_fname);
            
            for jj = 1 : length(obj.sample_start_t)
                
                str1 = sprintf('  Batch %i/%i\n', jj, length(obj.sample_start_t));
                fprintf(str1);
                
                start_t     = obj.sample_start_t(jj);
                end_t       = obj.sample_start_t(jj) + obj.batch_duration_s;
                    
                for ii = 1 : length(channel_ids)
                    
                    str2 = sprintf('   Channel %i/%i', ii, length(channel_ids));
                    fprintf(str2);
                    
                    channel_id = channel_ids(ii);
                    obj.power_raw(ii, jj) = obj.bandpower(channel_id, start_t, end_t);
                    fprintf(repmat('\b', 1, length(str2)));
                end
                
                fprintf(repmat('\b', 1, length(str1)));
            end
        end
        
        
        
        function compute_power_on_all_channels_legacy(obj)
        %%legacy version of computing raw power on all channels
            
            channel_ids = 3:obj.rec.n_saved_channels-2;
            channel_ids(ismember(channel_ids, obj.rec.reference_channel_ids)) = [];
            
            obj.power_raw = nan(obj.rec.n_saved_channels-1, length(obj.sample_start_t));
            
            for ii = 1 : length(channel_ids)
                
                channel_id = channel_ids(ii);
                
                for jj = 1 : length(obj.sample_start_t)
                    
                    start_t     = obj.sample_start_t(jj);
                    end_t       = obj.sample_start_t(jj) + obj.batch_duration_s;
                    
                    obj.power_raw(channel_ids(ii)+1, jj) = obj.bandpower(channel_id, start_t, end_t);
                end
            end
        end
        
        
        
        function interpolate_across_channels(obj)
        %%interpolate_across_channels Interpoates the power across channels at `interp_distance_um` resolution
        %
        %   interpolate_across_channels()
        %   uses the power computation from `compute_power_on_all_channels`,
        %   separates the data into electrode columns, and
        %   interpolates each column at obj.interp_distance_um resolution.
            
            channel_ids = obj.channel_ids_to_process;
            
            
            % remove bad and reference channels
            [power_raw, channel_ids] = obj.remove_unuseful_channels(obj.power_raw, channel_ids);
            
            % separate channels into electrode banks
            [power_raw_columns, column_channel_ids] = obj.separate_into_electrode_columns(power_raw, channel_ids);
            
            n_interpolated_points   = length(obj.interp_points_from_tip);
            n_batches               = size(obj.power_raw, 2);
            n_electrode_columns     = length(power_raw_columns);
            
            % preallocate
            obj.power_interp = nan(n_interpolated_points, n_electrode_columns, n_batches);
            
            for ii = 1 : n_batches
                
%                 power_interp = nan(n_interpolated_points, n_electrode_columns);
                
                % interpolate all banks of electrodes
                for jj = 1 : n_electrode_columns
                    
                    these_channels_from_tip_um = obj.rec.channel_from_tip_um(column_channel_ids{jj});
                    obj.power_interp(:, jj, ii) = interp1(these_channels_from_tip_um, power_raw_columns{jj}(:, ii), obj.interp_points_from_tip);
                end
                
%                 % average across electrode banks
%                 obj.power_interp(:, ii) = mean(power_interp, 2);
            end
        end
        
        
        
        function interpolate_across_channels_legacy(obj)
        %%legacy version of interpolation across channels to replicate old
        %   values
            
            good_channel_ids = 0:obj.rec.n_saved_channels-2;
            good_channel_ids([0;1;2;obj.rec.reference_channel_ids]+1)=[];
            
            obj.power_interp = nan(obj.rec.n_saved_channels-1, size(obj.power_raw, 2));
            
            for ii = 1 : size(obj.power_raw, 2)
                
                power = obj.power_raw(:, ii); %#ok<*PROP>
                obj.power_interp(:, ii) = ...
                    interp1(good_channel_ids+1, power(good_channel_ids+1), 1:length(power));
            end
        end
        
        
        
        function smooth_across_channels(obj)
        %%smooth_across_channels Smooth the interpolated profile after `interpolate_across_channels`
        %
        %   smooth_across_channels()
        %   smooth across the interpolated power profile from obj.interpolate_across_channels
        %   with a distance of obj.smoothing_distance_um
            
            obj.power_smooth = nan(size(obj.power_interp));
            n_smoothing_points = ceil(obj.smoothing_distance_um / obj.interp_distance_um);
            for ii = 1 : size(obj.power_interp, 2)
                for jj = 1 : size(obj.power_interp, 3)
                    obj.power_smooth(:, ii, jj) = smooth(obj.power_interp(:, ii, jj), n_smoothing_points);
                end
            end
        end
        
        
        
        function smooth_across_channels_legacy(obj)
        %%legacy version of smoothing across interpolated channels to replicate old
        %%values
            
            obj.power_smooth = nan(size(obj.power_interp));
            
            for ii = 1 : size(obj.power_smooth, 2)
                
                power = obj.power_interp(:, ii);
                obj.power_smooth(:, ii) = smooth(power, 10);
            end
        end
        
        
        
        function bp = bandpower(obj, channel_id, start_t, end_t)
        %%bandpower Computes the bandpower for a channel
        %
        %   bandpower(CHANNEL_ID, START_T, END_T)
        %   compute bandpower between START_T and END_T for channel ID
        %   CHANNEL_ID. Bandpower is computed between `lower_hz` and
        %   `upper_hz`.
        
            data = obj.rec.get_data_between_t(channel_id, [start_t, end_t], [], false);
            bp = bandpower(single(data), obj.rec.fs, [obj.lower_hz, obj.upper_hz]);
        end
        
        
        
        function bp = bandpower_legacy(obj, channel_id, start_t, end_t)
        %%legacy version of power computation
        
            data = obj.rec.get_data_between_t(channel_id, [start_t, end_t], [], true);
            bp = bandpower(single(data), obj.rec.fs, [obj.lower_hz, obj.upper_hz]);
        end
        
        
        
        function l5_from_tip = search_for_l5(obj)
        %%search_for_l5 Searchs the smoothed powerfor the peak
        %
        %   DISTANCE = search_for_l5()
        %   look for a peak in the high frequency power profile contained in
        %   property `power_smooth` and return the distance from the peak.
            
            % search for peak within bounds
            from_tip = obj.interp_points_from_tip;
            search_idx = find(from_tip > obj.search_above & from_tip < obj.search_below);
            
            avg_power = mean(mean(obj.power_smooth(:, :, obj.batches_to_use), 2), 3);
            
            [~, max_idx] = max(avg_power(search_idx));
            
            l5_from_tip = from_tip(search_idx(1) + max_idx - 1);
        end
        
        
        
        function offset = auto_offset(obj)
        %%auto_offset    
            offset = obj.ephys_l5 - obj.probe_track.mid_l5_visp();
        end
        
        
        
        function l5_from_tip = search_for_l5_legacy(obj)
        %legacy version of looking for L5 to replicate old values
            
            power = obj.power_smooth(1:2:end, obj.batches_to_use);
            idx = 1:2:obj.rec.n_saved_channels-1;  % -1 to account for trigger... but TODO
            from_tip = obj.rec.channel_from_tip_um(idx-1);  % -1 because this function is 0 indexed

            new_y = from_tip(1):from_tip(end);
            
            power = mean(bsxfun(@rdivide, power, max(power, [], 1)), 2);
            
            % interpolate the high frequency power at 1um resolution
            hf_power_interp = interp1(from_tip, power, new_y, 'spline');
            
            % search for peak within bounds
            search_idx = find(new_y > obj.search_above & new_y < obj.search_below);
            [~, max_idx] = max(hf_power_interp(search_idx));
            
            l5_from_tip = new_y(search_idx(1) + max_idx - 1);
        end
        
        
        
        function [power_columns, column_channel_ids] = separate_into_electrode_columns(obj, power, channel_ids)
        %%separate_into_electrode_columns Separates the full data into electrode columns
        %
        %   [POWER_COL, CHANNEL_IDS_COL] = separate_into_electrode_columns(POWER, CHANNEL_IDS)
        %   takes POWER (#channels x #batches array) computed on CHANNEL_IDS 
        %   (#channels x 1 vector) and separates it into
        %   electrode columns. Returns POWER_COL a cell array with each entry
        %   being the power data for an electrode column (#channel along
        %   column x # batches) and CHANNEL_IDS_COL  a cell array with the
        %   channels IDs along each column.
        
            % split power_raw into n_columns
            n_columns = obj.rec.n_electrode_columns_per_shank;
            
            % preallocate
            power_columns = cell(1, n_columns);
            column_channel_ids = cell(1, n_columns);
            
            % for each column
            for ii = 1 : n_columns
                
                % get channel IDs along column
                channel_ids_along_column = obj.rec.get_channel_ids_along_column(obj.shank_id, ii-1);
                
                % which of the channels in raw power are along this column
                channel_idx = ismember(channel_ids, channel_ids_along_column);
                
                % separate
                power_columns{ii} = power(channel_idx, :);
                column_channel_ids{ii} = channel_ids(channel_idx);
            end
        end
        
        
        
        function [power, channel_ids] = remove_unuseful_channels(obj, power, channel_ids)
        %%remove_unuseful_channels Removes reference and bad channels from the data
        %
        %   [POWER_RM, CHANNEL_IDS_RM] = remove_unuseful_channels(POWER, CHANNEL_IDS)
        %   takes POWER (#channels x #batches array) computed on CHANNEL_IDS 
        %   (#channels x 1 vector) and removes channels which are
        %   reference, identified as bad or otherwise flagged to remove.
        %       These channel IDs are specified in:
        %           SpikeGLXRecording.reference_channel_ids
        %           `bad_channel_ids`
        %           `channel_ids_to_remove`
        %
        %   Retruns POWER_RM (#channels after removal x # batches array)
        %   and CHANNEL_IDS_RM a vector of channels IDs after removal.
            
                channel_idx = ~ismember(channel_ids, [obj.rec.reference_channel_ids; obj.bad_channel_ids(:); obj.channel_ids_to_remove(:)]);
                channel_ids   = channel_ids(channel_idx);
                power         = power(channel_idx, :);
        end
        
        
        
        function power = set_unuseful_channels_to_nan(obj, power, channel_ids)
        %%set_unuseful_channels_to_nan Sets reference and bad channels to nan (rather than removing)
        %
        %   POWER_NAN = remove_unuseful_channels(POWER, CHANNEL_IDS)
        %   takes POWER (#channels x #batches array) computed on CHANNEL_IDS 
        %   (#channels x 1 vector) and sets channels which are
        %   reference, identified as bad or otherwise flagged to remove as
        %   NaN values.
        %       
        %   These channel IDs are specified in:
        %           SpikeGLXRecording.reference_channel_ids
        %           `bad_channel_ids`
        %           `channel_ids_to_remove`
        %
        %   Retruns POWER_NAN (#channels x # batches array) with the nan
        %   values.
            
                channel_idx             = ismember(channel_ids, [obj.rec.reference_channel_ids; obj.bad_channel_ids(:); obj.channel_ids_to_remove(:)]);
                power(channel_idx, :)   = nan;
        end
        
        
        
        function params = get_parameters(obj)
        %%get_parameters Return a structure with the analysis parameters
        %
        %   STRUCT = get_parameters()
        %   get parameters to save with the HF power profiles and offset
        %   file. Return these in a structure STRUCT.
        
            props = obj.props_to_save();
            
            for ii = 1 : length(props)
                params.(props{ii}) = obj.(props{ii});
            end
        end
        
        
        
        function load_saved_data(obj, mat_fname)
        %%load_saved_data Loads previously saved data
        %
        %   load_saved_data(FILENAME) loads in a .mat file which was
        %   previously saved. Inserts the parameters into the
        %   HighFrequencyPowerProfile object.
        
            props = obj.props_to_save();
            
            params = load(mat_fname, props{:});
            
            % make sure data came from this probe and shank originally
            assert(strcmp(params.probe_id, obj.probe_id) &  params.shank_id == obj.shank_id); %#ok<*PROPLC>
            
            dependent_props = {'bad_electrode_ids', ...
                               'bad_channel_ids', ...
                               'ephys_l5', ...
                               'anatomy_l5', ...
                               'delta_l5', ...
                               'sample_start_t', ...
                               'channel_ids_to_process', ...
                               'raw_channels_from_tip', ...
                               'interp_points_from_tip'};
            
            for ii = 1 : length(props)
                if any(strcmp(props{ii}, dependent_props))
                    continue
                end
                if strcmp(props{ii}, 'power_raw')
                    disp('');
                end
                obj.(props{ii}) = params.(props{ii});
            end
        end
        
        
        
        function props = props_to_save(obj)
        %%props_to_save Which properties of this class to save
        %
        %   PROPERTIES = props_to_save() returns a cell array of strings
        %   with the properties of the class to save.
        
            props = {'probe_id', ...
                     'shank_id', ...
                     'first_batch_s', ...
                     'batch_interval_s', ...
                     'batch_duration_s', ...
                     'n_batches_to_process', ...
                     'sample_start_t', ...
                     'lower_hz', ...
                     'upper_hz', ...
                     'smoothing_distance_um', ...
                     'interp_distance_um', ...
                     'channel_ids_to_remove', ...
                     'bad_channel_ids', ...
                     'channel_ids_to_process', ...
                     'search_above', ...
                     'search_below', ...
                     'batches_to_use', ...
                     'ephys_l5', ...
                     'anatomy_l5', ...
                     'delta_l5', ...
                     'power_raw', ...
                     'power_interp', ...
                     'power_smooth', ...
                     'raw_channels_from_tip', ...
                     'interp_points_from_tip'};
                 
        end
        
        
        
        function apply_parameters(obj, params)
        %%apply_parameters Applies parameters from a structure
        %
        %   apply_parameters(STRUCT)
        %   takes a structure STRUCT, which contains fields with the same
        %   names as the properties in this class and sets the properties
        %   to have the values that the fields have.
            
            if isempty(params); return; end
            
            props = fieldnames(params);
            for ii = 1 : length(props)
                obj.(props{ii}) = params.(props{ii});
            end
        end
        
        
        
        function hf_plot = interactive_plot(obj)
        %%interactive_plot Creates interactive plot to view the averaged power profile with anatomical boundaries overlaid
        %
        %   Wrapper around HighFrequencyPowerProfilePlot.plot_by_column_interactive
        %
        %   See also: HighFrequencyPowerProfilePlot
        
            hf_plot = HighFrequencyPowerProfilePlot(obj);
            hf_plot.plot_by_column_interactive();
        end
        
        
        
        function plot_raw_batches(obj)
        %%plot_raw_batches Plots the power across all channels separately for each batch in a separate subplot
        %
        %   Wrapper around HighFrequencyPowerProfilePlot.raw_batches
        %
        %   See also: HighFrequencyPowerProfilePlot
        
            hf_plot = HighFrequencyPowerProfilePlot(obj);
            hf_plot.raw_batches();
        end
        
        
        
        function overlay_raw_batches(obj)
        %%overlay_raw_batches Plots power across all channels for each batch on the same plot
        %
        %   Wrapper around HighFrequencyPowerProfilePlot.overlay_raw_batches 
        %
        %   See also: HighFrequencyPowerProfilePlot
        
            hf_plot = HighFrequencyPowerProfilePlot(obj);
            hf_plot.overlay_raw_batches();
        end
        
        
        
        function plot_interpolated_batches(obj)
        %%plot_interpolated_batches Similar to `plot_raw_batches` but plots the power interpolated across banks of electrodes instead of every channel
        %
        %   Wrapper around HighFrequencyPowerProfilePlot.interpolated_batches 
        %
        %   See also: HighFrequencyPowerProfilePlot
        
            hf_plot = HighFrequencyPowerProfilePlot(obj);
            hf_plot.interpolated_batches();
        end
        
        
        
        function h_fig = plot_summary(obj)
        %%plot_summary Create figure with 3 axes showing various features of the power profile
        %
        %   Wrapper around HighFrequencyPowerProfilePlot.plot_summary   
        %
        %   See also: HighFrequencyPowerProfilePlot
        
            hf_plot = HighFrequencyPowerProfilePlot(obj);
            hf_plot.plot_summary();
            h_fig = gcf;
        end
    end
end
