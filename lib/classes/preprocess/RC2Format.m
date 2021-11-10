classdef RC2Format < RC2Analysis
    
    properties
    end
    
    
    methods
        
        function obj = RC2Format()
            
            obj = obj@RC2Analysis();
        end
        
        
        
        function format(obj, probe_id)
            
            formatted_data.probe_id                             = probe_id;
            
            shank_ids                                           = obj.get_shank_ids(probe_id);
            
            for ii = 1 : length(shank_ids)
                formatted_data.anatomy(ii)                      = obj.format_anatomy(probe_id, shank_ids(ii)); %#ok<*NASGU>
            end
            
            formatted_data.clusters                             = obj.format_clusters(probe_id);
            formatted_data.sessions                             = obj.format_sessions(probe_id);
            [probe_t, formatted_data.n_triggers, trigger_t]     = obj.synchronize(probe_id, formatted_data.sessions); %#ok<*ASGLU>
            formatted_data.selected_clusters                    = obj.load.selected_clusters(probe_id);
            formatted_data.sessions                             = obj.insert_timebase(formatted_data.sessions, probe_t, trigger_t);
            
            % correct for acquisition mistakes
            formatted_data.sessions                             = data_cleansing(formatted_data.sessions);
            
            obj.save.formatted_data(probe_id, formatted_data)
        end
        
        
        
        function anatomy = format_anatomy(obj, probe_id, shank_id)
        %%reads in the probe track and offset.txt file and calculates the
        %%boundaries of each region on the probe
            
            % read probe file
            probe_track                 = obj.load_track(probe_id, shank_id);
            track_offset                = obj.load.track_offset(probe_id, shank_id);
            
            anatomy.probe_id            = probe_id;
            anatomy.offset              = track_offset;
            anatomy.region_boundaries   = [];
            anatomy.region_id           = [];
            anatomy.region_str          = [];
            
            if isempty(probe_track); return; end
            
            if ~isempty(track_offset)
                probe_track.offset      = track_offset;
            end
            
            % run 'adjusted' to account for the shift
            [anatomy.region_boundaries, anatomy.region_id, anatomy.region_str] = probe_track.region_boundaries_adjusted();
        end
        
        
        
        function probe_track = load_track(obj, probe_id, shank_id)
        %%loads track info and processes it a bit
            
            probe_track_tbl     = obj.load.track_csv(probe_id, shank_id);
            
            if isempty(probe_track_tbl); probe_track = []; return; end
            
            probe_track = ProbeTrack(probe_track_tbl, obj.get_probe_type_from_metadata(probe_id));
        end
        
        
        
        function probe_track = load_track_with_offset(obj, probe_id, shank_id)
        %%load the probe track file and apply the offset in the offset.txt
        %%file
            % read probe file
            probe_track         = obj.load_track(probe_id, shank_id);
            track_offset        = obj.load.track_offset(probe_id, shank_id);
            
            if ~isempty(track_offset)
                probe_track.offset = track_offset;
            end
        end
        
        
        
        function clusters = format_clusters(obj, probe_id)
        %%FORMAT_CLUSTERS
        %   CLUSTERS = FORMAT_CLUSTERS(probe_id)
        %
        %   Takes the files in the output of kilosort2/ecephys_spike_sorting or
        %   related fork, and produces a sequence of "cluster" objects with all the
        %   associated information.
        %
        %   Input:  ks_dir - directory of kilosort2 output
        %   Output: clusters - structure array of clusters with several properties
            
            % read the kilosort data
            params            = obj.load.params(probe_id);
            spike_clusters    = obj.load.ks2_npy(probe_id, 'spike_clusters');
            spike_templates   = obj.load.ks2_npy(probe_id, 'spike_templates');
            spike_times       = obj.load.ks2_npy(probe_id, 'spike_times');
            amplitudes        = obj.load.ks2_npy(probe_id, 'amplitudes');
            templates         = obj.load.ks2_npy(probe_id, 'templates');
            chan_map          = obj.load.ks2_npy(probe_id, 'channel_map');
            chan_pos          = obj.load.ks2_npy(probe_id, 'channel_positions');
            qm_table          = obj.load.metrics_csv(probe_id);
            waveform_fixed_table = obj.load.waveform_metrics_fixed_csv(probe_id);
            cluster_groups    = obj.load.cluster_groups(probe_id);
            ks_label          = obj.load.ks2_label(probe_id);
            
            spikeglx_meta       = obj.load.spikeglx_ap_metadata(probe_id);
            
            shank_ids         = obj.get_shank_ids(probe_id);
            
            % load available probe tracks
            for ii = length(shank_ids) : -1 : 1
                probe_track{ii}       = obj.load_track_with_offset(probe_id, shank_ids(ii));
            end
            
            
            % which quality metrics to unpack
            qm = {'firing_rate', 'presence_ratio', 'isi_viol', 'amplitude_cutoff', ...
                'isolation_distance', 'l_ratio', 'd_prime', 'nn_hit_rate', ...
                'nn_miss_rate', 'silhouette_score', 'max_drift', 'cumulative_drift'};
            
            % all clusters with spikes
            cluster_ids = sort(unique(spike_clusters));
            
            % iterate over clusters
            for i = length(cluster_ids) : -1 : 1
                
                % location in the total spike vector of this cluster
                spike_mask = spike_clusters == cluster_ids(i);
                
                clusters(i).id                  = cluster_ids(i); % cluster ID
                clusters(i).spike_times         = double(spike_times(spike_mask)) / params.sample_rate; % times of spikes in probe
                clusters(i).spike_sample_point  = spike_times(spike_mask); % sample points of spikes on probe
                clusters(i).amplitudes          = amplitudes(spike_mask); % amplitude of each spike
                
                % most numerous template before merging and splitting
                %   required to get template and KS properties which are not modified
                %   by phy... more for merged units (as split units will all have the
                %   same cluster_id beforehand.
                ori_cluster_id                  = mode(spike_templates(spike_mask));
                
                % get template of this cluster
                %   TODO: should we unwhiten the templates?
                template                        = permute(templates(ori_cluster_id+1, :, :), [2, 3, 1]); % +1 because id is 0 indexed
                
                % get amplitude of waveform on each channel
                amp                             = max(template, [], 1) - min(template, [], 1);
                
                % get channel of max amplitude
                [~, max_chan_good]              = max(amp);
                
                % index correctly
                max_chan                        = chan_map(max_chan_good);
                y_dist                          = chan_pos(max_chan_good, 2);
                [~, shank_id]                   = spikeglx_meta.electrode_id_from_channel_id(max_chan);
                shank_idx                       = find(shank_ids == shank_id);
                
                % store best channel and distance from the probe tip
                clusters(i).best_channel        = max_chan;
                clusters(i).shank_id            = shank_id;
                
                probe_type = obj.get_probe_type_from_metadata(probe_id);
                probe = Probe(probe_type);
                
                clusters(i).distance_from_probe_tip = probe.tip_length_um + y_dist;
                
                % find position in the probe track file
                if ~isempty(probe_track{shank_idx})
                    
                    [clusters(i).region_id, ...
                     clusters(i).region_str, ...
                     clusters(i).depth] = ...
                        probe_track{shank_idx}.get_region_of_point_from_tip_adjusted(clusters(i).distance_from_probe_tip);
                else
                    clusters(i).region_id   = nan;
                    clusters(i).region_str  = "";
                    clusters(i).depth       = nan;
                end
                
                % cluster groups info
                %   this is changed by Phy so we can use the current cluster
                if ~isempty(cluster_groups)
                    cluster_groups_idx = cluster_groups.cluster_id == cluster_ids(i);
                    clusters(i).class = cluster_groups.group{cluster_groups_idx};
                else
                    clusters(i).class = '';
                end
                
                % KS2 metrics
                if ~isempty(ks_label)
                    % use the original template as this is NOT changed after phy
                    %   takes the original cluster with most spikes contributing merged
                    %   cluster
                    %   if split, the original cluster will be the same label for all new
                    %   clusters.
                    ks_label_idx = ks_label.cluster_id == ori_cluster_id;
                    clusters(i).ks2_class = ks_label.KSLabel{ks_label_idx};
                else
                    clusters(i).ks2_class = nan;
                end
                
                % if quality metrics were read successfully
                for j = 1 : length(qm)
                    if ~isempty(qm_table)
                        % metric index
                        metric_idx = qm_table.cluster_id == cluster_ids(i);
                        clusters(i).(qm{j}) = qm_table.(qm{j})(metric_idx);
                    else
                        clusters(i).(qm{j}) = nan;
                    end
                end
                
                % additional metrics
                median_amp = median(clusters(i).amplitudes);
                min_amp = min(clusters(i).amplitudes);
                clusters(i).amplitude_ratio = median_amp/min_amp;
            end
            
            waveform_metrics = {'peak_channel', 'snr', 'duration', 'halfwidth', 'PT_ratio', ...
                'repolarization_slope', 'recovery_slope', 'amplitude', 'spread', ...
                'velocity_above', 'velocity_below'};
            
            % if the waveform fixed table doesn't exist
            if isempty(waveform_fixed_table)
                % create empty entries for each of the waveform properties
                for ii = 1 : length(clusters)
                    for jj = 1 : length(waveform_metrics)
                        clusters(ii).(waveform_metrics{jj}) = [];
                    end
                    clusters(ii).waveform_fixed = false;
                end
            else
                % otherwise fill in the waveform information for each
                % cluster
                for ii = 1 : length(clusters)
                    tbl_idx = waveform_fixed_table.cluster_id == clusters(ii).id;
                    for jj = 1 : length(waveform_metrics)
                        clusters(ii).(waveform_metrics{jj}) = waveform_fixed_table.(waveform_metrics{jj})(tbl_idx);
                    end
                    clusters(ii).waveform_fixed = true;
                end
            end
        end
        
        
        
        function sessions = format_sessions(obj, probe_id)
        %%formats data in RC2 bin files
            
            session_list = obj.get_session_ids_list(probe_id);
            
            for ii = 1 : length(session_list)
                
                sessions(ii) = obj.format_session(session_list{ii});
            end
        end
        
        
        
        function session = format_session(obj, session_id)
        %%format data in RC2 bin file and cameras
            
            % read the bin and cfg file
            [data, dt, chan_names, config] = obj.load.rc2_bin(session_id);
            camera_ids = obj.get_camera_ids(session_id);
            
            session = obj.empty_session_struct(chan_names, camera_ids);
            
            % get filename without full path or extension
            session.probe_id = obj.get_probe_id_from_session_id(session_id);
            session.session_id = session_id;
            
            session.fs = 1/dt;
            session.config = config;
            
            session.n_samples = size(data, 1);
            session.n_channels = size(data, 2);
            
            % unpack data matrix
            for j = 1 : length(chan_names)
                session.(chan_names{j}) = data(:, j);
            end
            
            for ii = 1 : length(camera_ids)
                session.(camera_ids{ii}) = obj.load.camera_csv(session_id, camera_ids{ii});
            end
            
            session.trials = obj.format_trials(session);
            session.n_trials = length(session.trials);
            
            session.rc2_t = (0:session.n_samples-1)' * dt;
        end
        
        
        
        function [probe_t, n_triggers, trigger_t] = synchronize(obj, probe_id, sessions)
        %%synchronize probe recording and RC2 session
            
            expected_min_interval = 2;  % s
            
            params              = obj.load.params(probe_id);
            trigger             = obj.load.trigger_mat(probe_id);
            trigger_ups         = find(diff(trigger) >= 1) / params.sample_rate;  % no need to add 1 to find output
            trigger_interval    = diff(trigger_ups);
            trigger_boundaries  = [1, find(trigger_interval > expected_min_interval) + 1, length(trigger_ups)+1];
            n_sessions          = length(trigger_boundaries)-1;
            
            assert(n_sessions == length(sessions), ['Number of AI recordings provided doesn''t ', ...
                'match expected number from trigger channel'])
            
            % for each expected AI recording
            n_triggers_per_session = nan(n_sessions, 1);
            for i = 1 : n_sessions
                n_triggers_per_session(i) = trigger_boundaries(i+1) - trigger_boundaries(i);
            end
            
            % preallocate
            probe_t             = cell(1, n_sessions);
            trigger_t           = cell(1, n_sessions);
            
            
            for ii = 1 : n_sessions
                
                % get times of triggers for this session
                trigger_idx         = trigger_boundaries(ii) + (0:n_triggers_per_session(ii) - 1)';
                trigger_t{ii}       = trigger_ups(trigger_idx);
                
                trigger_interval    = sessions(ii).config.nidaq.co.low_samps + sessions(ii).config.nidaq.co.high_samps;
                
                n_triggers_expected = ceil(sessions(ii).n_samples / trigger_interval);
                
                start_t             = trigger_ups(trigger_boundaries(ii));
                
                if n_triggers_expected > n_triggers_per_session(ii)
                    
                    end_t   = trigger_ups(trigger_boundaries(ii) + n_triggers_per_session(ii) - 1);
                    
                    % the number of samples acquired up to the rise of the very last
                    % trigger sent during *acquisition*
                    n_samples_accounted_for = (n_triggers_per_session(ii) - 1) * trigger_interval;
                else
                    
                    % get the time of the last trigger sent *during acquisition* measured
                    % on the probe
                    end_t = trigger_ups(trigger_boundaries(ii) + n_triggers_expected - 1);
                    
                    % the number of samples acquired up to the rise of the very last
                    % trigger sent during *acquisition*
                    n_samples_accounted_for = (n_triggers_expected - 1) * trigger_interval;
                    
                    % the number of samples after the rise of the very last trigger
                    err = sessions(ii).n_samples - n_samples_accounted_for;
                    
                    % make sure that this is less that the trigger interval
                    assert(err >= 0 && err < trigger_interval, '??');
                end
                
                f = sessions(ii).n_samples / n_samples_accounted_for;
                
                % the duration of the AI acquisition is f*(end_t - start_t)
                probe_t{ii} =  linspace(start_t, start_t + f * (end_t - start_t), sessions(ii).n_samples)';
                
                n_triggers(ii).rc       = n_triggers_expected;
                n_triggers(ii).probe    = n_triggers_per_session(ii);
            end
        end
    end
    
    
    
    methods (Static = true)
        
        function session = empty_session_struct(chan_names, camera_ids)
            
            session = struct('probe_id', [], ...
                             'session_id', [], ...
                             'fs', [], ...
                             'config', [], ...
                             'n_samples', [], ...
                             'n_channels', [], ...
                             'n_trials', [], ...
                             'trials', [], ...
                             'rc2_t', [], ...
                             'probe_t', [], ...
                             'camera_t', [], ...
                             'camera0', [], ...
                             'camera1', []);
            
            % unpack the channels
            for ii = 1 : length(chan_names)
                session.(chan_names{ii}) = [];
            end
            
            for ii = 1 : length(camera_ids)
                session.(camera_ids{ii}) = [];
            end
        end
        
        
        
        function trials = format_trials(session)
        %%takes a session structure and finds the 'trials' - which for now
        %%is very dependent on the structure of the solenoid signal.
            
            % make sure that the solenoid starts high in the session
            assert(session.solenoid(1) > 2.5, 'Solenoid does not start high');
            
            % get the sample points of solenoid up and down
            solenoid_down = find(diff(session.solenoid > 2.5) == -1) + 1;
            
            % TODO: should be solved at acquisition
            if ~isfield(session.config, 'prot')
                solenoid_down = 1;
                session.config.prot.type = 'Visual';
            end
            
            % make sure that the number of trials from the solenoid is equal to the number of protocols
            % in the config file
            assert(length(solenoid_down) == length(session.config.prot) || ...
                length(solenoid_down) == length(session.config.prot) - 1, ...
                'Number of trials from solenoid does not equal number of trials saved in config file.');
            
            % add an extra bound for the last trial
            solenoid_down(end+1) = length(session.solenoid);
            
            
            for ii = 1 : length(solenoid_down)-1
                
                protocol = session.config.prot(ii).type;
                
                switch protocol
                    case {'Coupled', 'EncoderOnly', 'CoupledMismatch', 'EncoderOnlyMismatch'}
                        chop_idx = (solenoid_down(ii) - 4*session.fs) : solenoid_down(ii+1);
                    case {'StageOnly', 'ReplayOnly'}
                        chop_idx = solenoid_down(ii) : solenoid_down(ii+1);
                    case 'Visual'
                        chop_idx = 1 : length(session.solenoid);
                end
                
                trials(ii).probe_id      = session.probe_id;
                trials(ii).session_id    = session.session_id;
                trials(ii).trial_id      = ii;
                trials(ii).config        = session.config.prot(ii);
                trials(ii).fs            = session.fs;
                trials(ii).protocol      = protocol;
                trials(ii).start_idx     = chop_idx(1);
                trials(ii).end_idx       = chop_idx(end);
            end
        end
        
        
        
        function sessions = insert_timebase(sessions, probe_t, trigger_t)
            
            for ii = 1 : length(sessions)
                sessions(ii).probe_t = probe_t{ii};
                sessions(ii).camera_t = trigger_t{ii};
            end
        end
    end 
end
