classdef CheckClusterQuality < handle
% CheckClusterQuality Class for plotting and checking the quality of a
% cluster.
%
%   CheckClusterQuality Properties:
%       isi_limit       - limits of the ISI histogram plot (default = 30ms)
%       isi_bin         - size of the bins in the ISI histogram plot (default = 0.25ms)
%       isi_thresh      - (default = 1.5ms)
%       isi_min         - (default = 0.166ms)
%       probe_fs        - sample rate of the probe recording (default = 3000Hz)
%       apply_jcolonell_isi_correction - whether to apply the correction
%                                       appearing in a recent version of
%                                       ecephys_spike_sorting (default = false)
%
%     Private:
%       ctl             - instance of RC2Preprocess
%       probe_id        - string with probe ID
%       session_bounds  - bounds in probe time of the different sessions
%       spike_times     - list of spike times for a cluster
%       spike_clusters  - ID of the cluster to which each spike belongs
%       amplitudes      - ampltiude of each spike
%       
%   CheckClusterQuality Methods:
%       plot                        - plots whole figure
%       amplitude_plot              - plots the spike amplitude over time
%       isi_plot                    - plots ISI histogram
%       compute_isi_distribution    - computes the ISIs from a list of spike times
%       compute_isi_info            - computes information to be printed about ISIs and amplitudes
%       amplitude_histogram_plot    - plot histogram of spike amplitudes
%       print_cluster_info          - prints info about ISI calculation
    

    properties
        
        apply_jcolonell_isi_correction = false
        isi_limit = 30; % ms
        isi_bin = 0.25; % ms
        isi_thresh = 1.5; % ms
        isi_min = 0.166; % ms
        probe_fs = 30000; % Hz
    end
    
    properties (SetAccess = private)
        
        ctl
        probe_id
        session_bounds
        
        spike_times
        spike_clusters
        amplitudes
    end
    
    
    
    methods
        
        function obj = CheckClusterQuality(ctl, probe_id)
            
            obj.ctl = ctl;
            obj.probe_id = probe_id;
            
            obj.session_bounds = obj.ctl.get_session_bounds(obj.probe_id);
            obj.spike_times = double(obj.ctl.load.ks2_npy(obj.probe_id, 'spike_times')) / obj.probe_fs;
            obj.spike_clusters = obj.ctl.load.ks2_npy(obj.probe_id, 'spike_clusters');
            obj.amplitudes = obj.ctl.load.ks2_npy(obj.probe_id, 'amplitudes');
        end
        
        
        
        function plot(obj, cluster_id)
        %%plot Plots three axes with cluster quality information
        %
        %   plot(CLUSTER_ID) produces three plots on one figure about
        %   cluster quality. On the left is the amplitude of the spikes of
        %   the cluster over time, with the session recordings delimited.
        %   In the middle is the distribution of spike amplitudes across
        %   the reocrding. On the right is the distribution of ISIs. Also
        %   printed to the command line are some details about the ISI
        %   calculation.
        
            figure('position', [137, 338, 1527, 526]);
            
            h_ax(1) = subplot(1, 3, 1);
            h_ax(2) = subplot(1, 3, 2); hold on;
            h_ax(3) = subplot(1, 3, 3);
            
            spike_idx = obj.spike_clusters == cluster_id;
            these_spike_times = obj.spike_times(spike_idx);
            these_amplitudes = obj.amplitudes(spike_idx);
            isis = obj.compute_isi_distribution(these_spike_times);
            info = obj.compute_isi_info(these_spike_times, these_amplitudes);
            
            obj.print_cluster_info(info);
            obj.amplitude_plot(h_ax(1), these_spike_times, these_amplitudes)
            obj.amplitude_histogram_plot(h_ax(2), these_amplitudes);
            title(h_ax(2), sprintf('Cluster %i', cluster_id));
            obj.isi_plot(h_ax(3), isis);
            title(h_ax(3), sprintf('ISI: %.4f', info.fp_rate));
        end
        
        
        
        function amplitude_plot(obj, h_ax, spike_times, amplitudes)
        %%amplitude_plot Plots the spike amplitude over time 
        %
        %   amplitude_plot(AXIS_HANDLE, SPIKE_TIMES, AMPLITUDES) takes the
        %   spike times (SPIKE_TIMES) and spike amplitudes (AMPLITUDES) for
        %   a cluster, and plots them on axis with handle AXIS_HANDLE.
        %   SPIKE_TIMES and AMPLITUDES should be # spikes x 1 vectors.
        
            % scatter amplitudes across time
            scatter(h_ax, spike_times, amplitudes, 10, 'k', 'fill');
            
            % format
            set(h_ax, 'PlotBoxAspectRatio', [3, 1, 1]);
            xlim(h_ax, [0, max(max(obj.spike_times), max(obj.session_bounds))]);
            ylim(h_ax, [0, max(get(h_ax, 'ylim'))]);
            xlabel(h_ax, 'Time (s)');
            ylabel(h_ax, 'Amplitude (a.u.)');
            
            % get enough colours
            cols = lines(length(obj.session_bounds)/2);
            
            % plot lines demarcating session boundaries
            for ii = 1 : length(obj.session_bounds)
                col = cols(ceil(ii/2), :);
                line(h_ax, obj.session_bounds(ii)*[1, 1], get(h_ax, 'ylim'), ...
                    'color', col, ...
                    'linestyle', '--', ...
                    'linewidth', 2)
            end
        end
        
        
        
        function isi_plot(obj, h_ax, isis)
        %%isi_plot Plot the ISI histogram
        %
        %   isi_plot(AXIS_HANDLE, ISIS) takes the list of inter-spike
        %   intervals (ISIS) and plots them as a histogram on axis with 
        %   handle AXIS_HANDLE.. `isi_thresh` is also plotted on the figure.
        %
        %   See also: compute_isi_distribution
        
            edges = -obj.isi_limit : obj.isi_bin : obj.isi_limit;
            histogram(h_ax, 1e3*[isis, -isis], edges, 'facecolor', 'k');
            
            line(h_ax, [obj.isi_thresh, obj.isi_thresh], get(gca, 'ylim'), 'color', 'r', 'linestyle', '--');
            line(h_ax, -[obj.isi_thresh, obj.isi_thresh], get(gca, 'ylim'), 'color', 'r', 'linestyle', '--');
            
            set(h_ax, 'box', 'off', 'plotboxaspectratio', [1, 1, 1]);
            xlabel(h_ax, 'ISI (ms)');
        end
        
        
        
        function isis = compute_isi_distribution(obj, spike_times)
        %%compute_isi_distribution Computes the ISIs from a list of spike times 
        %
        %   ISIS = compute_isi_distribution(SPIKE_TIMES) computes the ISIs
        %   from the SPIKE_TIMES, and returns it as ISIS.
        
            n_spikes = length(spike_times);
            isis = [];
            N = min(50000, n_spikes);
            I = randperm(n_spikes, N);
            spike_times = spike_times(I);
            n_spikes = N;
            for i = 1 : n_spikes-1
                idx = spike_times - spike_times(i) > 0 & ...
                    spike_times - spike_times(i) < obj.isi_limit;
                isis = [isis; spike_times(idx) - spike_times(i)];
            end
        end
        
        
        
        function info = compute_isi_info(obj, spike_times, amplitudes)
        %%compute_isi_info Computes information to be printed about ISIs
        %%and amplitudes
        %
        %   INFO = compute_isi_info(SPIKE_TIMES, AMPLITUDES) computes
        %   values about the ISI calculation and spike amplitudes. Returns
        %   as a structure INFO with fields:
        %
        %       n_spikes        - the number of spikes
        %       n_viol          - the number of ISI violations
        %       total_rate      - firing rate from first to last spike
        %       viol_rate       - the rate of ISI violations
        %       fp_rate         - ISI violation metric
        %       viol_expected   - number of expected violations from firing rate
        %       min_amp         - minimum spike amplitude
        %       med_amp         - median spike amplitude
        %       amp_ratio       - median / minimum amplitude
        %
        %   If `apply_jcolonell_isi_correction` is true, the `fp_rate`
        %   field is updated with the correction applied by a recent
        %   version of ecephys_spike_sorting.
        %
        %   See also: print_cluster_info
        
            % interspike intervals
            isis = diff(spike_times);
            
            % number of spikes for this cluster
            info.n_spikes = length(spike_times);
            
            % number of violations
            info.n_viol = sum(isis < obj.isi_thresh*1e-3);
            
            % amount of time in which violations could occur
            viol_time = 2*info.n_spikes*(obj.isi_thresh - obj.isi_min)*1e-3;
            
            % spike rate of this cluster
            info.total_rate = info.n_spikes/(max(obj.spike_times) - min(obj.spike_times));
            
            % violation rate
            info.viol_rate = info.n_viol/viol_time;
            
            % ISI violation metric
            info.fp_rate = info.viol_rate/info.total_rate;
            
            info.viol_expected = round(2*1e-3*(obj.isi_thresh-obj.isi_min)*info.n_spikes*info.total_rate);
            
            if obj.apply_jcolonell_isi_correction
                if info.fp_rate < 0.25
                    info.fp_rate = (1 - sqrt(1-4*info.fp_rate))/2;
                else
                    info.fp_rate = 1;
                end
            end
            
            info.min_amp = min(amplitudes);
            info.med_amp = median(amplitudes);
            info.amp_ratio = median(amplitudes)/min(amplitudes);
        end
    end
    
    
    
    methods (Static = true)
        
        function amplitude_histogram_plot(h_ax, amplitudes)
        %%amplitude_histogram_plot Plot histogram of spike amplitudes
        %
        %   amplitude_histogram_plot(AXIS_HANDLE, AMPLITUDES) takes the
        %   list of spike amplitudes and plots them as a histogram to
        %   assess the "amplitude cutoff" metric.
        
            % plot the histogram
            h = histogram(h_ax, amplitudes, 500);
            
            % plot a smooth curve on top
            bin_cent = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
            as = smooth(h.Values, 40);
            plot(h_ax, bin_cent, as, 'linewidth', 2);
            set(h_ax, 'box', 'off');
            xlabel(h_ax, 'Amplitude (a.u.)')
        end
        
        function print_cluster_info(info)
        %%print_cluster_info Prints the information about the clusters
        %
        %   print_cluster_info(INFO) prints the information produced by
        %   `compute_isi_info`.
        %
        %   See also: compute_isi_info
        
            fprintf('Firing rate: %.7f Hz\n', info.total_rate);
            fprintf('# spikes:    %i\n', info.n_spikes);
            fprintf('ISI viol:    %.7f\n', info.fp_rate);
            fprintf(' %i/(2*(isi_threhold-isi_min)*%i*%.5f)\n', info.n_viol, info.n_spikes, info.total_rate);
            fprintf(' # violations expected: %i\n', info.viol_expected)
            fprintf(' # violations observed: %i\n', info.n_viol);
            fprintf('Amplitude minimum:    %.2f\n', info.min_amp);
            fprintf('Amplitude median:    %.2f\n', info.med_amp);
            fprintf('Amplitude ratio:    %.2f\n', info.amp_ratio);
        end
    end
end
