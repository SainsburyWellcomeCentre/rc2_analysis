classdef FiringRate < handle
% FiringRate Class for handling cluster data such as spike times
%
%   FiringRate Properties:
%       spike_times         - #spikes x 1 vector of spike times (in probe timebase)
%       prepad              - when computing spike convolution over a period of time how much before to pad to avoid end effects, in seconds
%       postpad             - when computing spike convolution over a period of time how much after to pad to avoid end effects, in seconds
%       width               - width of the Gaussian to use for convolution, seconds
%       length              - length of the Gaussian window, seconds
%
%   FiringRate Methods:
%       get_convolution             - get continuous spike rate from a set of discrete spike times by convolving with Gaussian
%       get_fr_in_window            - return the firing rate of in a window
%       get_fr_in_multiple_windows  - return the firing rate in multiple windows
%       get_histogram               - return a histogram and bin edges for spiking in a window
%       psth                        - return a histogram and bin edges for spiking around a set of events
%       restrict_times              - restrict the spike times to a window
%       create_spike_train          - create a boolean array with a spike train from list of spike times
%       gauss_filter                - create a Gaussian kernel
%       sample2time                 - takes a sample point on a time base and returns the time
%       time2sample                 - takes a time in a timebase and converts to a sample point
%
%   See also: Cluster

    properties (SetAccess = private)
        
        spike_times
    end
    
    properties
        
        prepad = 1
        postpad = 1
        width = 0.02
        length = 1
    end
    
    
    
    methods
        
        function obj = FiringRate(spike_times)
        % FiringRate
        %
        %   FiringRate(SPIKE_TIMES) creates object for a vector of spike
        %   times (specified in seconds). SPIKE_TIMES is a #spikes x 1
        %   vector.
        
            obj.spike_times = spike_times;
        end
        
        
        
        function r = get_convolution(obj, T)
        %%get_convolution Get continuous spike rate from a set of discrete
        %%spike times by convolving with Gaussian
        %
        %   SPIKE_RATE = get_convolution(TIMEBASE) for a vector specifying
        %   a timebase (#samples x 1 vector), where the timebase is in
        %   seconds, compute the convolved spike rate. SPIKE_RATE is also a
        %   #samples x 1 vector.
        %
        %   e.g. If TIMEBASE is linspace(3, 5, 10e3)', SPIKE_RATE will be a
        %   10e3 x 1 vector with the convolved spike rate between 3 and 5
        %   seconds on the recording.
        
            fs = 1/(T(2) - T(1));
            spike_times = obj.restrict_times([T(1) - obj.prepad, T(end) + obj.postpad]); %#ok<*PROPLC>
            
            [spike_train, sr_t] = obj.create_spike_train(spike_times, [T(1) - obj.prepad, T(end) + obj.postpad], fs);
            gauss = obj.gauss_filter(obj.width, obj.length, fs);
            spike_train = conv(fs*spike_train, gauss, 'same');
            r = interp1(sr_t, spike_train, T);
        end
        
        
        
        function r = get_count(obj, T)
        %%TODO: UNUSED, REMOVE
        
            spike_times = obj.spike_times(obj.spike_times > T(1) & obj.spike_times < T(end));
            spike_points = round(length(T)*(spike_times - T(1)) / (T(end) - T(1))); %#ok<*CPROPLC>
            r = arrayfun(@(x)(sum(spike_points == x)), 1:length(T));
        end
        
        
        
        function [fr, dt, n_spikes] = get_fr_in_window(obj, T)
        %%get_fr_in_window Return the firing rate of in a window
        %
        %   [FIRING_RATE, TIME, N_SPIKES] = get_fr_in_window(TIME_LIMITS)
        %   calculates the number of spikes between the time limits
        %   specified in TIME_LIMITS (a 1x2 array [time_min, time_max]),
        %   and divides by the amount of time within the limits to get a
        %   firing rate, FIRING_RATE. The amount of time between the limits
        %   is also returned in TIME, and the number of spikes found in
        %   N_SPIKES.
        
            assert(length(T) == 2);
            assert(T(2) > T(1));
            
            n_spikes = sum(obj.spike_times >= T(1) & obj.spike_times < T(2));
            dt = T(2)-T(1);
            
            fr = n_spikes / dt;
        end
        
        
        
        function [fr, dt, n_spikes] = get_fr_in_multiple_windows(obj, T)
        %%get_fr_in_multiple_windows Return the firing rate in multiple
        %%windows
        %
        %   [FIRING_RATE, TIME, N_SPIKES] = get_fr_in_multiple_windows(TIME_LIMITS)
        %   calculates the number of spikes between the time limits
        %   specified in TIME_LIMITS (a Nx2 array with each row of the form
        %   [time_min, time_max], and each row is a different window), and
        %   divides by the amount of time within window to give an N x 1
        %   vector of firing rates, FIRING_RATE. The amount of time in each
        %   window is also returned in TIME, and the number of spikes found in
        %   each window returned in N_SPIKES.
        
            assert(size(T, 2) == 2);
            assert(sum(T(:, 2) > T(:, 1)) == size(T, 1));
            assert(~any(T(2:end, 1) < T(1:end-1, 2)));
            
            idx1 = bsxfun(@le, T(:, 1), obj.spike_times(:)');
            idx2 = bsxfun(@ge, T(:, 2), obj.spike_times(:)');
            
            n_spikes = sum(sum(idx1 & idx2, 2));
            dt = sum(T(:, 2)-T(:, 1));
            
            fr = n_spikes / dt;
        end
        
        
        
        function [hstgm, edges] = get_histogram(obj, T, w)
        %%get_histogram Return a histogram and bin edges for spiking in a
        %%window
        %
        %   [COUNTS, EDGES] = get_histogram(TIME_LIMITS, BIN_WIDTH)
        %   computes a histogram between two times in TIME_LIMITS (a 1x2
        %   array of form [time_min, time_max]), in bins of width
        %   BIN_WIDTH. Returns the COUNTS and EDGES of the bins.
        
            % restrict spike times
            st = obj.spike_times(obj.spike_times > T(1) & obj.spike_times < T(2));
            
            % get number of bins
            n_bins = ceil((T(2)-T(1))/w);
            
            % bin edges
            edges = T(1) : w : (T(1) + n_bins*w);
            
            % make sure that did it
            assert(n_bins == length(edges)-1);
            
            % get the histogram
            hstgm = histcounts(st, edges);
            
        end
        
        
        
        function [counts, bin_edges] = psth(obj, stimulus_t, bin_size, window_t)
        %%psth Return a histogram and bin edges for spiking around a set of
        %%events
        %
        %   [COUNTS, EDGES] = get_histogram(EVENT_TIMES, BIN_WIDTH, TIME_LIMITS)
        %   computes a peri-stimulus time histogram between two times in
        %   TIME_LIMITS (a 1x2 array of form [time_min, time_max] in
        %   seconds specifying the time around the events to pool spikes).
        %   EVENT_TIMES is a #events x 1 vector specifying a set of events
        %   of similar nature around which to examine spiking. BIN_WIDTH
        %   determines the width of the bins (seconds). 
        %   Returns the COUNTS and EDGES of the bins.
        
            n_stim = length(stimulus_t);
            
            before_t = window_t(1);
            after_t = window_t(2);
            
            n_bins_before = floor(before_t / bin_size);
            n_bins_after = floor(after_t / bin_size);
            
            bin_edges = sort([0:-bin_size:-bin_size*n_bins_before, bin_size:bin_size:bin_size*n_bins_after]);
            
            N = nan(n_stim, length(bin_edges)-1);
            
            for stim_i = 1 : n_stim
                for bin_i = 1 : length(bin_edges)-1
                    
                    T(1) = stimulus_t(stim_i) + bin_edges(bin_i);
                    T(2) = stimulus_t(stim_i) + bin_edges(bin_i+1);
                    
                    st = obj.spike_times(obj.spike_times > T(1) & obj.spike_times < T(2));
                    N(stim_i, bin_i) = sum(st);
                end
            end
            
            counts = sum(N, 1);
        end
        
        
        
        function st = restrict_times(obj, T)
        %%restrict_times Restrict the spike times to a window
        %
        %   SPIKE_TIMES = restrict_times(TIME_LIMITS)
        %   restricts the spike times to be between TIME_LIMITS (a 1x2
        %   vector of form [time_min, time_max]) and returns them in
        %   SPIKE_TIMES.
        
            st = obj.spike_times(obj.spike_times > T(1) & obj.spike_times < T(end));
        end
        
        
        
        function [spike_train, t] = create_spike_train(obj, spike_times, tlim, fs)
        %%create_spike_train Create a boolean array with a spike train from
        %%list of spike times
        %
        %   [SPIKE_TRAIN, TIMEBASE] = create_spike_train(SPIKE_TIMES, TIME_LIMITS, FS)
        %   takes a list of spike times, SPIKE_TIMES, and a 1x2 vector of start
        %   and end times, TIME_LIMITS, and computes a boolean spike train, sampled
        %   at FS Hz, with true in the spike location and false elsewhere.
        %   #samples x 1 array is returned in SPIKE_TRAIN, and the
        %   corresponding timebase returned in TIMEBASE (#samples x 1).
        %
        %   If spike times occur within the resultion of FS, they *will not be
        %   distinct*.
            
            % pretend that the trace starts at 0, these are the start and end indices
            % of a vector sampled at fs Hz.
            start_idx               = obj.time2sample(tlim(1), fs);
            end_idx                 = obj.time2sample(tlim(2), fs);
            
            % get the indices of each spike
            spike_idx               = obj.time2sample(spike_times, fs) - start_idx + 1;
            % this is the length of a vector starting at tlim(1), ending at tlim(2)
            % and sampled at 'fs'
            len                     = end_idx - start_idx + 1;
            
            % put the spikes in the vector
            spike_train             = false(len, 1);
            spike_train(spike_idx)  = 1;
            
            % the time base of the spike train
            t                       = obj.sample2time(start_idx:end_idx, fs);
        end
    end
    
    
    
    methods (Static = true)
        
        function b = gauss_filter(sigma, sz, fs)
        %%gauss_filter Create a Gaussian kernel
        %
        %   GAUSS = gauss_filter(SIGMA, SIZE, FS)
        %   creates a Gaussian kernel of length SIZE (seconds) and with Gaussian
        %   sigma SIGMA, sampled at FS Hz. Kernel is returned in Gauss
        %   (SIZE*FS x 1) vector.
        
            sigma = sigma*fs;
            sz = sz*fs;    % length of gaussFilter vector
            x = linspace(-sz / 2, sz / 2, sz);
            b = exp(-x .^ 2 / (2 * sigma ^ 2));
            b = b/sum(b); % normalize
        end
        
        
        
        function time = sample2time(sample, fs)
        %%sample2time Takes a sample point on a time base and returns the time
        %
        %   TIME = sample2time(SAMPLE, FS) assumes that a timebase with
        %   sampling rate FS Hz starts on the first sample at time 0, the
        %   second sample at time 1/FS, etc. It then converts sample number
        %   to the corresponding time, TIME.
        
            time = (double(sample) - 1) * (1/fs);
        end
        
        
        
        function sample = time2sample(time, fs)
        %%time2sample Takes a time in a timebase and converts to a sample
        %%point.
        %
        %   SAMPLE = time2sample(TIME, FS) assumes that a timebase with
        %   sampling rate FS Hz starts on the first sample at time 0, the
        %   second sample at time 1/FS, etc. It then converts a TIME in
        %   that timebase to a sample rate.
        
            sample = round((time * fs) + 1);
        end
    end
end
