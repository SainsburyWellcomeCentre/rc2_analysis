classdef FiringRate < handle
    
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
        
            obj.spike_times = spike_times;
        end
        
        
        
        function r = get_convolution(obj, T)
            
            fs = 1/(T(2) - T(1));
            spike_times = obj.restrict_times([T(1) - obj.prepad, T(end) + obj.postpad]); %#ok<*PROPLC>
            
            [spike_train, sr_t] = obj.create_spike_train(spike_times, [T(1) - obj.prepad, T(end) + obj.postpad], fs);
            gauss = obj.gauss_filter(obj.width, obj.length, fs);
            spike_train = conv(fs*spike_train, gauss, 'same');
            r = interp1(sr_t, spike_train, T);
        end
        
        
        
        function r = get_count(obj, T)
            
            spike_times = obj.spike_times(obj.spike_times > T(1) & obj.spike_times < T(end));
            spike_points = round(length(T)*(spike_times - T(1)) / (T(end) - T(1))); %#ok<*CPROPLC>
            r = arrayfun(@(x)(sum(spike_points == x)), 1:length(T));
        end
        
        
        
        function [fr, dt, n_spikes] = get_fr_in_window(obj, T)
            
            assert(length(T) == 2);
            assert(T(2) > T(1));
            
            n_spikes = sum(obj.spike_times >= T(1) & obj.spike_times < T(2));
            dt = T(2)-T(1);
            
            fr = n_spikes / dt;
        end
        
        
        
        function [fr, dt, n_spikes] = get_fr_in_multiple_windows(obj, T)
            
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
            
            st = obj.spike_times(obj.spike_times > T(1) & obj.spike_times < T(end));
        end
        
        
        
        function [spike_train, t] = create_spike_train(obj, spike_times, tlim, fs)
            %%[spike_train, t] = CREATE_SPIKE_TRAIN(spike_times, tlim, fs)
            %
            %   Takes a list of spike times ('spike_times') and a 1x2 vector of start
            %   and end times ('tlim'), and computes a (binary) spike train, sampled
            %   at 'fs' Hz, with 1s in the spike location and 0s elsewhere.
            %
            %   If spike times occur within the resultion of 'fs', they *will not be
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
            
            sigma = sigma*fs;
            sz = sz*fs;    % length of gaussFilter vector
            x = linspace(-sz / 2, sz / 2, sz);
            b = exp(-x .^ 2 / (2 * sigma ^ 2));
            b = b/sum(b); % normalize
        end
        
        
        
        function time = sample2time(sample, fs)
            
            time = (double(sample) - 1) * (1/fs);
        end
        
        
        
        function sample = time2sample(time, fs)
            
            sample = round((time * fs) + 1);
        end
    end
end
