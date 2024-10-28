classdef AccelerationBins < handle
    % AccelerationBins Class for computing acceleration bins from a set of trials
    %
    %   AccelerationBins Properties:
    %       trials          - a cell array of Trial objects used to compute bins
    %       prc_per_bin     - percentage of data in each bin (default: 5%)
    %       bin_edges       - the acceleration bin edges
    %       mode            - mode of analysis, either 'acc' for acceleration 
    %                         or 'dec' for deceleration
    %       n_bins          - number of bins (dependent property)
    %       bin_centers     - the center value of each bin (dependent property)
    %
    %   AccelerationBins Methods:
    %       acceleration_bounds   - computes the acceleration bin edges based 
    %                               on the mode and percentage per bin.
    %
    %   See also: ShuffleTuning, TuningTable, AccelerationTuningCurve
    %
    %   TODO: Ensure `bin_edges` recomputation when `prc_per_bin` changes.
    
        properties 
            trials           % Cell array of Trial objects
            prc_per_bin      % Percentage of data in each bin (default: 5%)
            bin_edges        % Acceleration bin edges based on trials data
            mode             % Mode for binning: 'acc' for acceleration, 'dec' for deceleration
        end
        
        properties (Dependent = true)
            n_bins           % Number of bins determined by prc_per_bin
            bin_centers      % Centers of the bins calculated from bin_edges
        end
        
    
        methods
            
            function obj = AccelerationBins(trials, mode)
            % AccelerationBins Constructor for AccelerationBins class
            %
            %   AccelerationBins(TRIALS, MODE) initializes the object with a
            %   cell array of Trial objects (TRIALS) and a specified binning mode (MODE).
            %   Mode can be 'acc' (acceleration) or 'dec' (deceleration).
                obj.mode = mode;
                obj.trials = trials;
                obj.prc_per_bin = 5;
                obj.bin_edges = obj.acceleration_bounds();
            end
            
    
            function set.prc_per_bin(obj, val)
            % SET.PRC_PER_BIN Adjust the percentage of data per bin
            %
            %   SET.PRC_PER_BIN(VAL) sets the property prc_per_bin to the value
            %   needed for the closest whole number of bins.
                n = ceil(100/val);
                obj.prc_per_bin = 100/n;
            end
            
            
            function val = get.n_bins(obj)
            % GET.N_BINS Compute the number of bins
            %
            %   Returns the number of bins based on the current prc_per_bin value.
                val = round(100/obj.prc_per_bin);
            end
            
            
            function val = get.bin_centers(obj)
            % GET.BIN_CENTERS Calculate the centers of the bins
            %
            %   Returns a vector of bin center values, derived from bin_edges.
                val = (obj.bin_edges(1:end-1) + obj.bin_edges(2:end))/2;
            end
            
            
            function bounds = acceleration_bounds(obj)
            % ACCELERATION_BOUNDS Compute the acceleration bin edges from trial data
            %
            %   BOUNDS = ACCELERATION_BOUNDS() calculates the acceleration bin 
            %   edges based on data in TRIALS. Only acceleration values that match 
            %   the selected mode ('acc' or 'dec') are included.
            %
            %   Output:
            %       bounds - a (#bins+1)x1 vector containing the bin edges.
                all_accelerations = [];
                
                for trial_i = 1 : length(obj.trials)
                    
                    acc = obj.trials{trial_i}.acceleration();
                    idx_motion = obj.trials{trial_i}.motion_mask();
                    selected = acc(idx_motion);
    
                    if obj.mode == "acc"
                        selected = selected(selected > 0);
                    elseif obj.mode == "dec"
                        selected = selected(selected < 0);
                    end
                    all_accelerations = [all_accelerations; selected];
                end
                
                prcbnd = 0 : obj.prc_per_bin : 100+eps;
                bounds = prctile(all_accelerations, prcbnd);
            end
        end
    end
    