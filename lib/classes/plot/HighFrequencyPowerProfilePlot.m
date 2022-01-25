classdef HighFrequencyPowerProfilePlot < handle
% HighFrequencyPowerProfilePlot Class for plotting the high-frequency power
% profile across the length of a probe shank.
%
%  HighFrequencyPowerProfilePlot Properties:
%
%        hf_power               - object of type HighFrequencyPowerProfile containing data to plot
%        clusters_from_tip_um   - list of depths of clusters from the tip of the probe
%        current_offset         - amount to offset anatomical boundaries
%
%        h_fig                  - handle to the most recent figure plotted
%        h_ax                   - handle to the most recent axis plotted
%        probe_track            - object of type ProbeTrack containing the information about anatomical regions
%        h_region_boundaries    - object of type RegionBoundariesPlot
%        h_offset_text          - handle to text on the plot regarding offset
%
%  HighFrequencyPowerProfilePlot Methods:
%
%       set_probe_track             - sets up the figure to use anatomical information in a ProbeTrack object
%       raw_batches                 - plots the power across all channels separately for each batch in a separate subplot
%       interpolated_batches        - similar to `raw_batches` but plots the power interpolated across 
%                                     banks of electrodes instead of every channel
%       overlay_raw_batches         - plots power across all channels for each batch on the same plot
%       plot_by_column_interactive  - creates interactive plot to view the averaged power profile with anatomical boundaries overlaid
%       plot_legacy                 - legacy function to create a plot in an old style (remove)
%       plot_summary                - create figure with 3 axes showing various features of the power profile
%
%     Internal:
%       plot_all_batches_to_use_and_columns - plot of power for batches to use, each electrode column is a separate line
%       overlay_cluster_histogram   - adds histogram of cluster depths to an axis
%       interactive_key_press       - internal function for interactive  shift 
%
%   See also: HighFrequencyPowerProfile, ProbeTrack, RegionBoundariesPlot

    properties
        
        hf_power
        clusters_from_tip_um
        current_offset
    end
    
    properties (SetAccess = private)
        
        h_fig
        h_ax
        
        probe_track
        
        h_region_boundaries
        h_offset_text
    end
    
    
    
    methods
        
        function obj = HighFrequencyPowerProfilePlot(hf_power)
        %%HighFrequencyPowerProfilePlot
        %
        %   HighFrequencyPowerProfilePlot(HF_POWER) creates the object
        %   where HF_POWER is an object of class HighFrequencyPowerProfile
        %   and contains data about the high-frequency power profile along
        %   a probe shank.
        %
        %   The HF_POWER object must be fully formed (i.e. all analysis has
        %   been run and the properties are complete).
        
            obj.hf_power = hf_power;
            obj.set_probe_track(hf_power.probe_track);
        end
        
        
        
        function set_probe_track(obj, probe_track)
        %%set_probe_track Sets the anatomical information to use
        %
        %   set_probe_track(PROBE_TRACK) takes in the anatomical
        %   information in PROBE_TRACK, which should be an object of type
        %   ProbeTrack and contain information about the anatomical regions
        %   in which a probe shank is located.
        %
        %   See also: ProbeTrack
        
            if isempty(probe_track); return; end
            obj.h_region_boundaries = RegionBoundariesPlot(probe_track);
        end
        
        
        
        function raw_batches(obj)
        %%raw_batches Plots the power profile for each batch (up to 20 batches)
        %
        %   raw_batches() plots the high-frequency power profile across the
        %   probe. Each batch is plotted on a separate axis, and for each
        %   batch the power on every channel is plotted. This is useful to
        %   selecting batches with clear cortical peaks.
            
            obj.h_fig = figure('position', [150, 140, 1400, 850], ...
                               'papersize', [50, 40], ...
                               'renderer', 'painters');
            
            for ii = 1 : min(length(obj.hf_power.sample_start_t), 20)
                
                subplot(4, 5, ii);
                
                channel_ids = obj.hf_power.channel_ids_to_process;
                from_tip_um = obj.hf_power.raw_channels_from_tip;
                
                % channels to ignore
                ignore_idx = ismember(channel_ids, obj.hf_power.rec.reference_channel_ids) | ...
                             ismember(channel_ids, obj.hf_power.bad_channel_ids) | ...
                             ismember(channel_ids, obj.hf_power.channel_ids_to_remove);
                
                power = obj.hf_power.power_raw(:, ii);
                power(ignore_idx) = nan;
                
                plot(gca, from_tip_um, power, 'color', 'k', 'linewidth', 0.75);
                
                if ~isempty(obj.h_region_boundaries)
                    obj.h_region_boundaries.plot_boundary_lines(gca);
                    obj.h_region_boundaries.plot_middle_of_VISp5(gca)
                end
                
                box off
                
                set(gca, 'plotboxaspectratio', [2, 1, 1]);
                
                title(sprintf('%i-%is', obj.hf_power.sample_start_t(ii), ...
                                        obj.hf_power.sample_start_t(ii) + obj.hf_power.batch_duration_s));
            end
        end
        
        
        
        function interpolated_batches(obj)
        %%interpolated_batches Plots the power profile for each batch (up to 20 batches)
        %
        %   interpolated_batches() plots the high-frequency power profile across the
        %   probe. Each batch is plotted on a separate axis, and for each
        %   batch the interpolated power along each electrode columns is
        %   shown. Similar to `raw_batches`.
        
            obj.h_fig = figure('position', [150, 140, 1400, 850], ...
                               'papersize', [50, 40], ...
                               'renderer', 'painters');
            
            for ii = 1 : min(length(obj.hf_power.sample_start_t), 20)
                
                subplot(4, 5, ii);
                
                from_tip_um = obj.hf_power.interp_points_from_tip();
                
                plot(gca, from_tip_um, obj.hf_power.power_interp(:, :, ii), 'linewidth', 0.75);
                
                if ~isempty(obj.h_region_boundaries)
                    obj.h_region_boundaries.plot_boundary_lines(gca);
                    obj.h_region_boundaries.plot_middle_of_VISp5(gca)
                end
                
                box off
                
                set(gca, 'plotboxaspectratio', [2, 1, 1]);
                
                title(sprintf('%i-%is', obj.hf_power.sample_start_t(ii), ...
                                        obj.hf_power.sample_start_t(ii) + obj.hf_power.batch_duration_s));
            end
        end
        
        
        
        function overlay_raw_batches(obj)
        %%overlay_raw_batches Plots the power profile for each batch (up to 20 batches)
        %
        %   overlay_raw_batches() plots the high-frequency power profile across the
        %   probe. Each batch is plotted on the same axis and the power on
        %   all channels is shown for each batch.
        
            obj.h_fig = figure('position', [150, 140, 1400, 850], ...
                'papersize', [50, 40], ...
                'renderer', 'painters');
            
            power_raw = ...
                obj.hf_power.set_unuseful_channels_to_nan(obj.hf_power.power_raw, obj.hf_power.channel_ids_to_process);
            
            from_tip_um = obj.hf_power.rec.channel_from_tip_um(obj.hf_power.channel_ids_to_process);
            
            plot(gca, from_tip_um, power_raw, 'linewidth', 0.75);
            
            if ~isempty(obj.h_region_boundaries)
                obj.h_region_boundaries.plot_boundary_lines(gca);
                obj.h_region_boundaries.plot_middle_of_VISp5(gca)
                obj.h_region_boundaries.print_region_labels(gca);
            end
            
            box off
            
            set(gca, 'plotboxaspectratio', [2, 1, 1]);
            
        end
        
        
        
        function plot_by_column_interactive(obj, type)
        %%plot_by_column_interactive Create interactive plot to view the averaged power profile with anatomical boundaries overlaid
        %
        %   plot_by_column_interactive() creates an interactive plot
        %   showing the interpolated power (along electrode columns) along the probe shank, 
        %   for all the batches selected (selected batches are in the `batches_to_use` property of `hf_power`). 
        %   The average power profile is also plotted along with an indication of the peak HF power in cortex.
        %
        %   Anatomical boundaries from `probe_track` are also initially overlaid as if there
        %   were no offset. The user can then shift those boundaries using
        %   the left and right arrow keys (1um shift) or `a` and `d` keys
        %   on the keyboard (10um shift), to check whether aligning mid L5
        %   from anatomy and HF power cortical peak makes sense or if extra
        %   criteria are required to detect the peak.
            
            VariableDefault('type', 'raw');
        
            obj.h_fig = figure('position', [150, 140, 1400, 850], ...
                               'papersize', [50, 40], ...
                               'renderer', 'painters');
                           
            obj.h_ax = axes(obj.h_fig, 'nextplot', 'add');
            
            if strcmp(type, 'raw')
                obj.plot_all_batches_to_use_and_columns(obj.h_ax);
            end
            
            from_tip_um = obj.hf_power.interp_points_from_tip();
            plot(obj.h_ax, from_tip_um, mean(mean(obj.hf_power.power_smooth(:, :, obj.hf_power.batches_to_use), 2), 3), 'color', 'b', 'linewidth', 2);
            
            obj.current_offset = 0;
            
            set(obj.h_ax, 'ylim', [0, max(get(obj.h_ax, 'ylim'))]);
            
            obj.h_offset_text = text(max(get(obj.h_ax, 'xlim')), ...
                                     max(get(obj.h_ax, 'ylim')), ...
                                     sprintf('Current offset: %ium', obj.current_offset), ...
                                     'horizontalalignment', 'right', ...
                                     'verticalalignment', 'top');
            
            if ~isempty(obj.h_region_boundaries)
                obj.h_region_boundaries.plot_boundary_lines(obj.h_ax);
                obj.h_region_boundaries.plot_middle_of_VISp5(obj.h_ax);
                obj.h_region_boundaries.print_region_labels(obj.h_ax);
            end
            
            plot(obj.h_ax, obj.hf_power.ephys_l5*[1, 1], get(obj.h_ax, 'ylim'), 'color', 'b');
            
            obj.overlay_cluster_histogram(obj.h_ax);
            
            title_x = mean(get(obj.h_ax, 'xlim'));
            title_y = max(get(obj.h_ax, 'ylim')) + 0.1*range(get(obj.h_ax, 'ylim'));
            
            text(obj.h_ax, title_x, title_y, sprintf('%s, shank %i', obj.hf_power.probe_id, obj.hf_power.shank_id), ...
                  'interpreter', 'none', ...
                  'fontsize', 14, ...
                  'fontweight', 'bold', ...
                  'horizontalalignment', 'center');
            
            set(obj.h_ax, 'plotboxaspectratio', [2, 1, 1]);
            set(obj.h_fig, 'keypressfcn', @(x, y)obj.interactive_key_press(x, y));
        end
        
        
        
        function plot_legacy(obj, hf_power, probe_track, clusters_from_tip)
            
            % gather info needed to make plot
            ephys_l5 = hf_power.search_for_l5();
            
            power_profile = hf_power.power_smooth;
            channels_from_probe_tip = hf_power.interp_points_from_tip();
            
            norm_power_profiles = bsxfun(@rdivide, power_profile, max(power_profile, [], 1));
            norm_power = norm_power_profiles/max(norm_power_profiles);
            
            if ~isempty(probe_track)
                [boundaries, ~, region_str] = probe_track.region_boundaries();
                anatomy_l5 = probe_track.mid_l5_visp();
                probe_track.offset = ephys_l5 - anatomy_l5;
                boundaries_adjusted = probe_track.region_boundaries_adjusted();
            else
                boundaries = 0;
                boundaries_adjusted = boundaries;
            end
            
            % one figure
            obj.h_fig = figure('position', [680, 84, 1100, 894]);
            
            % left most plot axis
            subplot(1, 3, 1); hold on;
            
            % draw the HF power profile line
            plot(norm_power_profiles, channels_from_probe_tip, 'color', 'k', 'linewidth', 2);
            
            % some formatting
            set(gca, 'ylim', [0, 2500], 'xlim', [0, 1.2])
            
            xlabel(sprintf('Norm. power (%i-%iHz)', hf_power.lower_hz, hf_power.upper_hz));
            ylabel('Distance from probe tip (\mum)')
            
            set(gca, 'plotboxaspectratio', [1, 3, 1])
            box off;
            title('HF power');
            
            % plot line for middle of L5 from ephys
            line(get(gca, 'xlim'), ephys_l5*[1, 1], 'color', 'r', 'linewidth', 2);
            % text for distance of line from tip
            text(min(get(gca, 'xlim')), ephys_l5, sprintf('peak power: %i', ephys_l5), ...
                'verticalalignment', 'bottom', 'horizontalalignment', 'left', 'color', 'r');
            
            if ~isempty(probe_track)
                % plot line for middle of L5 from anatomy
                line(get(gca, 'xlim'), anatomy_l5*[1, 1], 'color', 'b', 'linewidth', 2);
                % text for distance of line from tip
                text(min(get(gca, 'xlim')), anatomy_l5, sprintf('mid VISp5: %i', round(anatomy_l5)), ...
                    'verticalalignment', 'top', 'horizontalalignment', 'left', 'color', 'b');
            end
            
            
            % plot the boundaries
            for i = 1 : length(boundaries)
                % boundary
                line(get(gca, 'xlim'), boundaries(i)*[1, 1], 'color', 'k', 'linestyle', '--');
                % text in middle of boundary
                if i < length(boundaries)
                    y = sum(boundaries([i, i+1]))/2;
                    text(max(get(gca, 'xlim')), y, region_str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
                end
            end
            
            % middle plot axis
            subplot(1, 3, 2); hold on;
            
            % again plot HF power profile line
            plot(norm_power, channels_from_probe_tip, 'color', 'k', 'linewidth', 2); %p(1:2:end, :),
            
            % formatting
            set(gca, 'ylim', [0, 2500], 'xlim', [0, 1.2])
            ylabel('Distance from probe tip (\mum)')
            xlabel(sprintf('Norm. power (%i-%iHz)', hf_power.lower_hz, hf_power.upper_hz));
            set(gca, 'plotboxaspectratio', [1, 3, 1])
            box off
            title('HF power (boundaries corrected)');
            
            % plot line for middle of L5 from ephys (also will be equal to mid L5 from
            % corrected boundaries)
            line(get(gca, 'xlim'), ephys_l5*[1, 1], 'color', 'r', 'linewidth', 2);
            
            % plot the corrected boundaries
            for i = 1 : length(boundaries_adjusted)
                line(get(gca, 'xlim'), boundaries_adjusted(i)*[1, 1], 'color', 'k', 'linestyle', '--');
                if i < length(boundaries_adjusted)
                    y = sum(boundaries_adjusted([i, i+1]))/2;
                    text(max(get(gca, 'xlim')), y, region_str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
                end
            end
            
            % right most plot axis
            subplot(1, 3, 3); hold on;
            
            % histogram of number of clusters found at different depths
            histogram(clusters_from_tip, 40, 'orientation', 'horizontal')
            
            % some formatting
            set(gca, 'plotboxaspectratio', [1, 3, 1])
            set(gca, 'ylim', [0, 2500])
            xlabel('# units')
            box off
            title('MUA (boundaries corrected)')
            
            % plot line for middle of L5 from ephys
            line(get(gca, 'xlim'), ephys_l5*[1, 1], 'color', 'r', 'linewidth', 2);
            
            % plot the corrected boundari
            for i = 1 : length(boundaries_adjusted)
                line(get(gca, 'xlim'), boundaries_adjusted(i)*[1, 1], 'color', 'k', 'linestyle', '--');
                if i < length(boundaries_adjusted)
                    y = sum(boundaries_adjusted([i, i+1]))/2;
                    text(max(get(gca, 'xlim')), y, region_str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
                end
            end
            
            % for saving set the paper orientation to landscape
            set(gcf, 'paperorientation', 'landscape');
        end
        
        
        
        function plot_summary(obj)
        %%plot_summary Create summary figure of the high-frequency power profile data
        %
        %   plot_summary() creates a figure with 3 axis. Left shows the
        %   high-frequency power profile from individual batches (normalised to
        %   peak for each batch), along with the unshifted boundaries from
        %   the anatomy. Peak HF power is shown (red) and mid L5 from the
        %   anatomy is shown in blue.
        %
        %   Middle axis shows the average high-frequency power profile
        %   along with the anatomical boundaries shifted so that the peak
        %   HF power and mid L5 match.
        %
        %   Right axis shows the distribution of clusters along the probe 
        %   along with the anatomical boundaries shifted so that the peak
        %   HF power and mid L5 match.
        
            % gather info needed to make plot
            ephys_l5 = obj.hf_power.ephys_l5;
            anatomy_l5 = obj.hf_power.anatomy_l5;
            
            power_profile = squeeze(mean(obj.hf_power.power_smooth(:, :, obj.hf_power.batches_to_use), 2));
            channels_from_probe_tip = obj.hf_power.interp_points_from_tip();
            
            norm_power_profiles = bsxfun(@rdivide, power_profile, max(power_profile, [], 1));
            norm_power = norm_power_profiles/max(norm_power_profiles);
            
            if ~isempty(obj.hf_power.probe_track)
                [boundaries, ~, region_str] = obj.hf_power.probe_track.region_boundaries();
                obj.hf_power.probe_track.offset = obj.hf_power.delta_l5;
                boundaries_adjusted = obj.hf_power.probe_track.region_boundaries_adjusted();
            else
                boundaries = 0;
                boundaries_adjusted = boundaries;
            end
            
            % one figure
            obj.h_fig = figure('position', [680, 84, 1100, 894]);
            
            % left most plot axis
            subplot(1, 3, 1); hold on;
            
            % draw the HF power profile line
            plot(norm_power_profiles, channels_from_probe_tip, 'color', 'k', 'linewidth', 2);
            
            % some formatting
            set(gca, 'ylim', [0, 2500], 'xlim', [0, 1.2])
            
            xlabel(sprintf('Norm. power (%i-%iHz)', obj.hf_power.lower_hz, obj.hf_power.upper_hz));
            ylabel('Distance from probe tip (\mum)')
            
            set(gca, 'plotboxaspectratio', [1, 3, 1])
            box off;
            title('HF power');
            
            % plot line for middle of L5 from ephys
            line(get(gca, 'xlim'), ephys_l5*[1, 1], 'color', 'r', 'linewidth', 2);
            % text for distance of line from tip
            text(min(get(gca, 'xlim')), ephys_l5, sprintf('peak power: %i', ephys_l5), ...
                'verticalalignment', 'bottom', 'horizontalalignment', 'left', 'color', 'r');
            
            if ~isempty(obj.hf_power.probe_track)
                % plot line for middle of L5 from anatomy
                line(get(gca, 'xlim'), anatomy_l5*[1, 1], 'color', 'b', 'linewidth', 2);
                % text for distance of line from tip
                text(min(get(gca, 'xlim')), anatomy_l5, sprintf('mid VISp5: %i', round(anatomy_l5)), ...
                    'verticalalignment', 'top', 'horizontalalignment', 'left', 'color', 'b');
            end
            
            
            % plot the boundaries
            for i = 1 : length(boundaries)
                % boundary
                line(get(gca, 'xlim'), boundaries(i)*[1, 1], 'color', 'k', 'linestyle', '--');
                % text in middle of boundary
                if i < length(boundaries)
                    y = sum(boundaries([i, i+1]))/2;
                    text(max(get(gca, 'xlim')), y, region_str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
                end
            end
            
            % middle plot axis
            subplot(1, 3, 2); hold on;
            
            % again plot HF power profile line
            plot(norm_power, channels_from_probe_tip, 'color', 'k', 'linewidth', 2); %p(1:2:end, :),
            
            % formatting
            set(gca, 'ylim', [0, 2500], 'xlim', [0, 1.2])
            ylabel('Distance from probe tip (\mum)')
            xlabel(sprintf('Norm. power (%i-%iHz)', obj.hf_power.lower_hz, obj.hf_power.upper_hz));
            set(gca, 'plotboxaspectratio', [1, 3, 1])
            box off
            title('HF power (boundaries corrected)');
            
            % plot line for middle of L5 from ephys (also will be equal to mid L5 from
            % corrected boundaries)
            line(get(gca, 'xlim'), ephys_l5*[1, 1], 'color', 'r', 'linewidth', 2);
            
            % plot the corrected boundaries
            for i = 1 : length(boundaries_adjusted)
                line(get(gca, 'xlim'), boundaries_adjusted(i)*[1, 1], 'color', 'k', 'linestyle', '--');
                if i < length(boundaries_adjusted)
                    y = sum(boundaries_adjusted([i, i+1]))/2;
                    text(max(get(gca, 'xlim')), y, region_str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
                end
            end
            
            % right most plot axis
            subplot(1, 3, 3); hold on;
            
            % histogram of number of clusters found at different depths
            histogram(obj.hf_power.clusters_from_tip_um, 40, 'orientation', 'horizontal')
            
            % some formatting
            set(gca, 'plotboxaspectratio', [1, 3, 1])
            set(gca, 'ylim', [0, 2500])
            xlabel('# units')
            box off
            title('MUA (boundaries corrected)')
            
            % plot line for middle of L5 from ephys
            line(get(gca, 'xlim'), ephys_l5*[1, 1], 'color', 'r', 'linewidth', 2);
            
            % plot the corrected boundari
            for i = 1 : length(boundaries_adjusted)
                line(get(gca, 'xlim'), boundaries_adjusted(i)*[1, 1], 'color', 'k', 'linestyle', '--');
                if i < length(boundaries_adjusted)
                    y = sum(boundaries_adjusted([i, i+1]))/2;
                    text(max(get(gca, 'xlim')), y, region_str{i}, 'verticalalignment', 'middle', 'horizontalalignment', 'right');
                end
            end
            
            FigureTitle(gcf, sprintf('%s, shank %i', obj.hf_power.probe_id, obj.hf_power.shank_id));
            
            % for saving set the paper orientation to landscape
            set(gcf, 'paperorientation', 'landscape');
        end
    end
    
    
    methods (Access = protected)
        
        function plot_all_batches_to_use_and_columns(obj, h_ax)
        %%plot_all_batches_to_use_and_columns    
            power_raw = ...
                obj.hf_power.set_unuseful_channels_to_nan(obj.hf_power.power_raw, obj.hf_power.channel_ids_to_process);

            power_raw = power_raw(:, obj.hf_power.batches_to_use);
            
            % separate the data into electrode columns
            [power_raw_columns, column_channel_ids] = ...
                obj.hf_power.separate_into_electrode_columns(power_raw, obj.hf_power.channel_ids_to_process);
            
            cols = lines(size(power_raw_columns, 2));
            
            for ii = 1 : size(power_raw_columns, 2)
                
                from_tip_um = obj.hf_power.rec.channel_from_tip_um(column_channel_ids{ii});
                
                for jj = 1 : size(power_raw_columns{ii}, 2)
                    plot(h_ax, from_tip_um, power_raw_columns{ii}(:, jj), 'color', cols(ii, :), 'linewidth', 0.75);
                end
            end
        end
        
        
        
        function overlay_cluster_histogram(obj, h_ax)
            
            [n, edges] = histcounts(obj.hf_power.clusters_from_tip_um, 0:25:2000);
            if sum(n) == 0
                return
            end
            h = histogram(h_ax, 'bincounts', max(get(h_ax, 'ylim'))*n/max(n)/2, 'binedges', edges);
            set(h, 'facealpha', 0.6);
        end
        
        
        
        function interactive_key_press(obj, ~, key_data)
        %%moves the boundaries upon keypress
            
            if strcmp(key_data.Key, 'leftarrow')
                obj.current_offset = obj.current_offset - 1;
            elseif strcmp(key_data.Key, 'rightarrow')
                obj.current_offset = obj.current_offset + 1;
            elseif strcmp(key_data.Key, 'a')
                obj.current_offset = obj.current_offset - 10;
            elseif strcmp(key_data.Key, 'd')
                obj.current_offset = obj.current_offset + 10;
            end
            
            if ~isempty(obj.h_region_boundaries)
                obj.h_region_boundaries.apply_offset(obj.current_offset);
            end
            set(obj.h_offset_text, 'string', sprintf('Current offset = %i', obj.current_offset));
        end
    end
end
