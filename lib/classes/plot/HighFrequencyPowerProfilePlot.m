classdef HighFrequencyPowerProfilePlot < handle
    
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
        
        boundaries
        region_str
        anat_l5
    end
    
    
    
    methods
        
        function obj = HighFrequencyPowerProfilePlot()
        end
        
        
        
        function set_probe_track(obj, probe_track)
            obj.h_region_boundaries = RegionBoundariesPlot(probe_track);
        end
        
        
        
        function batches(obj)
        %%plots the power profile for each batch (up to 20 batches)
            
            obj.h_fig = figure('position', [150, 140, 1400, 850], ...
                               'papersize', [50, 40], ...
                               'renderer', 'painters');
            
            for ii = 1 : min(length(obj.hf_power.sample_start_t), 20)
                
                subplot(4, 5, ii);
                
                channel_ids = obj.hf_power.channel_ids_to_process;
                from_tip_um = obj.hf_power.raw_channels_from_tip;
                
                % channels to ignore
                ignore_idx = ismember(channel_ids, obj.hf_power.rec.reference_channel_ids) | ...
                             ismember(channel_ids, obj.hf_power.bad_channel_ids);
                
                power = obj.hf_power.power_raw(:, ii);
                power(ignore_idx) = nan;
                
                plot(gca, from_tip_um, power, 'color', 'k', 'linewidth', 0.75);
                
                obj.h_region_boundaries.plot_boundary_lines(gca);
                obj.h_region_boundaries.plot_middle_of_VISp5(gca)
                
                box off
                
                set(gca, 'plotboxaspectratio', [2, 1, 1]);
                
                title(sprintf('%i-%is', obj.hf_power.sample_start_t(ii), ...
                                        obj.hf_power.sample_start_t(ii) + obj.hf_power.batch_duration_s));
            end
        end
        
        
        
        function plot_by_column_interactive(obj, type)
        %%separates the raw power on each channel into separate columns and
        %%plots them all
            
            VariableDefault('type', 'processed');
        
            obj.h_fig = figure('position', [150, 140, 1400, 850], ...
                               'papersize', [50, 40], ...
                               'renderer', 'painters');
            obj.h_ax = axes(obj.h_fig, 'nextplot', 'add');
            
            obj.current_offset = 0;
            
            if strcmp(type, 'raw')
                obj.plot_all_batches_to_use_and_columns(obj.h_ax);
            else
                obj.plot_all_batches_to_use_avg_columns(obj.h_ax);
            end
            
            obj.current_offset = 0;
            obj.h_region_boundaries.plot_boundary_lines(obj.h_ax);
            obj.h_region_boundaries.plot_middle_of_VISp5(obj.h_ax);
            obj.h_region_boundaries.print_region_labels(obj.h_ax);
            
            set(obj.h_ax, 'plotboxaspectratio', [2, 1, 1]);
            set(obj.h_fig, 'keypressfcn', @(x, y)obj.interactive_key_press(x, y));
        end
        
        
        
        function plot_all_batches_to_use_and_columns(obj, h_ax)
            
            power_raw = ...
                obj.hf_power.set_unuseful_channels_to_nan(obj.hf_power.power_raw, obj.hf_power.channel_ids_to_process);

            power_raw = power_raw(:, obj.hf_power.batches_to_use);
            
            % separate the data into electrode columns
            [power_raw_columns, column_channel_ids] = ...
                obj.hf_power.separate_into_electrode_columns(power_raw, obj.hf_power.channel_ids_to_process);
                
            for ii = 1 : size(power_raw_columns, 2)
                
                from_tip_um = obj.hf_power.rec.channel_from_tip_um(column_channel_ids{ii});
                
                for jj = 1 : size(power_raw_columns{ii}, 2)
                    plot(h_ax, from_tip_um, power_raw_columns{ii}(:, jj), 'linewidth', 0.75);
                end
            end
        end
        
        
        
        function plot_all_batches_to_use_avg_columns(obj, h_ax)
            
            power_interp = obj.hf_power.power_interp(:, obj.hf_power.batches_to_use);
            power_norm = bsxfun(@rdivide, power_interp, max(power_interp, [], 1));
            from_tip_um = obj.hf_power.interp_points_from_tip();
            plot(h_ax, from_tip_um, power_norm, 'linewidth', 0.75);
            plot(h_ax, from_tip_um, mean(power_norm, 2), 'color', 'k', 'linewidth', 2)
            plot(h_ax, from_tip_um, smooth(mean(power_norm, 2), 100), 'color', [0.5, 0.5, 0.5], 'linewidth', 1)
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
            
            obj.h_region_boundaries.apply_offset(obj.current_offset);
%             set(obj.h_offset_text, 'string', sprintf('Current offset = %i', obj.current_offset));
        end
    end
end