classdef RasterDisplayFigure < handle
% RasterDisplayFigure Class for handling the plotting of raster data, along
% with other axes such as running velocity and continuous spike rates
%
%   RasterDisplayFigure Properties:
%       raster_marker_type - style of the dots on the raster
%       
%   RasterDisplayFigure Methods:
%       fill_data       - fills the axes with the raster data
%       sync_sections   - synchronize the limits of the axes
%       x_lim           - update all x-axes
%       y_lim           - update all y-axes
%
%   See also: Raster, RasterDisplaySection, TracePlot, PSTH, RasterData

    properties
        
        h_fig
        h_section = RasterDisplaySection.empty()
        h_title
        h_section_title
        h_raster = Raster.empty()
        h_motion = TracePlot.empty()
        h_rates = TracePlot.empty()
        h_psth = PSTH.empty()
        n = 0
        raster_marker_type = 'dot'
        max_spike_rate = -inf;
    end
    
    
    
    methods
        
        function obj = RasterDisplayFigure(n_prot)
        %%RasterDisplayFigure
        %
        %   RasterDisplayFigure(NUMBER_OF_SECTIONS) creates a figure with
        %   NUMBER_OF_SECTIONS sections. Each section will be controlled by
        %   the class RasterDisplaySection, and is composed of 3 axes for
        %   raster data, velocity data and firing rate data. Currently
        %   NUMBER_OF_SECTIONS can only be one of 2, 3, 4, or 6.
        
            % setup a nice figure
            obj.h_fig = figure('Visible', 'off');
            
            % get screen limitations
            set(0, 'units', 'centimeters');
            
            rt = get(0, 'screensize');
            H = min(29.7, rt(4));
            
            % set the properties
            set(obj.h_fig, 'units', 'centimeters', ...
                           'paperunits', 'centimeters', ...
                           'paperorientation', 'portrait', ...
                           'papertype', 'a4', ...
                           'papersize', [21, H], ...
                           'paperposition', [0, 0, 21, H], ...
                           'outerposition', [0, 0, 21, H], ...
                           'color', 'w');
            
            title_height = 1.5;
                       
            % work out the position of the axes
            if n_prot == 1
                
                section_size = obj.h_fig.InnerPosition(3:4) - [0, title_height];
                x_i = 1;
                y_i = 1;
                section_position = [0, 0, section_size(1), section_size(2)];
                obj.h_section(x_i, y_i) = RasterDisplaySection(obj.h_fig);
                obj.h_section(x_i, y_i).set_position(section_position);
                
            elseif n_prot == 2
                
                section_size = (obj.h_fig.InnerPosition(3:4) - [0, title_height])./[2, 3];
                c = 0;
                for y_i = 1
                    for x_i = 1 : 2
                        c = c + 1;
                        if c > 3; break; end
                        section_position = [(x_i - 1) * section_size(1), ...
                            (3 - y_i) * section_size(2), ...
                            section_size(1), ...
                            section_size(2)];
                        obj.h_section(x_i, y_i) = RasterDisplaySection(obj.h_fig);
                        obj.h_section(x_i, y_i).set_position(section_position);
                    end
                end
                
            elseif n_prot == 3
                
                section_size = (obj.h_fig.InnerPosition(3:4) - [0, title_height])./[2, 3];
                c = 0;
                for y_i = 1 : 2
                    for x_i = 1 : 2
                        c = c + 1;
                        if c > 3; break; end
                        section_position = [(x_i - 1) * section_size(1), ...
                            (3 - y_i) * section_size(2), ...
                            section_size(1), ...
                            section_size(2)];
                        obj.h_section(x_i, y_i) = RasterDisplaySection(obj.h_fig);
                        obj.h_section(x_i, y_i).set_position(section_position);
                    end
                end
            elseif n_prot == 4
                
                section_size = (obj.h_fig.InnerPosition(3:4) - [0, title_height])./[2, 3];
                c = 0;
                for y_i = 1 : 2
                    for x_i = 1 : 2
                        c = c + 1;
                        if c > 4; break; end
                        section_position = [(x_i - 1) * section_size(1), ...
                            (3 - y_i) * section_size(2), ...
                            section_size(1), ...
                            section_size(2)];
                        obj.h_section(x_i, y_i) = RasterDisplaySection(obj.h_fig);
                        obj.h_section(x_i, y_i).set_position(section_position);
                    end
                end
                
            elseif n_prot == 6
                
                section_size = (obj.h_fig.InnerPosition(3:4) - [0, title_height])./[2, 3];
                
                for x_i = 1 : 2
                    for y_i = 1 : 3
                        
                        section_position = [(x_i - 1) * section_size(1), ...
                                            (3 - y_i) * section_size(2), ...
                                            section_size(1), ...
                                            section_size(2)];
                        obj.h_section(x_i, y_i) = RasterDisplaySection(obj.h_fig);
                        obj.h_section(x_i, y_i).set_position(section_position);
                    end
                end
                
            end
        end
        
        
        
        function fill_data(obj, x_i, y_i, spike_times, velocity_traces, spike_rates, common_t, title_str)
        %%fill_data Puts all data on the axes.
        %
        %   fill_data(X, Y, SPIKE_TIMES, VELOCITY_TRACES, SPIKE_RATES, COMMON_T, TITLE_STRING)
        %   
        %   Args:
        %       X - column of the section
        %       Y - row of the section
        %       SPIKE_TIMES - cell array with each entry containing the
        %                     spike times for each row of the raster
        %       VELOCITY_TRACES - # samples x # raster rows matrix with
        %                         velocity traces to show
        %       SPIKE_RATES - # samples x # raster rows matrix with firing
        %                     rates traces to show
        %       COMMON_T - # samples x 1 vector with the time base of the
        %                  traces in VELOCITY_TRACES and SPIKE_RATES
        
            h_raster = Raster(spike_times, obj.h_section(x_i, y_i).h_ax(1), obj.raster_marker_type); %#ok<*PROPLC>
            h_motion = TracePlot(velocity_traces, common_t, obj.h_section(x_i, y_i).h_ax(2));
            
            h_raster.ylabel('Bout #');
            h_raster.xlabel('');
            h_motion.add_traces();
            h_motion.add_sd();
            h_motion.ylim([0, nan]);
            h_motion.ylabel('cm/s');
            h_motion.xlabel('');
            
            % Add firing rate lines on the same axis with right y-axis
            yyaxis(obj.h_section(x_i, y_i).h_ax(2), 'right');
            hold(obj.h_section(x_i, y_i).h_ax(2), 'on');
            
            % Calculate max and median firing rates across all trials
            max_firing_rate = max(spike_rates, [], 2);
            median_firing_rate = median(spike_rates, 2);
            
            % Apply slight smoothing (moving average with window of 5 samples)
            window_size = 5;
            max_firing_rate_smoothed = movmean(max_firing_rate, window_size, 'omitnan');
            median_firing_rate_smoothed = movmean(median_firing_rate, window_size, 'omitnan');
            
            % Plot max firing rate (thin black line)
            plot(obj.h_section(x_i, y_i).h_ax(2), common_t, max_firing_rate_smoothed, 'k-', 'LineWidth', 1);
            % Plot median firing rate (thick black line)
            plot(obj.h_section(x_i, y_i).h_ax(2), common_t, median_firing_rate_smoothed, 'k-', 'LineWidth', 2);
            
            ylabel(obj.h_section(x_i, y_i).h_ax(2), 'Hz');
            ylim(obj.h_section(x_i, y_i).h_ax(2), [0, inf]);
            
            % Switch back to left y-axis to keep velocity properties
            yyaxis(obj.h_section(x_i, y_i).h_ax(2), 'left');
            
            for i = 1 : 2
                set(obj.h_section(x_i, y_i).h_ax(i), 'xlim', common_t([1, end]));%[-0.2, 0.5])
            end
            
            title(h_raster.h_ax, title_str, 'interpreter', 'none')
            
            % Set ylabels after any axis modifications
            ylabel(obj.h_section(x_i, y_i).h_ax(2), 'cm/s', 'interpreter', 'none', 'fontsize', 8);
            
            obj.h_raster{x_i, y_i} = h_raster;
            obj.h_motion{x_i, y_i} = h_motion;
            
            obj.max_spike_rate = max([obj.max_spike_rate, max(max_firing_rate_smoothed)]);
            
        end
        
        
        
        function sync_sections(obj)
        %%sync_sections Loops through axes and makes sure axes of the same
        %%type have same y limits 
        
            m = inf*[1, 1];
            M = -inf*[1, 1];
            for x_i = 1 : size(obj.h_section, 1)
                for y_i = 1 : size(obj.h_section, 2)
                    if isempty(obj.h_section(x_i, y_i)); continue; end
                    
                    for i = 1 : length(obj.h_section(x_i, y_i).h_ax)
                        if isempty(obj.h_section(x_i, y_i).h_ax(i)); continue; end
                        yl = get(obj.h_section(x_i, y_i).h_ax(i), 'ylim');
                        m(i) = min([m(i), yl]);
                        M(i) = max([M(i), yl]);
                    end
                end
            end
            
            for x_i = 1 : size(obj.h_section, 1)
                for y_i = 1 : size(obj.h_section, 2)
                    if isempty(obj.h_section(x_i, y_i)); continue; end
                    for i = 1 : length(obj.h_section(x_i, y_i).h_ax)
                        if isempty(obj.h_section(x_i, y_i).h_ax(i)); continue; end
                        set(obj.h_section(x_i, y_i).h_ax(i), 'ylim', [m(i), M(i)]);
                    end
                end
            end
        end
        
        
        
        function x_lim(obj, val)
        %%x_lim Loop through axes of the same type and give them the same
        %%x-axis limits
        %
        %   x_lim(VALUE) where value is a 1 x 2 vector containing the
        %   limits
        
            for x_i = 1 : size(obj.h_section, 1)
                for y_i = 1 : size(obj.h_section, 2)
                    if isempty(obj.h_section(x_i, y_i)); continue; end
                    for i = 1 : length(obj.h_section(x_i, y_i).h_ax)
                        if isempty(obj.h_section(x_i, y_i).h_ax(i)); continue; end
                        set(obj.h_section(x_i, y_i).h_ax(i), 'xlim', val);
                    end
                end
            end
        end
        
        
        
        function y_lim(obj, val, idx)
        %%y_lim Loop through axes of the same type and give them the same
        %%y-axis limits
        %
        %   y_lim(VALUE, INDEX) where VALUE is a 1 x 2 vector containing the
        %   limits. INDEX is the index of axis within the section that is
        %   to be updated (e.g. if there are 3 axes in each section INDEx
        %   can be 1 2 or 3).
        
            if val(2) <= val(1)
                return
            end
            
            for x_i = 1 : size(obj.h_section, 1)
                for y_i = 1 : size(obj.h_section, 2)
                    if isempty(obj.h_section(x_i, y_i)); continue; end
                    if length(obj.h_section(x_i, y_i).h_ax) < idx; continue; end
                    set(obj.h_section(x_i, y_i).h_ax(idx), 'ylim', val);
                end
            end
        end
    end
end