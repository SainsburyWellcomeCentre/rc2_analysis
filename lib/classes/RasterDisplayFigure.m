classdef RasterDisplayFigure < handle
    
    properties
        h_fig
        h_section = RasterDisplaySection.empty()
        h_title
        h_section_title
        h_raster = Raster.empty()
        h_motion = TracePlot.empty()
        h_rates = TracePlot.empty()
        n = 0
        raster_marker_type = 'dot'
    end
    
    methods
        
        function obj = RasterDisplayFigure(n_prot)
            
            % setup a nice figure
            obj.h_fig = figure();
            
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
            
            title_height = 2;
                       
            % work out the position of the axes
            if n_prot == 1
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
            
            h_raster = Raster(spike_times, obj.h_section(x_i, y_i).h_ax(1), obj.raster_marker_type); %#ok<*PROPLC>
            h_motion = TracePlot(velocity_traces, common_t, obj.h_section(x_i, y_i).h_ax(2));
            h_rates = TracePlot(spike_rates, common_t, obj.h_section(x_i, y_i).h_ax(3));
%             h_psth = PSTH(spike_times, obj.h_section(x_i, y_i).h_ax(3));
            
            h_rates.mean_colour('k');
            h_raster.ylabel('Bout #');
            h_raster.xlabel('');
            h_motion.add_traces();
            h_motion.add_sd();
            h_motion.ylim([0, nan]);
            h_motion.ylabel('cm/s');
            h_motion.xlabel('');
            h_rates.mean_colour('k');
            h_rates.add_sem([0.6, 0.6, 0.6]);
            h_rates.ylabel('Hz');
            
            for i = 1 : 3
                set(obj.h_section(x_i, y_i).h_ax(i), 'xlim', common_t([1, end]));%[-0.2, 0.5])
            end
            
            title(h_raster.h_ax, title_str)
            
            obj.h_raster{x_i, y_i} = h_raster;
            obj.h_motion{x_i, y_i} = h_motion;
            obj.h_rates{x_i, y_i} = h_rates;
            
        end
        
        
        function sync_sections(obj)
            
            m = inf*[1, 1, 1];
            M = -inf*[1, 1, 1];
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