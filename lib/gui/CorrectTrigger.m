classdef CorrectTrigger < handle
    
    properties (SetAccess = private)
        
        ctl
        probe_id
        
        original_trigger
        current_trigger
        
        app
        
        h_line
        n_lines = 0
        h_trigger
        
        remove_between = {}
        
        active = true
    end
    
    
    
    methods
        
        function obj = CorrectTrigger(ctl, probe_id)
            
            obj.ctl = ctl;
            obj.probe_id = probe_id;
            
            trigger = obj.ctl.load.trigger_mat(probe_id);
            
            obj.original_trigger = trigger;
            obj.current_trigger = trigger;
            
            obj.app = CorrectTriggerGUI(obj);
            
            obj.h_trigger = plot(obj.app.trigger_axes, 1:10:length(trigger), trigger(1:10:end));
            
            title(obj.app.trigger_axes, probe_id, 'interpreter', 'none');
            
            set(obj.app.ui_figure, 'closerequestfcn', @(x, y)obj.close_figure(x, y));
            
            fprintf('click and drag to remove triggers\n');
        end
        
        
        
        function add_line(obj)
            
            figure(obj.app.ui_figure);
            
            % create the line
            obj.n_lines = obj.n_lines + 1;
            obj.h_line{obj.n_lines} = drawline(obj.app.trigger_axes);
            
            % max and min of line
            low_val = round(min(obj.h_line{obj.n_lines}.Position(:, 1)));
            high_val = round(max(obj.h_line{obj.n_lines}.Position(:, 1)));
            
            % don't allow single points (e.g. if user has clicked on Add more than
            % once before placing the line)
            if isequal(single(low_val), single(high_val))
                   obj.remove_last_added_line();
                   return
            end
            
            % make sure at least one point exists in the data
            if high_val < 1 || low_val > length(obj.current_trigger)
                   obj.remove_last_added_line();
                   return
            end
            
            % add a callback for deletion
            n = obj.n_lines;
            obj.h_line{obj.n_lines}.UserData.OnCleanup = onCleanup(@(x)obj.delete_line(n));
            
            % store the indices between which we want to remove
            obj.remove_between{obj.n_lines} = [max(low_val, 1), min(high_val, length(obj.current_trigger))];
            
            % update the line
            obj.h_line{obj.n_lines}.Position = ...
                [obj.remove_between{obj.n_lines}(:), mean(obj.h_line{obj.n_lines}.Position(:, 2))*[1;1]];
            
            % do the removal and update
            obj.update_trigger();
        end
        
        
        
        function remove_last_added_line(obj)
            
            delete(obj.h_line{obj.n_lines});
            obj.h_line(obj.n_lines) = [];
            obj.n_lines = obj.n_lines - 1;
        end
        
        
        
        function update_trigger(obj)
            
            obj.current_trigger = obj.original_trigger;
            
            for ii = 1 : length(obj.h_line)
                if ~isempty(obj.remove_between{ii})
                    obj.remove_from_trigger(obj.remove_between{ii}(1), obj.remove_between{ii}(2));
                end
            end
        end
        
        
        
        function remove_from_trigger(obj, start_idx, end_idx)
            
            % remove from trigger, set to baseline value
            obj.current_trigger(start_idx:end_idx) = obj.ctl.get_trigger_baseline_val(obj.probe_id);
            
            % update display
            set(obj.h_trigger, 'ydata', obj.current_trigger(1:10:end));
        end
        
        
        
        function delete_line(obj, n)
            
            obj.h_line{n} = [];
            obj.remove_between{n} = [];
            obj.update_trigger();
        end
        
        
        
        function save(obj)
            
            if isequal(obj.original_trigger, obj.current_trigger)
                warning('original and new triggers are the same');
                return
            end
            
            % first save original trigger
            obj.ctl.save.original_trigger_mat(obj.probe_id, obj.original_trigger);
            % save new trigger
            obj.ctl.save.trigger_mat(obj.probe_id, obj.current_trigger);
            % save details of points removed
            points = [];
            for ii = 1 : length(obj.remove_between)
                if ~isempty(obj.remove_between{ii})
                    points(end+1, :) = obj.remove_between{ii};
                end
            end
            obj.ctl.save.trigger_points_removed(obj.probe_id, points);
        end
        
        
        
        function close_figure(obj, ~, ~)
            
            obj.active = false;
            delete(obj.app.ui_figure);
        end
    end
end
