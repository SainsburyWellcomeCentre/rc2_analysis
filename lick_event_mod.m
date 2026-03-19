% Plot lick events per trial for each probe, Long vs Short (with raw traces)
%
% Outputs (per probe):
%   (A) Raw lick traces - Long  (up to 20 trials), exported as 2 pages (10/page)
%   (B) Raw lick traces - Short (up to 20 trials), exported as 2 pages (10/page)
%
% Output (all probes combined):
%   One combined figure (rows = probes) with 4 columns:
%     (1) Lick Events vs. Time (Long=Blue, Short=Red)     xlim: 0–100 s
%     (2) Lick Events vs. Position (Long only, Blue)      xlim: 0–120 cm, ticks every 10
%     (3) Lick Events vs. Position (Short only, Red)      xlim: 60–120 cm, ticks every 10
%     (4) Lick Events vs. Time-to-Goal (Long=Blue, Short=Red) xlim: 0–100 %
%
% KEY FIX (x-axis numbers missing):
%   - Do NOT manipulate XTickLabel at all (setting it to [] makes it "manual"
%     and it often stays blank in exported/joined PDFs).
%   - Instead, hide/show the entire x-axis using ax.XAxis.Visible = 'off'/'on'.
%     Bottom occupied row: XAxis.Visible='on' (numbers appear).
%     Other rows: XAxis.Visible='off' (clean).
%
% Export fixes (raw trace pages):
%   - Paginate raw traces (10 per page) so A4 export is readable
%   - Force tiledlayout to fill the printable page (OuterPosition)
%   - Avoid yline(Label=...) (blows up in PDF). Use small text() instead
%   - PaperPositionMode='auto' for consistent export geometry
%
% Assumptions:
%   - Trial.lick is voltage signal
%   - Trial.probe_t is time vector (same length as lick)
%   - Trial.motion_mask() exists
%   - Trial.position() exists (returns position vector aligned to probe_t)

%%
close all;

% ---------------------------
% Configuration
% ---------------------------
experiment_groups   = {'ambient_light'};
trial_group_labels  = {'RT'};

show_figs           = true;
save_figs           = true;
overwrite           = true;
figure_dir          = {'lick_events', 'ambient_light'};

% Lick detection
lick_threshold      = 1.0;   % V
min_lick_interval   = 0.05;  % s

% Long/Short classification threshold (cm)
position_threshold_cm = 90;

% Raw traces selection/export
n_raw_per_group     = 20;    % total per group
raw_per_page        = 10;    % exported per page

% Raw trace axes
raw_ylim_v = [0 3];
raw_yticks = 0:1:3;

% Raw trace layout (10 per page -> 5x2)
raw_tile_rows = 5;
raw_tile_cols = 2;

% Requested fixed axis limits for rasters
time_xlim_s     = [0 100];
pos_long_xlim   = [0 120];
pos_short_xlim  = [60 120];
ttg_xlim_pct    = [0 100];

% Tick/grid requirements
time_xticks_s       = 0:10:100;
pos_long_xticks_cm  = 0:10:120;
pos_short_xticks_cm = 60:10:120;

%%
ctl = RC2Analysis();
probe_ids = ctl.get_probe_ids(experiment_groups{:});
ctl.setup_figures(figure_dir, save_figs);

fprintf('Found %d probe(s) for experiment group(s): %s\n', length(probe_ids), strjoin(experiment_groups, ', '));
if isempty(probe_ids)
    warning('No probes found for these experiment_groups. Nothing to plot.');
    return;
end

% ---------------------------
% Cache trials per probe
% ---------------------------
probe_cache = struct('probe_id', {}, 'trials', {});

for pid = 1:length(probe_ids)
    fprintf('Caching probe %d/%d: %s\n', pid, length(probe_ids), probe_ids{pid});
    data = ctl.load_formatted_data(probe_ids{pid});

    all_trials = {};
    for kk = 1:length(trial_group_labels)
        trials = data.get_trials_with_trial_group_label(trial_group_labels{kk});
        all_trials = [all_trials; trials]; %#ok<AGROW>
    end

    probe_cache(pid).probe_id = probe_ids{pid};
    probe_cache(pid).trials   = all_trials;
end

% ---------------------------
% (A) Raw lick traces per probe: up to 20 Long + up to 20 Short
% Exported as pages (10 per page) with correct x-axis numbers
% ---------------------------
fprintf('\nCreating raw lick trace figures (export layout fixed)...\n');

for pid = 1:numel(probe_cache)
    probe_id   = probe_cache(pid).probe_id;
    all_trials = probe_cache(pid).trials;

    [idx_long, idx_short] = classify_long_short_trials(all_trials, position_threshold_cm);

    sel_long  = pick_n_indices(idx_long,  n_raw_per_group);
    sel_short = pick_n_indices(idx_short, n_raw_per_group);

    if ~isempty(sel_long)
        export_raw_trace_pages(ctl, probe_id, all_trials, sel_long, ...
            'Long', [0 0 1], lick_threshold, min_lick_interval, ...
            raw_per_page, raw_tile_rows, raw_tile_cols, raw_ylim_v, raw_yticks, ...
            show_figs, save_figs);
    end

    if ~isempty(sel_short)
        export_raw_trace_pages(ctl, probe_id, all_trials, sel_short, ...
            'Short', [1 0 0], lick_threshold, min_lick_interval, ...
            raw_per_page, raw_tile_rows, raw_tile_cols, raw_ylim_v, raw_yticks, ...
            show_figs, save_figs);
    end

    fprintf('  %s: Long=%d (selected %d), Short=%d (selected %d)\n', ...
        probe_id, numel(idx_long), numel(sel_long), numel(idx_short), numel(sel_short));
end

% ---------------------------
% (B) Combined lick-event raster figure (all probes)
% ---------------------------
fprintf('\nCreating combined lick event raster figure...\n');

n_probes = numel(probe_ids);
fig_combined = ctl.figs.a4figure('landscape');

if show_figs
    set(fig_combined, 'Visible', 'on');
    figure(fig_combined);
end

for pid = 1:n_probes
    probe_id   = probe_cache(pid).probe_id;
    all_trials = probe_cache(pid).trials;
    n_trials   = numel(all_trials);

    fprintf('  Raster probe %d/%d: %s (%d trials)\n', pid, n_probes, probe_id, n_trials);

    [idx_long, idx_short] = classify_long_short_trials(all_trials, position_threshold_cm);

    lick_times_per_trial     = cell(n_trials, 1);
    lick_positions_per_trial = cell(n_trials, 1);
    lick_ttg_per_trial       = cell(n_trials, 1);

    for trial_idx = 1:n_trials
        trial = all_trials{trial_idx};

        try
            lick_signal = trial.lick(:);
            trial_time  = trial.probe_t(:);
        catch
            continue;
        end
        if numel(trial_time) < 2 || numel(lick_signal) ~= numel(trial_time)
            continue;
        end
        trial_time_rel = trial_time - trial_time(1);

        lick_events_abs = detect_lick_events(lick_signal, lick_threshold, trial_time, min_lick_interval);
        if isempty(lick_events_abs)
            continue;
        end

        lick_times = lick_events_abs - trial_time(1);
        lick_times_per_trial{trial_idx} = lick_times;

        try
            pos_full = safe_trial_position(trial, []);
            if isempty(pos_full) || numel(pos_full) ~= numel(trial_time_rel)
                lick_positions_per_trial{trial_idx} = nan(size(lick_times));
            else
                lick_positions = interp1(trial_time_rel, pos_full(:), lick_times, 'linear', nan);
                lick_positions_per_trial{trial_idx} = lick_positions;
            end
        catch
            lick_positions_per_trial{trial_idx} = nan(size(lick_times));
        end

        lick_ttg_per_trial{trial_idx} = compute_ttg_at_licks(trial, lick_times);
    end

    % Column 1: Time
    subplot(n_probes, 4, (pid-1)*4 + 1);
    hold on;
    plot_raster_onegroup(lick_times_per_trial, idx_long,  [0 0 1]);
    plot_raster_onegroup(lick_times_per_trial, idx_short, [1 0 0]);
    if pid == 1; title('Lick Events vs. Time'); end
    ylabel(sprintf('%s\nTrial #', probe_id), 'Interpreter','none');
    if pid == n_probes; xlabel('Time from trial start (s)'); end
    ylim([0.5, n_trials + 0.5]);
    xlim(time_xlim_s);
    xticks(time_xticks_s);
    grid on; box on;
    hold off;

    % Column 2: Position (Long)
    subplot(n_probes, 4, (pid-1)*4 + 2);
    hold on;
    plot_raster_onegroup(lick_positions_per_trial, idx_long, [0 0 1]);
    if pid == 1; title('Lick Events vs. Position (Long)'); end
    ylabel('Trial #');
    if pid == n_probes; xlabel('Position (cm)'); end
    ylim([0.5, n_trials + 0.5]);
    xlim(pos_long_xlim);
    xticks(pos_long_xticks_cm);
    grid on; box on;
    hold off;

    % Column 3: Position (Short)
    subplot(n_probes, 4, (pid-1)*4 + 3);
    hold on;
    plot_raster_onegroup(lick_positions_per_trial, idx_short, [1 0 0]);
    if pid == 1; title('Lick Events vs. Position (Short)'); end
    ylabel('Trial #');
    if pid == n_probes; xlabel('Position (cm)'); end
    ylim([0.5, n_trials + 0.5]);
    xlim(pos_short_xlim);
    xticks(pos_short_xticks_cm);
    grid on; box on;
    hold off;

    % Column 4: TTG
    subplot(n_probes, 4, (pid-1)*4 + 4);
    hold on;
    plot_raster_onegroup(lick_ttg_per_trial, idx_long,  [0 0 1]);
    plot_raster_onegroup(lick_ttg_per_trial, idx_short, [1 0 0]);
    if pid == 1; title('Lick Events vs. Time-to-Goal'); end
    ylabel('Trial #');
    if pid == n_probes; xlabel('Time-to-goal (%)'); end
    ylim([0.5, n_trials + 0.5]);
    xlim(ttg_xlim_pct);
    xticks(0:10:100);
    grid on; box on;
    hold off;
end

FigureTitle(fig_combined, sprintf('Lick Events (Long=Blue, Short=Red) - %s', strjoin(experiment_groups, ', ')));

if show_figs, drawnow; end
if save_figs, ctl.figs.save_fig_to_join(); end

% ---------------------------
% Join figures
% ---------------------------
if save_figs
    fprintf('\nJoining figures into single PDF...\n');
    fname = sprintf('lick_events_%s.pdf', strjoin(experiment_groups, '_'));
    ctl.figs.join_figs(fname, overwrite);

    if ~show_figs
        ctl.figs.clear_figs();
    end
end

fprintf('\n=== All probes processed ===\n');

%% ---------------------------
% Helper functions
% ---------------------------

function export_raw_trace_pages(ctl, probe_id, all_trials, sel_idx, group_name, group_color, ...
                                lick_threshold, min_lick_interval, raw_per_page, nrows, ncols, ...
                                ylim_v, yticks_v, show_figs, save_figs)

    n_total = numel(sel_idx);
    n_pages = ceil(n_total / raw_per_page);

    for p = 1:n_pages
        i1 = (p-1)*raw_per_page + 1;
        i2 = min(p*raw_per_page, n_total);
        page_idx = sel_idx(i1:i2);

        fig = ctl.figs.a4figure('portrait');

        % Export/print consistency
        set(fig, 'Color','w', 'Renderer','painters');
        set(fig, 'PaperPositionMode','auto');

        if show_figs
            set(fig, 'Visible', 'on');
            figure(fig);
        end

        tl = tiledlayout(fig, nrows, ncols, 'TileSpacing','compact', 'Padding','none');

        % Make the tile area use the page
        tl.OuterPosition = [0.03 0.06 0.94 0.88];

        sgtitle(tl, sprintf('Raw Lick Signals (%s, page %d/%d, n=%d) - %s', ...
            group_name, p, n_pages, n_total, probe_id), 'Interpreter','none', 'FontWeight','bold');

        n_tiles_used  = numel(page_idx);
        last_row_used = ceil(n_tiles_used / ncols);

        for r = 1:n_tiles_used
            trial_idx = page_idx(r);
            trial = all_trials{trial_idx};

            [t_rel, lick_signal, tt_abs] = get_trial_lick_and_time(trial);

            ax = nexttile(tl, r);
            plot(ax, t_rel, lick_signal, 'k-', 'LineWidth', 0.8); hold(ax,'on');

            % Threshold line (no giant yline label)
            yline(ax, lick_threshold, 'r--', 'LineWidth', 1.0);

            % Small threshold text (export-stable)
            text(ax, 0.99, 0.92, sprintf('Threshold = %.1fV', lick_threshold), ...
                'Units','normalized', 'HorizontalAlignment','right', ...
                'Color','r', 'FontWeight','bold', 'FontSize', 8);

            % Lick events as vertical lines
            lick_events_abs = detect_lick_events(lick_signal, lick_threshold, tt_abs, min_lick_interval);
            if ~isempty(lick_events_abs)
                lick_events_rel = lick_events_abs - tt_abs(1);
                for lt = lick_events_rel(:)'
                    xline(ax, lt, '-', 'Color', group_color, 'LineWidth', 0.6, 'Alpha', 0.5);
                end
            end

            ylim(ax, ylim_v);
            yticks(ax, yticks_v);

            % Correct row/col computation (ROW-MAJOR)
            row_i = ceil(r / ncols);
            col_i = mod(r-1, ncols) + 1;

            % Y labels only on left column (hide y-axis on right column)
            if col_i == 1
                ylabel(ax, 'Lick (V)');
                ax.YAxis.Visible = 'on';
            else
                ax.YAxis.Visible = 'off';
            end

            % ---- FIX: Do NOT touch XTickLabel. Hide/show the entire XAxis. ----
            if row_i == last_row_used
                xlabel(ax, 'Time (s)');
                ax.XAxis.Visible = 'on';     % ensures tick numbers appear in PDF
            else
                ax.XAxis.Visible = 'off';    % cleaner; preserves auto ticks internally
            end

            title(ax, sprintf('%s - Trial %d (ID: %s)', group_name, trial_idx, safe_trial_id(trial, trial_idx)), ...
                'Interpreter','none');

            set(ax, 'FontSize', 8);
            grid(ax, 'on'); box(ax, 'on');
            hold(ax,'off');
        end

        if show_figs, drawnow; end
        if save_figs, ctl.figs.save_fig_to_join(); end
    end
end

function [idx_long, idx_short] = classify_long_short_trials(all_trials, position_threshold_cm)
    n_trials = numel(all_trials);
    is_long = false(n_trials,1);

    for t = 1:n_trials
        trial = all_trials{t};
        try
            mm = trial.motion_mask();
            pos = safe_trial_position(trial, mm);
            if isempty(pos)
                is_long(t) = false;
            else
                is_long(t) = max(pos) > position_threshold_cm;
            end
        catch
            is_long(t) = false;
        end
    end

    idx_long  = find(is_long);
    idx_short = find(~is_long);
end

function pos = safe_trial_position(trial, motion_mask)
    pos = [];
    try
        if nargin >= 2 && ~isempty(motion_mask)
            try
                pos = trial.position(motion_mask);
                return;
            catch
                pfull = trial.position();
                if numel(pfull) == numel(motion_mask)
                    pos = pfull(motion_mask);
                else
                    pos = pfull(:);
                end
                return;
            end
        else
            pos = trial.position();
        end
    catch
        pos = [];
    end
end

function sel = pick_n_indices(idx, n)
    idx = idx(:);
    if isempty(idx)
        sel = [];
        return;
    end
    if numel(idx) <= n
        sel = idx;
        return;
    end
    ii = round(linspace(1, numel(idx), n));
    ii(ii < 1) = 1;
    ii(ii > numel(idx)) = numel(idx);
    sel = idx(ii);
end

function [t_rel, lick_signal, tt_abs] = get_trial_lick_and_time(trial)
    lick_signal = trial.lick(:);
    tt_abs = trial.probe_t(:);
    t_rel = tt_abs - tt_abs(1);
end

function tid = safe_trial_id(trial, fallback_idx)
    tid = sprintf('%d', fallback_idx);
    try
        if isprop(trial,'trial_id') && ~isempty(trial.trial_id)
            tid = sprintf('%d', trial.trial_id);
        elseif isfield(trial,'trial_id') && ~isempty(trial.trial_id)
            tid = sprintf('%d', trial.trial_id);
        end
    catch
    end
end

function plot_raster_onegroup(x_per_trial, idx_trials, rgb)
    for k = 1:numel(idx_trials)
        t = idx_trials(k);
        x = x_per_trial{t};
        if ~isempty(x)
            x = x(:);
            x = x(isfinite(x));
            if ~isempty(x)
                plot(x, t*ones(size(x)), '.', 'Color', rgb, 'MarkerSize', 4);
            end
        end
    end
end

function lick_ttg = compute_ttg_at_licks(trial, lick_times_rel)
    lick_ttg = nan(numel(lick_times_rel),1);
    try
        mm = trial.motion_mask();
        tt = trial.probe_t(:);

        if isempty(mm) || numel(mm) ~= numel(tt)
            return;
        end

        motion_times = tt(mm);
        if numel(motion_times) < 2
            return;
        end

        motion_rel = motion_times - motion_times(1);
        T = motion_rel(end);
        if T <= 0
            return;
        end

        ttg_abs  = T - motion_rel;
        ttg_norm = 100 * (ttg_abs / T);

        for i = 1:numel(lick_times_rel)
            lt = lick_times_rel(i);
            [~, nearest_motion_idx] = min(abs(motion_rel - lt));
            if isempty(nearest_motion_idx) || ~isfinite(nearest_motion_idx)
                lick_ttg(i) = nan;
            else
                lick_ttg(i) = ttg_norm(nearest_motion_idx);
            end
        end
    catch
    end
end

function lick_events = detect_lick_events(lick_signal, threshold, time_vector, min_interval)
    above_threshold = lick_signal > threshold;
    rising_edges = find(diff([0; above_threshold]) == 1);

    if isempty(rising_edges)
        lick_events = [];
        return;
    end

    lick_events = time_vector(rising_edges);

    if numel(lick_events) > 1
        intervals = diff(lick_events);
        keep_idx = [true; intervals >= min_interval];
        lick_events = lick_events(keep_idx);
    end
end
