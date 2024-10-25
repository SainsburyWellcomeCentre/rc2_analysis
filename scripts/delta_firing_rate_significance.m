ctl                 = RC2Analysis();

% trial_group_labels = {{'V_RVT', 'V_RV', 'V'}, {'VT_RVT', 'VT_RV', 'VT'}};
trial_group_labels = {{'RVT'}, {'RV'}}; 

probe_ids   = ctl.get_probe_ids('visual_flow');%, 'passive_same_luminance');

pop_fr_delta = {};
for ll = 1 : length(trial_group_labels)
    c           = 0;
    x_med       = [];
    y_med       = [];
    direction   = [];

    for ii = 1 : length(probe_ids)  
        data   = ctl.load_formatted_data(probe_ids{ii});
        clusters  = data.VISp_clusters();

        for jj = 1 : length(clusters)
            c = c + 1;
            [~, ~, direction(c), x_med(c), y_med(c)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{ll});
        end
    end

    pop_fr_delta{ll}  = y_med - x_med;
end

[p, h, stats] = signrank(pop_fr_delta{1}, pop_fr_delta{2})
prctile(pop_fr_delta{1}, [25 50 75])
prctile(pop_fr_delta{2}, [25 50 75])
