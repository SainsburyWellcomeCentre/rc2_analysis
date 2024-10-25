% Print a variety of statistics related to mismatch trials.

ctl                 = RC2Analysis();
probe_ids           = ctl.get_probe_ids('mismatch_nov20');
trial_group_label   = 'RVT_gain_down';


% Data

mm                  = MismatchAnalysis();
mm.method           = 'anova';

c                   = 0;
p_val               = [];
direction           = [];
avg_baseline        = [];
avg_response        = [];
baseline_normal     = [];
response_normal     = [];
spike_class         = {};

for ii = 1 : length(probe_ids)
    this_data       = ctl.load_formatted_data(probe_ids{ii});
    clusters        = this_data.VISp_clusters();
    
    trials          = this_data.get_trials_with_trial_group_label(trial_group_label);
    
    % remove trials with small mismatch period
    trials          = remove_invalid_mm_trials(trials);
    
    for jj = 1 : length(clusters)
        
        c = c + 1;
        avg_baseline(c) = mm.get_avg_baseline_fr(clusters(jj), trials);
        avg_response(c) = mm.get_avg_response_fr(clusters(jj), trials);
        [~, p_val(c), direction(c)] = mm.is_response_significant(clusters(jj), trials);
        baseline_normal(c) = mm.is_baseline_normal(clusters(jj), trials);
        response_normal(c) = mm.is_response_normal(clusters(jj), trials);
        spike_class{c} = clusters(jj).spiking_class;
    end
end


% Print
fprintf('\n\nFigure 2G, mismatch response R:VF\n');
print_unity_plot_stats(avg_baseline, avg_response, direction, spike_class)
fprintf(' Fraction in which we reject normality in baseline: %.2f%% (%i/%i)\n', 100*sum(~baseline_normal)/length(baseline_normal), sum(~baseline_normal), length(baseline_normal));
fprintf(' Fraction in which we reject normality in response: %.2f%% (%i/%i)\n', 100*sum(~response_normal)/length(response_normal), sum(~response_normal), length(response_normal));




%% No change trials / matched
bsl_limits = [-0.4, 0]; % pre-slip window
rsp_limits = [0, 0.4]; % slip window

load('C:\Users\lee\Documents\mvelez\mvelez_ms_figures\each_figure\figure_s4\running_around_mismatch_matched_trials', ...
                                               'lib_git', ...                                               
                                               'common_t', ...
                                               'decrease_traces', ...
                                               'no_change_traces', ...
                                               'matched_decrease_traces', ...
                                               'decrease_traces_fr', ...
                                               'no_change_traces_fr', ...
                                               'matched_decrease_traces_fr');
                                           
bsl_idx = common_t >= bsl_limits(1) & common_t < bsl_limits(2);
rsp_idx = common_t >= rsp_limits(1) & common_t < rsp_limits(2);

traces = [matched_decrease_traces_fr{:}];

bsl_fr = mean(traces(bsl_idx, :), 1);
rsp_fr = mean(traces(rsp_idx, :), 1);

perc_pre_slip = prctile(bsl_fr, [50 25 75])    
perc_slip = prctile(rsp_fr, [50 25 75])  

delta_fr = rsp_fr - bsl_fr;
pop_delta_fr_avg = median(delta_fr);
pop_delta_fr_q = prctile(delta_fr, [50 25 75])
delta_fr_p_val = signrank(rsp_fr, bsl_fr)


