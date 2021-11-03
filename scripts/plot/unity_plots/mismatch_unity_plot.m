%%script for producing a mismatch unity plot

%%
% experiment_groups       = {'mismatch_nov20'};
% trial_group_labels      = {'RVT_gain_up', 'RV_gain_up'};
% marker_style            = {'o', 'o'};



%%
experiment_groups       = {'mismatch_darkness_oct21'};
trial_group_labels      = {'RT_gain_up'};
marker_style            = {'o'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'mismatch_unity_plot', 'mismatch_darkness_oct21'};



%%
ctl                     = RC2Analysis();
ctl.setup_figures(figure_dir, save_figs);

probe_ids               = ctl.get_probe_ids(experiment_groups{:});
mm                      = MismatchAnalysis();

x_all                   = cell(1, length(trial_group_labels));
y_all                   = cell(1, length(trial_group_labels));
p_val                   = cell(1, length(trial_group_labels));
direction               = cell(1, length(trial_group_labels));

for ii = 1 : length(probe_ids)
    
    data        = ctl.load_formatted_data(probe_ids{ii});
    clusters    = data.selected_clusters;
    
    for jj = 1 : length(trial_group_labels)
        
        trials      = data.get_trials_with_trial_group_label(trial_group_labels{jj});
        
        for kk = 1 : length(clusters)
            
            x_all{jj}(end+1) = mm.get_avg_baseline_fr(clusters(kk), trials);
            y_all{jj}(end+1) = mm.get_avg_response_fr(clusters(kk), trials);
            [~, p_val{jj}(end+1), direction{jj}(end+1)] = mm.is_response_significant(clusters(kk), trials);
        end
    end
end



%%
h_fig                   = a4figure();
plot_array              = PlotArray(3, 2);
u                       = UnityPlotPopulation.empty();


for ii = 1 : length(trial_group_labels)
    
    pos         = plot_array.get_position(ii);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    u(end+1)   = UnityPlotPopulation(x_all{ii}, ...
        y_all{ii}, ...
        p_val{ii}, ...
        direction{ii}, ...
        h_ax);
    
    u(end).marker_style = marker_style{ii};
    
    u(end).plot();
    
    u(end).xlabel('Baseline FR (Hz)');
    u(end).ylabel('Response FR (Hz)');
    
    u(end).title(trial_group_labels{ii});
    u(end).add_histogram(1);
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for ii = 1 : length(u)
    u(ii).xlim([m, M]);
end


ctl.figs.save_fig('unity_plot', overwrite);
