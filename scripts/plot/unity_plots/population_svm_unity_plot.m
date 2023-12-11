%%script for producing a unity plot for stationary vs. motion for each probe recording and trial group


% experiment_groups       = {'visual_flow'};
% trial_group_labels      = {'RVT', 'RV', {'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};
% marker_style            = {'o', 'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'visual_flow', 'stationary_vs_motion', 'population'};


% experiment_groups       = {'visual_flow','mismatch_nov20','mismatch_jul21'};
% trial_group_labels      = {{'RVT','RVT_gain_up'}, {'RV','RV_gain_up'}};
% marker_style            = {'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'visual_flow', 'stationary_vs_motion', 'population'};

%%
experiment_groups       = {'darkness','mismatch_darkness_oct21'};
trial_group_labels      = {{'RT', 'RT_gain_up'}, 'R', {'T_bank', 'T_RT', 'T_R', 'T'}};
marker_style            = {'o', 'o', 'o'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'unity_plots', 'darkness', 'stationary_vs_motion', 'population'};



%%
% experiment_groups       = {'mismatch_nov20','mismatch_jul21'};
% trial_group_labels      = {'RVT_gain_up', 'RV_gain_up', 'R', 'T'};
% marker_style            = {'o', 'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'mismatch_jul21', 'stationary_vs_motion', 'population'};



%%
% experiment_groups       = {'mismatch_darkness_oct21'};
% trial_group_labels      = {'R', 'T', 'RT_gain_up'};
% marker_style            = {'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'mismatch_darkness_oct21', 'stationary_vs_motion', 'population'};



%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

% preallocate cell arrays
x_all                   = cell(1, length(probe_ids));
y_all                   = cell(1, length(probe_ids));
p_val                   = cell(1, length(probe_ids));
direction               = cell(1, length(probe_ids));

for ii = 1 : length(probe_ids)
    
    data                = ctl.load_formatted_data(probe_ids{ii});
    clusters            = data.VISp_clusters;
    cluster_ids         = data.VISp_cluster_ids;
    
    for jj = 1 : length(trial_group_labels)
        
        % skip if the trial group label is not in the experiment
        if ~data.check_trial_group(trial_group_labels{jj})
            continue
        end
        
        % loop over clusters and get the motion and stationary firing rates
        % for each trial
        for kk = 1 : length(cluster_ids)
            
            x = data.stationary_fr_for_trial_group(cluster_ids(kk), trial_group_labels{jj});
            y = data.motion_fr_for_trial_group(cluster_ids(kk), trial_group_labels{jj});
            
            x_all{ii}{jj}(kk) = median(x);
            y_all{ii}{jj}(kk) = median(y);
            
            [~, p_val{ii}{jj}(kk), direction{ii}{jj}(kk)] = ...
                data.is_stationary_vs_motion_significant(cluster_ids(kk), trial_group_labels{jj});
        end
    end
end



%% Unity plot
ctl.setup_figures(figure_dir, save_figs);

for ii = 1 : length(probe_ids)
        
    h_fig                   = ctl.figs.a4figure();
    plot_array              = PlotArray(3, 2);
    u                       = UnityPlotPopulation.empty();

    for jj = 1 : length(trial_group_labels)

        pos         = plot_array.get_position(jj);
        h_ax        = axes('units', 'centimeters', 'position', pos);

        u(end+1)   = UnityPlotPopulation(x_all{ii}{jj}, ...
                                        y_all{ii}{jj}, ...
                                        p_val{ii}{jj}, ...
                                        direction{ii}{jj}, ...
                                        h_ax);

        u(end).marker_style = marker_style{jj};

        u(end).plot();

        u(end).xlabel('Baseline FR (Hz)');
        u(end).ylabel('Response FR (Hz)');

        u(end).title(trial_group_labels{jj});
%             u(end).add_histogram(1);
    end

    % sync across all axes
    m               = min([u(:).min]);
    M               = max([u(:).max]);

    for kk = 1 : length(u)
        u(kk).xlim([m, M]);
    end

     % give the page a title
    FigureTitle(h_fig, probe_ids{ii});
        
    ctl.figs.save_fig(probe_ids{ii});
end




%% pooled
h_fig                   = ctl.figs.a4figure();
plot_array              = PlotArray(3, 2);
u                       = UnityPlotPopulation.empty();


for jj = 1 : length(trial_group_labels)
    
    x_pooled    = cellfun(@(x)(x{jj}), x_all, 'UniformOutput', false);
    y_pooled    = cellfun(@(x)(x{jj}), y_all, 'UniformOutput', false);
    p_pooled    = cellfun(@(x)(x{jj}), p_val, 'UniformOutput', false);
    direction_pooled = cellfun(@(x)(x{jj}), direction, 'UniformOutput', false);
    
    pos         = plot_array.get_position(jj);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    u(end+1)   = UnityPlotPopulation([x_pooled{:}], ...
        [y_pooled{:}], ...
        [p_pooled{:}], ...
        [direction_pooled{:}], ...
        h_ax);
    
    u(end).marker_style = marker_style{jj};
   
    u(end).plot();
    
    u(end).xlabel('Baseline FR (Hz)');
    u(end).ylabel('Response FR (Hz)');
    
    u(end).title(trial_group_labels{jj});
end

% sync across all axes
m               = min([u(:).min]);
M               = max([u(:).max]);

for kk = 1 : length(u)
    u(kk).xlim([m, M]);
end

% give the page a title
FigureTitle(h_fig, 'pooled');

ctl.figs.save_fig('pooled');
