% Plot stationary vs. motion unity plot for an experiment group and trial
% group as well as a modulation index vs. depth plot
%
%   Specify options:
%
%       experiment_groups:      Will generate unity plots/MI vs. depth plots including all
%                               clusters in the probe recordings 
%                               in the specified experiment group. e.g. one of:
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21'
%                               Should be a cell array of strings with each
%                               entry an experiment group.
%
%       trial_group_labels:     Will generate a unity plot/MI vs. depth for each of the
%                               trial groups specified. 
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%                               e.g. {'R', {'T_bank', 'T_RT', 'T_R'}, 'RT'}
%                               will create three unity/MIvdepth plots, the first
%                               combining across 'R' (running) trials, the
%                               second combining across any trial of
%                               'T_bank', 'T_RT' or 'T_R' type, and the
%                               third combining across 'RT'
%                               (running+translation) trials.
%
%       marker_type:            the type of marker to use on the plot (any
%                               of those accepted by matlab) e.g. 'o' or
%                               'v'. Should be a cell array of the same
%                               length as `trial_group_labels`.
%
%       save_figs:              true or false, whether to save the figures to pdf
%
%       overwrite:              true or false. If figure pdf's already exist,
%                               whether to overwrite 
%       
%       figure_dir:             cell array of strings specifying which
%                               directory to save pdf's. The directory will
%                               be relative to the directory specified by
%                               path_config.figure_dir (in
%                               `path_config.m`), so that {'one', 'two',
%                               'three'} will save .pdfs to:
%                               <path_config.figure_dir>\one\two\three\      
%
%       mi_figure_dir:          same as figure_dir but where to store the
%                               MI vs. depth plots
%
% If `save_figs` is true, one pdf is saved for each probe recording, and a
% separate one for clusters from all probe recordings pooled.


%% experiments for visual flow
% experiment_groups       = {'visual_flow'};
% trial_group_labels      = {'RVT', 'RV', {'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};
% marker_style            = {'o', 'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'visual_flow', 'stationary_vs_motion', 'population'};
% mi_figure_dir           = {'mi_vs_depth', 'visual_flow', 'stationary_vs_motion', 'population'};

% experiment_groups       = {'visual_flow','mismatch_nov20','mismatch_jul21'};
% trial_group_labels      = {{'RVT','RVT_gain_up'}, {'RV','RV_gain_up'}};
% marker_style            = {'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'visual_flow', 'stationary_vs_motion', 'population'};
% mi_figure_dir           = {'mi_vs_depth', 'visual_flow', 'stationary_vs_motion', 'population'};


%% experiments performed in darkness
% experiment_groups       = {'darkness','mismatch_darkness_oct21'};
% trial_group_labels      = {{'T_bank', 'T_RT', 'T_R', 'T'}};
% marker_style            = {'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'darkness', 'stationary_vs_motion', 'population'};
% mi_figure_dir           = {'mi_vs_depth', 'darkness', 'stationary_vs_motion', 'population'};


%% mismatch experiments
experiment_groups       = {'mismatch_jul21'};
trial_group_labels      = {'RVT_gain_up', 'RV_gain_up', 'R', 'T'};
marker_style            = {'o', 'o', 'o', 'o'};
save_figs               = true;
overwrite               = true;
figure_dir              = {'unity_plots', 'mismatch_jul21', 'stationary_vs_motion', 'population'};
mi_figure_dir           = {'mi_vs_depth', 'mismatch_jul21', 'stationary_vs_motion', 'population'};



%%
% experiment_groups       = {'mismatch_darkness_oct21'};
% trial_group_labels      = {'R', 'T', 'RT_gain_up'};
% marker_style            = {'o', 'o', 'o'};
% save_figs               = true;
% overwrite               = true;
% figure_dir              = {'unity_plots', 'mismatch_darkness_oct21', 'stationary_vs_motion', 'population'};
% mi_figure_dir           = {'mi_vs_depth', 'mismatch_darkness_oct21', 'stationary_vs_motion', 'population'};



%%
ctl                     = RC2Analysis();
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

% preallocate cell arrays
x_all                   = cell(1, length(probe_ids));
y_all                   = cell(1, length(probe_ids));
p_val                   = cell(1, length(probe_ids));
direction               = cell(1, length(probe_ids));
modulation_index        = cell(1, length(probe_ids));
relative_depth          = cell(1, length(probe_ids));
layer                   = cell(1, length(probe_ids));
anatomies               = cell(1, length(probe_ids));
anatomy_idx             = cell(1, length(probe_ids));
cortical_position       = cell(1, length(probe_ids));

for ii = 1 : length(probe_ids)
    
    data                = ctl.load_formatted_data(probe_ids{ii});
    clusters            = data.VISp_clusters;
    cluster_ids         = data.VISp_cluster_ids;
    anatomies{ii}       = data.anatomy;
    
    shank_ids           = cellfun(@(x)(x.shank_id), anatomies{ii});
    
    % get layer and relative position in layer for each cluster
    for jj = 1 : length(cluster_ids)
        
        % which of the anatomies is this cluster part of
        idx = find(shank_ids == clusters(jj).shank_id);
        
        % get the cortical position in this anatomy
        cortical_position{ii}(jj) = anatomies{ii}{idx}.from_tip_to_from_pia(clusters(jj).distance_from_probe_tip);
        
        % store the anatomy index
        anatomy_idx{ii}(jj) = idx;
        
        % layer and relative depth of cluster in layer
        [relative_depth{ii}(jj), layer{ii}{jj}] = data.get_relative_layer_depth_of_cluster(cluster_ids(jj));
    end
    
    
    for jj = 1 : length(trial_group_labels)
        
        % skip if the trial group label is not in the experiment
        if ~data.check_trial_group(trial_group_labels{jj})
            continue
        end
        
        % loop over clusters and get the motion and stationary firing rates
        % for each trial
        for kk = 1 : length(cluster_ids)
            
            [~, p, d, x, y] = data.is_stationary_vs_motion_significant(cluster_ids(kk), trial_group_labels{jj});
            
            modulation_index{ii}{jj}(kk) = (y - x)/(y + x);
            
            x_all{ii}{jj}(kk) = x;
            y_all{ii}{jj}(kk) = y;
            p_val{ii}{jj}(kk) = p;
            direction{ii}{jj}(kk) = d;
        end
    end
end

avg_anatomy = AverageAnatomy(anatomies);
boundary_position = avg_anatomy.average_VISp_boundaries_from_pia();
for ii = length(probe_ids) : -1 : 1
    averaged_cortical_position{ii} = avg_anatomy.from_pia_using_relative_position(relative_depth{ii}, layer{ii});
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
        
    ctl.figs.save_fig(probe_ids{ii}, overwrite);
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

ctl.figs.save_fig('pooled', overwrite);





%% MI vs. depth plot
ctl.setup_figures(mi_figure_dir, save_figs);

for ii = 1 : length(probe_ids)
    
    % get the number of shanks in this recording
    shank_ids               = cellfun(@(x)(x.shank_id), anatomies{ii});
    
    for shank_ii = 1 : length(shank_ids)
        
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(3, 2);
        mi                      = MIDepthPlot.empty();
        
        % which of the clusters are on this shank
        clust_idx = anatomy_idx{ii} == shank_ii;
        
        for jj = 1 : length(trial_group_labels)
            
            pos         = plot_array.get_position(jj);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            mi(end+1)   = MIDepthPlot(modulation_index{ii}{jj}(clust_idx), ...
                                      p_val{ii}{jj}(clust_idx), ...
                                      direction{ii}{jj}(clust_idx), ...
                                      cortical_position{ii}(clust_idx), ...
                                      anatomies{ii}{shank_ii}.VISp_boundaries_from_pia_flat, ...
                                      anatomies{ii}{shank_ii}.VISp_layers, ...
                                      h_ax);
            
            mi(end).plot();
            
            if jj == length(trial_group_labels)
                mi(end).print_layers();
                mi(end).xlabel('MI');
            end
            
            mi(end).title(trial_group_labels{jj});
            
            mi(end).title(trial_group_labels{jj});
        end
        
        % give the page a title
        title_str = sprintf('%s, shank %i', probe_ids{ii}, shank_ids(shank_ii));
        FigureTitle(h_fig, title_str);
        
        ctl.figs.save_fig_to_join()
    end
    
    ctl.figs.join_figs(sprintf('%s.pdf', probe_ids{ii}), overwrite);
    ctl.figs.clear_figs();
end




%% pooled
h_fig                   = ctl.figs.a4figure();
plot_array              = PlotArray(3, 2);
u                       = MIDepthPlot.empty();

% pool all cortical positions
cortical_pos_pooled    = cat(1, averaged_cortical_position{:});
 

for jj = 1 : length(trial_group_labels)
    
    mi_pooled    = cellfun(@(x)(x{jj}), modulation_index, 'UniformOutput', false);
    y_pooled    = cellfun(@(x)(x{jj}), y_all, 'UniformOutput', false);
    p_pooled    = cellfun(@(x)(x{jj}), p_val, 'UniformOutput', false);
    direction_pooled = cellfun(@(x)(x{jj}), direction, 'UniformOutput', false);
    
    pos         = plot_array.get_position(jj);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    u(end+1)   = MIDepthPlot([mi_pooled{:}], ...
                             [p_pooled{:}], ...
                             [direction_pooled{:}], ...
                             cortical_pos_pooled, ...
                             avg_anatomy.average_VISp_boundaries_from_pia, ...
                             avg_anatomy.VISp_layers, ...
                             h_ax);
    
    u(end).marker_style = marker_style{jj};
   
    u(end).plot();
    
    u(end).xlabel('Baseline FR (Hz)');
    u(end).ylabel('Response FR (Hz)');
    
    u(end).title(trial_group_labels{jj});
end

% give the page a title
FigureTitle(h_fig, 'pooled');

ctl.figs.save_fig('pooled', overwrite);
