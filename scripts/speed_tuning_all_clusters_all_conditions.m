% Plot speed tuning curves for all cortical (VISp) clusters in all
% experiments for the trial types in which the activity of that cluster was
% recorded. 
%
% If specified, a .pdf is created for each experiment containing the
% tuning curve information for each VISp cluster in that experiment.
%
%   Specify options:
%
%       experiment_groups:      Will generate plots for all (VISp) clusters
%                               and all probe recordings 
%                               in the specified experiment groups.
%                                   'darkness',
%                                   'visual_flow',
%                                   'mismatch_nov20',
%                                   'mismatch_jul21',
%                                   'mismatch_darkness_oct21'
%                               Should be a cell array of strings with each
%                               entry an experiment group
%
%       trial_group_labels:     Will generate plots with speed tuning
%                               curves for the trials specified. 
%                               Should be a cell array, with each entry
%                               either a string specifying a trial group,
%                               or a cell array of strings specifying
%                               multiple trial groups.
%                               e.g. {'R', {'T_bank', 'T_RT', 'T_R'}, 'RT'}
%                               will plot three sets of tuning curves, the first
%                               for all 'R' (running) trials, the
%                               second for all trials of any of 
%                               'T_bank', 'T_RT' or 'T_R' type, and the
%                               third for all 'RT'
%                               (running+translation) trials.
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
% If `save_figs` is true, one pdf will be created for each probe recording,
% and contain a A4 page for each cluster, containing the raster data for
% several conditions.

% experiment_groups       = {'darkness', ...
%                            'visual_flow', ...
%                            'mismatch_nov20', ...
%                            'mismatch_jul21', ...
%                            'mismatch_darkness_oct21'};
% 
% trial_group_labels      = {{'RVT', 'RVT_gain_up'}...
%                            {'RV', 'RV_gain_up'}...
%                            {'VT', 'VT_RVT', 'VT_RV'}, ...
%                            {'V', 'V_RVT', 'V_RV'}, ...
%                            {'RT', 'RT_gain_up'}, ...
%                            'R', ...
%                            {'T_bank', 'T_RT', 'T_R', 'T'}};

% experiment_groups       = {'darkness','mismatch_darkness_oct21'};
% 
% trial_group_labels      = {{'T_bank', 'T_RT', 'T_R', 'T'}};

experiment_groups   = {'mismatch_nov20', 'mismatch_jul21'};
trial_group_labels   = {'RVT_gain_up'};

%V_VT_idx                = [3, 4];
                       
save_figs               = true;
overwrite               = true;
figure_dir              = {'all_tuning_curves'};



%% collect data

% main controller object
ctl                     = RC2Analysis();

% get all probe_ids occurring in the specified experiment groups
probe_ids               = ctl.get_probe_ids(experiment_groups{:});

% setup directory to save figures in if 'save_figs' is true
ctl.setup_figures(figure_dir, save_figs);

% arrays and cell arrays in which to store information for each experiment,
% VISp cluster and trial group type
tuning          = {};
p_svm           = [];
direction       = [];
probe_id        = {};
cluster_id      = [];
cluster_region  = {};
cluster_in_V_or_VT = [];

% loop over experiments
for ii = 1 : length(probe_ids)
    
    % load data for the experiment
    data        = ctl.load_formatted_data(probe_ids{ii});
    
    % return array of Cluster objects containing info about the VISp
    % clusters
    clusters    = data.VISp_clusters();
    
    % creates tuning curves for all clusters (not just VISp) for each of
    % the trials specified in trial_group_labels
    % will be cell array of length equal to 'trial_group_labels', each
    % entry is itself a structure array with each entry the tuning curve
    % info for a cluster
    tuning_curves = data.create_tuning_curves(trial_group_labels);
    
    % loop over all VISp clusters
    for jj = 1 : length(clusters)
        
        % loop over trial group types
        for kk = 1 : length(trial_group_labels)
            
            fprintf('Collecting for %i/%i, %i/%i, %i/%i\n', ...
                    kk, length(trial_group_labels), ...
                    jj, length(clusters), ...
                    ii, length(probe_ids));
            
            % if tuning_curves{kk} is empty there were no trials of the
            % specified type in the experiment
            if isempty(tuning_curves{kk})
                
                % store nan's and empty array 
                p_svm(ii, jj, kk) = nan;
                direction(ii, jj, kk) = nan;
                tuning{ii}{jj}{kk} = [];
            else
                
                % for this cluster, and this collection of trial types
                % return statistics of whether activity during stationary
                % and motion periods was significantly different
                [~, p_svm(ii, jj, kk), direction(ii, jj, kk)] = data.is_stationary_vs_motion_significant(clusters(jj).id, trial_group_labels{kk});
                
                % find the tuning curve information for the current cluster
                % in 'tuning_curves' and re-store it in 'tuning'
                cluster_idx = [tuning_curves{kk}(:).cluster_id] == clusters(jj).id;
                tuning{ii}{jj}{kk} = tuning_curves{kk}(cluster_idx);
            end
        end
        
        % extra info to plot on the figure
        probe_id{ii, jj} = probe_ids{ii};
        cluster_id(ii, jj) = clusters(jj).id;
        cluster_region{ii, jj} = clusters(jj).region_str;
        
        % whether cluster appears in V or VT trials
%         is_cluster_in_V_or_VT = cellfun(@(x)(~isempty(x)), tuning_curves);
%         is_cluster_in_V_or_VT = any(is_cluster_in_V_or_VT(V_VT_idx));
%         cluster_in_V_or_VT(ii, jj) = is_cluster_in_V_or_VT;
    end
end



%% plotting

% loop over experiments
for ii = 1 : length(probe_ids)  
    
    % loop over VISp clusters in that experiment
    for jj = 1 : length(tuning{ii})
        
        % create A4 figure and object controlling the array of axes on the
        % figure
        h_fig                   = ctl.figs.a4figure();
        plot_array              = PlotArray(4, 4);
        
        % cell arrays which will contain objects controlling the main
        % tuning curve plot and the shuffled tuning curve histogram info
        tuning_curve_plot       = {};
        shuff_histogram         = {};
        
        % loop over trial group types
        for kk = 1 : length(trial_group_labels)
            
            % if there was no information in an experiment for a particular
            % trial type then skip this trial type
            if isempty(tuning{ii}{jj}{kk})
                continue
            end
            
            % find where the axis should go and create it
            pos         = plot_array.get_position(2*kk-1);
            h_ax        = axes('units', 'centimeters', 'position', pos);
            
            % object controlling the tuning curve plot
            tuning_curve_plot{kk} = TuningCurvePlot(h_ax);
            
            % determine the colour of the plot from the stationary vs.
            % motion statistics (increase/decrease/no change)
            if p_svm(ii, jj, kk) < 0.05 && direction(ii, jj, kk) == 1
                % increase
                main_col = [0.85, 0.32, 0.1];
            elseif p_svm(ii, jj, kk) < 0.05 && direction(ii, jj, kk) == -1
                % decrease
                main_col = [0, 0.45, 0.75];
            else
                % no change
                main_col = 'k';
            end
            
            % plot the tuning curve and give axis a title
            tuning_curve_plot{kk}.main_col = main_col;
            tuning_curve_plot{kk}.plot(tuning{ii}{jj}{kk});
            title(gca, trial_group_labels{kk}, 'interpreter', 'none');
            
            % axis labels on first axis
            if kk == 1
                tuning_curve_plot{kk}.xlabel('Speed (cm/s)');
                tuning_curve_plot{kk}.ylabel('Firing rate (Hz)');
            end
             
            % find where the shuffled histogram axis should go and create
            % it 
            pos         = plot_array.get_position(2*kk);
            h_ax2        = axes('units', 'centimeters', 'position', pos);
            
            % object controlling the shuffled histogram plot
            shuff_histogram{kk} = TuningCurveHistogram(h_ax2);
            shuff_histogram{kk}.plot(tuning{ii}{jj}{kk});
        end
        
        % loop over all axes and find the min X, max X and max Y values
        mx              = inf;
        Mx              = -inf;
        My              = -inf;
        
        for kk = 1 : length(tuning_curve_plot)
            if isempty(tuning{ii}{jj}{kk})
                continue
            end
            mx = min(tuning_curve_plot{kk}.xmin, mx);
            Mx = max(tuning_curve_plot{kk}.xmax, Mx);
            My = max(tuning_curve_plot{kk}.ymax, My);
        end
        
        % again loop over axes and set the X and Y limits to the same value
        for kk = 1 : length(tuning_curve_plot)
            if isempty(tuning{ii}{jj}{kk})
                continue
            end
            tuning_curve_plot{kk}.xlim([mx, Mx]);
            tuning_curve_plot{kk}.ylim([0, My]);
        end

        % give whole figure a title (probe_id, cluster # and cluster
        % region)
        FigureTitle(h_fig, sprintf('%s, Cluster %i, %s', ...
                probe_id{ii, jj}, ...
                cluster_id(ii, jj), ...
                cluster_region{ii, jj}));
        
        % temporarily save the figure to a .pdf to join later
        ctl.figs.save_fig_to_join();
    end
    
    % join temporarily saved figures and save to a new .pdf
    fname = sprintf('%s.pdf', probe_ids{ii});
    ctl.figs.join_figs(fname, overwrite);
    % delete temporarily saved .pdfs and clear MATLAB figure windows
    ctl.figs.clear_figs();
end



%% plot average across clusters

restrict_to_V_VT_clusters = false;

% get an average tuning curve for each trial type
trial_type_tuning_store = cell(1, length(trial_group_labels));
stationary_store = cell(1, length(trial_group_labels));

% loop over experiments
for ii = 1 : length(probe_ids)  
    
    % loop over VISp clusters in that experiment
    for jj = 1 : length(tuning{ii})
        
        % loop over trial group types
        for kk = 1 : length(trial_group_labels)
            
            % if there was no information in an experiment for a particular
            % trial type then skip this trial type
            if isempty(tuning{ii}{jj}{kk})
                continue
            end
            
            % average tuning curve for this cluster
            avg_tuning = nanmean(tuning{ii}{jj}{kk}.tuning, 2);
            avg_stationary = mean(tuning{ii}{jj}{kk}.stationary_fr);
            
            if restrict_to_V_VT_clusters
                if cluster_in_V_or_VT(ii, jj)
                    trial_type_tuning_store{kk} = [trial_type_tuning_store{kk}, avg_tuning];
                    stationary_store{kk} = [stationary_store{kk}, avg_stationary];
                end
            else
                trial_type_tuning_store{kk} = [trial_type_tuning_store{kk}, avg_tuning];
                stationary_store{kk} = [stationary_store{kk}, avg_stationary];
            end
        end
    end
end

% create A4 figure and object controlling the array of axes on the
% figure
h_fig                   = ctl.figs.a4figure();
plot_array              = PlotArray(4, 4);

% cell arrays which will contain objects controlling the main
% tuning curve plot and the shuffled tuning curve histogram info
tuning_curve_plot       = {};
h_ax                    = {};

% loop over trial group types
for kk = 1 : length(trial_group_labels)

    % find where the axis should go and create it
    pos         = plot_array.get_position(2*kk-1);
    h_ax{kk}    = axes('units', 'centimeters', 'position', pos);
    hold on;
    
    if restrict_to_V_VT_clusters
        if isempty(trial_type_tuning_store{kk})
            continue
        end
    end
    
    % get mean and std across clusters for each speed bin (these don't
    % necessarily correspond to the same underlying speed across different
    % animals)
    avg_fr = mean(trial_type_tuning_store{kk}, 2);
    sd_fr = std(trial_type_tuning_store{kk}, [], 2);
    n = sum(~isnan(trial_type_tuning_store{kk}), 2);
    
    % plot the errorbar
    h_ebar = errorbar(h_ax{kk}, avg_fr, sd_fr./sqrt(n));
    
    h_ebar.Marker = 'o';
    h_ebar.MarkerSize = 3;
    h_ebar.Color = 'k';
    h_ebar.CapSize = 0;
    
    % mean and std of firing rate across clusters for stationary data
    avg_stat_fr = mean(stationary_store{kk}, 2);
    sd_stat_fr = std(stationary_store{kk}, [], 2);
    n_stat = sum(~isnan(stationary_store{kk}), 2);
    
    % plot the stationary data
    stat_col = [0.5, 0.5, 0.5];
    scatter(h_ax{kk}, 0, avg_stat_fr, 5, stat_col);
    line(h_ax{kk}, [0, 0], [-1, 1]*sd_stat_fr/n_stat, 'color', stat_col)
    
    % display the trial type as title
    title(gca, trial_group_labels{kk}, 'interpreter', 'none', ...
            'fontsize', 8);
    
    % axis labels on last axis
    if kk == length(trial_group_labels)
        xlabel('Bin #');
        ylabel('Firing rate (Hz)');
    end
    
    % print the number of clusters
    n_bins = size(trial_type_tuning_store{kk}, 1);
    n_clusters = size(trial_type_tuning_store{kk}, 2);
    text(h_ax{kk}, n_bins, 0, sprintf('n = %i', n_clusters), ...
        'horizontalalignment', 'right', 'verticalalignment', 'bottom');
    
end

% loop over all axes and find the min X, max X and max Y values
mx              = inf;
Mx              = -inf;
My              = -inf;

for kk = 1 : length(trial_group_labels)
    mx = min(h_ax{kk}.XLim(1), mx);
    Mx = max(h_ax{kk}.XLim(2), Mx);
    My = max(h_ax{kk}.YLim(2), My);
end

% again loop over axes and set the X and Y limits to the same value
for kk = 1 : length(trial_group_labels)
    h_ax{kk}.XLim = [-1, Mx];
    h_ax{kk}.YLim = [0, My];
    h_ax{kk}.Box = 'off';
end

% figure title
FigureTitle(h_fig, 'Averaged tuning curves');

% save to a new .pdf
if restrict_to_V_VT_clusters
    ctl.figs.save_fig('average_tuning_each_condition_restrict_clusters.pdf');
else
    ctl.figs.save_fig('average_tuning_each_condition.pdf');
end

