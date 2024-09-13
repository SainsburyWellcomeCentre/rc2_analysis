% Plot acceleration tuning curves for all cortical (VISp) clusters
% in the selected experiment_groups and conditions (trial_group_labels)

% Using the acceleration table with 20 bins calculated without splitting acceleration and deceleration


% Initialize parameters
experiment_groups       = {'darkness', 'mismatch_darkness_oct21'};
trial_group_labels      = {{'T_bank', 'T_RT', 'T_R'}, 'T'};

restricted              = false;
save_figs               = true;
overwrite               = true;
figure_dir              = {'all_tuning_curves'};

% main controller object
ctl                     = RC2Analysis();

% setup directory to save figures in if 'save_figs' is true
ctl.setup_figures(figure_dir, save_figs);

% get an average tuning curve for each trial type
tuning_store = [];
stationary_store = [];

for ll = 1 : length(experiment_groups)
    % get probe ids
    probe_ids           = ctl.get_probe_ids(experiment_groups{ll});
    % loop over experiments
    for ii = 1 : length(probe_ids)  

        data        = ctl.load_formatted_data(probe_ids{ii});
        clusters    = data.VISp_clusters();

        % loop over VISp clusters in that experiment
        for jj = 1 : length(clusters)
            %load tuning
            tuning      = {};
            % Trial group labels depend on experimental group
            tuning      = data.load_tuning_curves_acceleration(clusters(jj).id, trial_group_labels{ll});
    
            % average tuning curve for this cluster
            avg_tuning_all = nanmean(tuning{1, 1}.tuning, 2); % bins of acceleration and deceleration calculated together
            avg_stationary = mean(tuning{1, 1}.stationary_fr);

            tuning_store = [tuning_store, avg_tuning_all];
            stationary_store = [stationary_store, avg_stationary];
        end
    end
end

%%

% create A4 figure and object controlling the array of axes on the
% figure
h_fig                   = ctl.figs.a4figure();
plot_array              = PlotArray(4, 4);

% cell arrays which will contain objects controlling the main
% tuning curve plot and the shuffled tuning curve histogram info
tuning_curve_plot       = {};
h_ax                    = {};

% loop over trial group types
% find where the axis should go and create it
pos         = plot_array.get_position(2);
h_ax    = axes('units', 'centimeters', 'position', pos);
hold on;


% get mean and std across clusters for each speed bin (these don't
% necessarily correspond to the same underlying speed across different
% animals)
avg_fr = mean(tuning_store, 2);
sd_fr = std(tuning_store, [], 2);
n = sum(~isnan(tuning_store), 2);

% plot the errorbar
h_ebar = errorbar(h_ax, avg_fr, sd_fr./sqrt(n));

h_ebar.Marker = 'o';
h_ebar.MarkerSize = 3;
h_ebar.Color = 'k';
h_ebar.CapSize = 0;

% mean and std of firing rate across clusters for stationary data
avg_stat_fr = mean(stationary_store, 2);
sd_stat_fr = std(stationary_store, [], 2);
n_stat = sum(~isnan(stationary_store), 2);

% plot the stationary data
stat_col = [0.5, 0.5, 0.5];
scatter(h_ax, 0, avg_stat_fr, 5, stat_col);
line(h_ax, [0, 0], [-1, 1]*sd_stat_fr/n_stat, 'color', stat_col)

% print the number of clusters
n_bins = size(tuning_store, 1);
n_clusters = size(tuning_store, 2);
text(h_ax, n_bins, 0, sprintf('n = %i', n_clusters), ...
    'horizontalalignment', 'right', 'verticalalignment', 'bottom');
    
%find the min X, max X and max Y values
mx              = inf;
Mx              = -inf;
My              = -inf;

mx = min(h_ax.XLim(1), mx);
Mx = max(h_ax.XLim(2), Mx);
My = max(h_ax.YLim(2), My);

% set the X and Y limits to the same value
h_ax.XLim = [-1, Mx];
h_ax.YLim = [0, My];
h_ax.Box = 'off';

% figure title
FigureTitle(h_fig, 'Averaged tuning curves');

% save to a new .pdf
if restricted
    ctl.figs.save_fig('all_acceleration_average_tuning_each_condition_restrict_clusters.pdf', true);
else
    ctl.figs.save_fig('all_acceleration_average_tuning_each_condition.pdf', true);
end

%% Anova
[p,tbl,stats] = anova1(tuning_store.')
