clear all
close all

% heatmaps showing transition from baseline to motion for all conditions

experiment      = 'mismatch_nov20+visual_flow';

baseline_t      = [-0.4, 0];
response_t      = [0, 0.4];
heatmap_padding = [-1, 1];


%%
config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, 'motion_onset');

if strcmp(experiment, 'mismatch_nov20+visual_flow')
    probe_fnames            = experiment_details('visual_flow', 'protocols');
    probe_fnames            = [probe_fnames, experiment_details('mismatch_nov20', 'protocols')];
else
    probe_fnames            = experiment_details(experiment, 'protocols');
end

if strcmp(experiment, 'visual_flow')
    protocol_ids        = [1, 2, 7, 8];
    protocol_labels     = VisualFlowExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'darkness')
    protocol_ids        = 1:3;
    protocol_labels     = DarknessExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'mismatch_nov20')
    protocol_ids        = [2, 4];
    protocol_labels     = MismatchExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'passive')
    protocol_ids        = 1:3;
    protocol_labels     = PassiveExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'head_tilt')
    protocol_ids        = 1:3;
    protocol_labels     = HeadTiltExperiment.protocol_label(protocol_ids);
elseif strcmp(experiment, 'mismatch_nov20+visual_flow')
    protocol_labels     = {'MVT', 'MV'};
end


n_narrow = [];
n_wide = [];
n_all = [];

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    narrow_clusters            = data.VISp_clusters([], 'narrow');
    wide_clusters             = data.VISp_clusters([], 'wide');
    all_clusters        = data.VISp_clusters([], 'any');
    
    n_narrow(probe_i) = length(narrow_clusters);
    n_wide(probe_i) = length(wide_clusters);
    n_all(probe_i) = length(all_clusters);
end

%%
figure, plot(n_narrow./n_all, 'o')
set(gca, 'xtick', 1:length(n_narrow), 'xticklabel', probe_fnames, 'ticklabelinterpreter', 'none')
for i = 1 : length(n_narrow)
    text(i, n_narrow(i)/n_all(i), sprintf('%i/%i', n_narrow(i), n_all(i)), 'rotation', 90);
end
box off
    
    
    