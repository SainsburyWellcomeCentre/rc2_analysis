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


n_fs = [];
n_rs = [];
n_all = [];

for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    fs_clusters            = data.VISp_clusters([], 'FS');
    rs_clusters             = data.VISp_clusters([], 'RS');
    all_clusters        = data.VISp_clusters([], 'any');
    
    n_fs(probe_i) = length(fs_clusters);
    n_rs(probe_i) = length(rs_clusters);
    n_all(probe_i) = length(all_clusters);
end

%%
figure, plot(n_fs./n_all, 'o')
set(gca, 'xtick', 1:length(n_fs), 'xticklabel', probe_fnames, 'ticklabelinterpreter', 'none')
for i = 1 : length(n_fs)
    text(i, n_fs(i)/n_all(i), sprintf('%i/%i', n_fs(i), n_all(i)), 'rotation', 90);
end
box off
    
    
    