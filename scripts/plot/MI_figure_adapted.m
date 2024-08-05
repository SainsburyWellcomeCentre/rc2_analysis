% Make a modulation index (MI) plot containing cortex layer information and a histogram.
% REQUIRES adding mvelez_ms_figures to PATH

% Initialize controller and variables
ctl                         = RC2Analysis();

experiment_groups           = {{'visual_flow'}};
trial_group_labels_1        = {{'V_RVT', 'V_RV'}};
trial_group_labels_2       = {{'VT_RVT', 'VT_RV'}};

cols                        = get_colours();

% Initialize objects for plotting 
modulation_index_full = [];
direction_full = [];
avg_anatomy_full = Anatomy.empty();
averaged_cortical_position_full = [];

% For every experiment group
for ii = 1 : length(experiment_groups)
    % Get the relevant data and calculate MI with a dedicated method
    probe_ids = ctl.get_probe_ids(experiment_groups{ii});
    
    [modulation_index, direction, avg_anatomy, averaged_cortical_position] = ...
        modulation_index_data([], experiment_groups{ii}, trial_group_labels_1{ii}, trial_group_labels_2{ii}, inf);

    % Fill in the objects required for plotting
    modulation_index_full = [modulation_index_full modulation_index];
    direction_full = [direction_full direction];
    
    if ii == 1
        avg_anatomy_full.VISp_layers = avg_anatomy.VISp_layers;
        avg_anatomy_full.anatomy_array = avg_anatomy.anatomy_array;
    else
        avg_anatomy_full.anatomy_array = [avg_anatomy_full.anatomy_array avg_anatomy.anatomy_array];
    end
    
    averaged_cortical_position_full = [averaged_cortical_position_full; averaged_cortical_position];
end

% Plot
figure(1);
h_ax = axes();
hold on
fmt.x_label = {'Modulation index', 'cond1 vs. cond2'};
modulation_index_plot(h_ax, modulation_index_full, direction_full, avg_anatomy_full, averaged_cortical_position_full, fmt);
