% combine clusters from two 

mm_motion_csv_fname = 'C:\Users\Lee\Documents\mvelez\figures\mismatch_nov20\population_unity_plots\any\population_motion_vs_stationary.csv';
visual_flow_motion_csv_fname = 'C:\Users\Lee\Documents\mvelez\figures\visual_flow\population_unity_plots\any\population_motion_vs_stationary.csv';

mm_motion_table = readtable(mm_motion_csv_fname);
visual_flow_motion_table = readtable(visual_flow_motion_csv_fname);



%%
config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir('mismatch_nov20+visual_flow', 'population_unity_plots', 'any');


%%  motion vs. stationary
h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);
u                       = UnityPlotPopulation.empty();

protocols               = {[1, 1], [3, 2]};
protocol_labels         = {'MVT', 'MV'};

for prot_i = 1 : length(protocols)
    
    idx_1 = mm_motion_table.protocol_id == protocols{prot_i}(1);
    idx_2 = visual_flow_motion_table.protocol_id == protocols{prot_i}(2);
    
    x_all = [mm_motion_table.stationary_fr(idx_1); visual_flow_motion_table.stationary_fr(idx_2)];
    y_all = [mm_motion_table.motion_fr(idx_1); visual_flow_motion_table.motion_fr(idx_2)];
    p_all = [mm_motion_table.p_val_signrank(idx_1); visual_flow_motion_table.p_val_signrank(idx_2)];
    is_increase = [mm_motion_table.is_increase(idx_1); visual_flow_motion_table.is_increase(idx_2)];
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    u(end+1)   = UnityPlotPopulation(x_all, ...
        y_all, ...
        p_all, ...
        is_increase, ...
        h_ax);
    
    u(end).plot();
    
    u(end).xlabel('Stationary (Hz)');
    u(end).ylabel('Motion (Hz)');
    
    u(end).title(protocol_labels{prot_i});
    u(end).add_histogram(1);
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for prot_i = 1 : length(u)
    u(prot_i).xlim([0, 25]);%[m, M])
end

FigureTitle(h_fig, 'population, motion vs. stationary');
figs.save_fig('population_motion_vs_stationary.pdf');



%%  all vs. all
mm_motion_csv_fname = 'C:\Users\Lee\Documents\mvelez\figures\mismatch_nov20\population_unity_plots\any\population_all_vs_all.csv';
visual_flow_motion_csv_fname = 'C:\Users\Lee\Documents\mvelez\figures\visual_flow\population_unity_plots\any\population_all_vs_all.csv';

mm_motion_table = readtable(mm_motion_csv_fname);
visual_flow_motion_table = readtable(visual_flow_motion_csv_fname);

h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);
u                       = UnityPlotPopulation.empty();

protocols_x             = [4, 2];
protocols_y             = [2, 1];
        
idx_1       = mm_motion_table.protocol_x_id == protocols_x(1) & mm_motion_table.protocol_y_id == protocols_y(1);
idx_2       = visual_flow_motion_table.protocol_x_id == protocols_x(2) & visual_flow_motion_table.protocol_y_id == protocols_y(2);

x_all       = [mm_motion_table.x_fr(idx_1); visual_flow_motion_table.x_fr(idx_2)];
y_all       = [mm_motion_table.y_fr(idx_1); visual_flow_motion_table.y_fr(idx_2)];
p_all       = [mm_motion_table.p_val_signrank(idx_1); visual_flow_motion_table.p_val_signrank(idx_2)];
is_increase = [mm_motion_table.is_increase(idx_1); visual_flow_motion_table.is_increase(idx_2)];

pos         = plot_array.get_position(1);
h_ax        = axes('units', 'centimeters', 'position', pos);

u   = UnityPlotPopulation(x_all, y_all, p_all, is_increase, h_ax);

u.plot();

u.xlabel('MV');
u.ylabel('MVT');
u.add_histogram(1);

m               = min([u(:).min]);
M               = max([u(:).max]);

for prot_i = 1 : length(u)
    u(prot_i).xlim([m, M])
end

FigureTitle(h_fig, 'population, all vs. all');
figs.save_fig('population_all_vs_all.pdf');