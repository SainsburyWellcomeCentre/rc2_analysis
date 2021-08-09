% plot overlap of running significant and translation signficant responses
experiment              = 'visual_flow';

config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = true;
figs.set_figure_subdir(experiment, 'overlap_plots');

probe_fnames            = experiment_details(experiment, 'protocols');

if strcmp(experiment, 'visual_flow')
    protocols           = 1:6;
    protocol_labels     = VisualFlowExperiment.protocol_label(protocols);
elseif strcmp(experiment, 'darkness')
    protocols           = 1:3;
    protocol_labels     = DarknessExperiment.protocol_label(protocols);
end


%% test each protocol motion against stationary

x_all_stat                   = cell(length(protocols), 1);
y_all_mot                   = cell(length(protocols), 1);
p_all_stat_v_mot                   = cell(length(protocols), 1);

for probe_i = 1 : length(probe_fnames)
    
%     data                = load_formatted_data(probe_fnames{probe_i}, config);
    this_data           = get_data_for_recording_id(data, probe_fnames{probe_i});
    
    clusters            = this_data.VISp_clusters;
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(this_data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(this_data, config);
    end
    
    for prot_i = 1 : length(protocols)
        
        for cluster_i = 1 : length(clusters)
            
            x           = exp_obj.trial_stationary_fr(clusters(cluster_i).id, protocols(prot_i));
            y           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_i));
            
            x_all_stat{prot_i}(end+1) = nanmedian(x);
            y_all_mot{prot_i}(end+1) = nanmedian(y);
            p_all_stat_v_mot{prot_i}(end+1) = signrank(x, y);
            
            if p_all_stat_v_mot{prot_i}(end) < 0.05 & x_all_stat{prot_i}(end) == y_all_mot{prot_i}(end)
                probe_i, clusters(cluster_i).id
            end
        end
    end
end

% overlap signficantly
if strcmp(experiment, 'visual_flow')
    sig_and_MV_VT_overlap = p_all_stat_v_mot{2} < 0.05 & p_all_stat_v_mot{4} < 0.05;
%     sig_and_MV_VT_overlap = p_all{2, 6} < 0.05 & p_all{4, 6} < 0.05;
elseif strcmp(experiment, 'darkness') 
    sig_and_M_T_overlap = p_all_stat_v_mot{2} < 0.05 & p_all_stat_v_mot{3} < 0.05;
end



%% now test each protocol against each other protocol
if strcmp(experiment, 'visual_flow')
    protocols           = 1:6;
    protocol_labels     = VisualFlowExperiment.protocol_label(protocols);
elseif strcmp(experiment, 'darkness')
    protocols           = 1:2;
    protocol_labels     = DarknessExperiment.protocol_label(protocols);
end

x_all                   = cell(length(protocols));
y_all                   = cell(length(protocols));
is_increase             = cell(length(protocols));
p_all                   = cell(length(protocols));

for probe_i = 1 : length(probe_fnames)
    
    this_data           = get_data_for_recording_id(data, probe_fnames{probe_i});
    clusters            = this_data.VISp_clusters;
    
    if strcmp(experiment, 'visual_flow')
        exp_obj         = VisualFlowExperiment(this_data, config);
    elseif strcmp(experiment, 'darkness')
        exp_obj         = DarknessExperiment(this_data, config);
    end
    
    for prot_y = 1 : length(protocols)-1
        for prot_x = prot_y+1 : length(protocols)
            
            for cluster_i = 1 : length(clusters)
                
                x           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_x));
                y           = exp_obj.trial_motion_fr(clusters(cluster_i).id, protocols(prot_y));
                
                [x_all{prot_y, prot_x}(end+1), ...
                    y_all{prot_y, prot_x}(end+1), ...
                    p_all{prot_y, prot_x}(end+1), ...
                    is_increase{prot_y, prot_x}(end+1)] = compare_groups_with_signrank(x, y);
                
            end
        end
    end
end



%% venn diagram
if strcmp(experiment, 'visual_flow')
    
%     A = p_all{2, 6} < 0.05;
%     B = p_all{4, 6} < 0.05;
    
    A = p_all_stat_v_mot{2} < 0.05;
    B = p_all_stat_v_mot{4} < 0.05;
    
    figure
    h_ax = axes();
    s = SimpleVenn(h_ax);
    s.A = sum(A);
    s.B = sum(B);
    s.A_and_B = sum(A & B);
    s.N = length(A);
    s.name_A = 'MV-V';
    s.name_B = 'VT-V';
    s.plot();
    
    print('venn_MV-V_and_VT-V.pdf', '-bestfit', '-dpdf');
    
elseif strcmp(experiment, 'darkness')
    
    A = p_all_stat_v_mot{2} < 0.05;
    B = p_all_stat_v_mot{3} < 0.05;
    
    figure
    h_ax = axes();
    s = SimpleVenn(h_ax);
    s.A = sum(A);
    s.B = sum(B);
    s.A_and_B = sum(A & B);
    s.N = length(A);
    s.name_A = 'M';
    s.name_B = 'T';
    s.plot();
    
    print('venn_M_and_T.pdf', '-bestfit', '-dpdf');
end


%% plot all against all, highlight the M/T overlap
h_fig                   = figs.a4figure();
plot_array             = PlotArray(6, 6);
u                       = UnityPlotPopulation.empty();

for prot_y = 1 : length(protocols)-1
    for prot_x = prot_y+1 : length(protocols)
        
        sp_idx      = (prot_y-1)*length(protocols) + prot_x;
        pos         = plot_array.get_position(sp_idx);
        h_ax        = axes('units', 'centimeters', 'position', pos);
        
        u(end+1)   = UnityPlotPopulation(x_all{prot_y, prot_x}, ...
            y_all{prot_y, prot_x}, ...
            p_all{prot_y, prot_x}, ...
            is_increase{prot_y, prot_x}, ...
            h_ax);
        
        u(end).plot();
        
        % highlight the ones with significant overlap
        if strcmp(experiment, 'visual_flow')
            u(end).highlight(sig_and_MV_VT_overlap);
        elseif strcmp(experiment, 'darkness')
            u(end).highlight(sig_and_M_T_overlap);
        end
        
        u(end).xlabel(protocol_labels{prot_x});
        u(end).ylabel(protocol_labels{prot_y});
        u(end).add_histogram(1);
    end
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for prot_i = 1 : length(u)
    u(prot_i).xlim([m, M])
end

FigureTitle(h_fig, 'population, all vs. all');
% figs.save_fig('population_all_vs_all.pdf');