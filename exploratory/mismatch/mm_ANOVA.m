config                  = RC2AnalysisConfig();

figs                    = RC2Figures(config);
figs.save_on            = false;
figs.set_figure_subdir('mismatch_nov20', 'mm_ANOVA');

csvs                    = CSVManager(config);
csvs.save_on            = true;
csvs.set_csv_fulldir(figs.curr_dir);

probe_fnames            = experiment_details('mismatch_nov20', 'protocol');

protocols               = MismatchExperiment.protocol_ids;
protocol_labels         = MismatchExperiment.protocol_label;

cluster_id              = cell(length(protocols), 1);
probe_name              = cell(length(protocols), 1);
protocol_id             = cell(length(protocols), 1);

x_all                   = cell(length(protocols), 1);
y_all                   = cell(length(protocols), 1);
p_all                   = cell(length(protocols), 1);


for probe_i = 1 : length(probe_fnames)
    
    probe_i
    
%     data                = load_formatted_data(probe_fnames{probe_i}, config);
    this_data           = get_data_for_recording_id(data, probe_fnames{probe_i});
    mm                  = MismatchExperiment(this_data, config);
    clusters            = this_data.VISp_clusters;
    
    for cluster_i = 1 : length(clusters)
        
        for prot_i = 1 : length(protocols)
            
            [baseline, response, response_ctl] = mm.windowed_mm_responses(clusters(cluster_i), prot_i);
            
            probe_name{prot_i}{end+1, 1} = probe_fnames{probe_i};
            protocol_id{prot_i}(end+1, 1) = protocols(prot_i);
            cluster_id{prot_i}(end+1, 1) = clusters(cluster_i).id;
            
            x_all{prot_i}(end+1, 1) = nanmean(baseline(:));
            y_all{prot_i}(end+1, 1) = nanmean(response(:));
            
            p = mm_do_ANOVA(baseline', response');
            p_ctl = mm_do_ANOVA(baseline', response_ctl');
%             [~, pt] = ttest2(sum(baseline', 2), sum(response', 2));
            
            if p_ctl(1) < 0.05
                p_all{prot_i}(end+1, 1) = nan;
%                 pt_all{prot_i}(end+1) = nan;
            else
                p_all{prot_i}(end+1, 1) = p(1);
%                 pt_all{prot_i}(end+1) = pt;
            end
        end
    end
end

csvs.create_table(  'probe_name',   cat(1, probe_name{:}), ...
                    'cluster_id',   cat(1, cluster_id{:}), ...
                    'protocol_id',  cat(1, protocol_id{:}), ...
                    'baseline_fr',  cat(1, x_all{:}), ...
                    'response_fr',  cat(1, y_all{:}), ...
                    'p_val_anova',  cat(1, p_all{:}))
csvs.save('mm_ANOVA_unity_plot');


%%
h_fig                   = figs.a4figure();
plot_array              = PlotArray(3, 2);
u                       = UnityPlotPopulation.empty();

for prot_i = 1 : length(protocols)
    
    pos         = plot_array.get_position(prot_i);
    h_ax        = axes('units', 'centimeters', 'position', pos);
    
    u(end+1)   = UnityPlotPopulation(x_all{prot_i}, ...
        y_all{prot_i}, ...
        p_all{prot_i}, ...
        x_all{prot_i} < y_all{prot_i}, ...
        h_ax);
    
    if ismember(prot_i, [1, 3])
        u(end).marker_style = 'v';
    end
    
    u(end).plot();
    
    u(end).xlabel('Baseline FR (Hz)');
    u(end).ylabel('Response FR (Hz)');
    
    u(end).title(protocol_labels{prot_i});
    u(end).add_histogram(1);
end

m               = min([u(:).min]);
M               = max([u(:).max]);

for prot_i = 1 : length(u)
    u(prot_i).xlim([m, M])
end

FigureTitle(h_fig, '');
figs.save_fig('mm_ANOVA_unity_plot.pdf');