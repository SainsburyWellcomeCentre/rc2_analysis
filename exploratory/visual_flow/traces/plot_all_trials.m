clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'visual_flow';  % 'darkness', 'visual_flow', 'head_tilt', 'mismatch_nov20'
combination = 'protocols';

config = RC2AnalysisConfig();
figs = RC2Figures(config);

% where to save figure
figs.save_on = true;
figs.set_figure_subdir('visual_flow', 'all_trials');

% get details of the experiment
probe_fnames = experiment_details(experiment, combination);

trial_types = {'Coupled', 'EncoderOnly', 'ReplayOnly', 'StageOnly'};


for probe_i = 1 : length(probe_fnames)
    
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    vf = VisualFlowExperiment(data, config);
    
    base_t = {};
    M_trace = {};
    V_trace = {};
    T_trace = {};
    trial_id = {};
    
    for type_i = 1 : length(trial_types)
        
        these_trials = vf.trials_of_type(trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            
            base_t{type_i}{trial_i} = these_trials(trial_i).probe_t - ...
                these_trials(trial_i).probe_t(1);
            
            M_trace{type_i}{trial_i} = these_trials(trial_i).filtered_teensy;
            V_trace{type_i}{trial_i} = these_trials(trial_i).multiplexer_output;
            T_trace{type_i}{trial_i} = these_trials(trial_i).stage;
            trial_id{type_i}{trial_i} = these_trials(trial_i).id;
        end
    end
    
    
    %% plot each trial
    for type_i = 1 : length(base_t)
        
        
        
        for trial_i = 1 : length(base_t{type_i})
            
            section_n = mod(trial_i-1, 4) + 1;
            sp_n = (ceil(section_n/2) - 1)*8 + mod(section_n-1, 2) + 1;
            
            if section_n == 1
                h_fig = figs.a4figure();
            end
            
            yL = 0;
            h_ax(1) = subplot(8, 2, sp_n);
            plot(base_t{type_i}{trial_i}, M_trace{type_i}{trial_i});
            set(gca, 'plotboxaspectratio', [3, 1, 1]);
            box off;
            title(sprintf('Trial %i', trial_id{type_i}{trial_i}));
            yL = max(yL, max(get(gca, 'ylim')));
            h_ax(2) = subplot(8, 2, sp_n+2);
            plot(base_t{type_i}{trial_i}, V_trace{type_i}{trial_i});
            set(gca, 'plotboxaspectratio', [3, 1, 1])
            box off;
            yL = max(yL, max(get(gca, 'ylim')));
            h_ax(3) = subplot(8, 2, sp_n+4);
            plot(base_t{type_i}{trial_i}, T_trace{type_i}{trial_i});
            set(gca, 'plotboxaspectratio', [3, 1, 1])
            box off;
            yL = max(yL, max(get(gca, 'ylim')));
            
            for ax_i = 1 : 3
                set(h_ax(ax_i), 'ylim', [-5, yL]);
            end
            
            
            if section_n == 4 || trial_i == length(base_t{type_i})
                FigureTitle(h_fig, sprintf('%s, %s', probe_fnames{probe_i}, trial_types{type_i}));
                figs.save_fig_to_join();
            end
            
        end
        
        figs.join_figs(sprintf('%s_%s_all_trials.pdf', probe_fnames{probe_i}, trial_types{type_i}));
        figs.clear_figs();
    end
end






