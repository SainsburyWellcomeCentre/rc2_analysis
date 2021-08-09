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
figs.set_figure_subdir('visual_flow', 'all_trials_with_cluster_raster');

% get details of the experiment
probe_fnames = experiment_details(experiment, combination);

trial_types = {'Coupled', 'EncoderOnly'};
fname_suffix = {'MVT', 'MV'};

for probe_i = 1 : length(probe_fnames)
    
    % all data for this probe recording
    data                = load_formatted_data(probe_fnames{probe_i}, config);
    
    % get visual flow experiment object
    vf = VisualFlowExperiment(data, config);
    
    % get selected clusters in primary visual cortex
    clusters = data.VISp_clusters;
    
    base_t = {};
    M_trace = {};
    V_trace = {};
    T_trace = {};
    trial_id = {};
    spike_times = {};
    spike_conv = {};
    
    base_t_rep = {};
    M_trace_rep = {};
    V_trace_rep = {};
    T_trace_rep = {};
    trial_id_rep = {};
    bouts_rep = {};
    protocol_type_rep = {};
    spike_times_rep = {};
    spike_conv_rep = {};
    
    for type_i = 1 : length(trial_types)
        
        these_trials = vf.trials_of_type(trial_types{type_i});
        
        for trial_i = 1 : length(these_trials)
            
            base_t{type_i}{trial_i} = these_trials(trial_i).probe_t - ...
                these_trials(trial_i).probe_t(1);
            M_trace{type_i}{trial_i} = these_trials(trial_i).filtered_teensy;
            V_trace{type_i}{trial_i} = these_trials(trial_i).multiplexer_output;
            T_trace{type_i}{trial_i} = these_trials(trial_i).stage;
            trial_id{type_i}{trial_i} = these_trials(trial_i).id;
            
            % start and end times of this trial in probe time
            start_time = these_trials(trial_i).probe_t(1);
            end_time = these_trials(trial_i).probe_t(end);
            
            % for each cluster get spike times during this trial
            for cluster_i = 1 : length(clusters)
                
                % which spikes occur on this trial
                valid_spikes_mask = clusters(cluster_i).spike_times > start_time & ...
                    clusters(cluster_i).spike_times < end_time;
                
                % store the spike times, relative to base_t
                spike_times{type_i}{trial_i}{cluster_i} = ...
                    clusters(cluster_i).spike_times(valid_spikes_mask) - start_time;
                
                % get convolved spike rate
                fr = FiringRate(clusters(cluster_i).spike_times);
                spike_conv{type_i}{trial_i}{cluster_i} = ...
                    fr.get_convolution(these_trials(trial_i).probe_t);
            end
            
            
            % get repeats of trials
            repeats = vf.get_replay_of(these_trials(trial_i));
            
            % match
            for rep_i = 1 : length(repeats)
                
                %
                offset = vf.get_offset(these_trials(trial_i), repeats(rep_i));
                
                n_remaining = length(repeats(rep_i).probe_t) - offset;
                n_to_plot = min(n_remaining, length(these_trials(trial_i).probe_t));
                
                base_t_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).probe_t(offset+(0:n_to_plot-1)) - ...
                    repeats(rep_i).probe_t(offset);
                M_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).filtered_teensy(offset+(0:n_to_plot-1));
                V_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).multiplexer_output(offset+(0:n_to_plot-1));
                T_trace_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).stage(offset+(0:n_to_plot-1));
                trial_id_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).id;
                protocol_type_rep{type_i}{trial_i}{rep_i} = repeats(rep_i).protocol;
                
                % start and end times of this trial in probe time
                start_time = repeats(rep_i).probe_t(offset);
                end_time = repeats(rep_i).probe_t(offset+n_to_plot-1);
                
                % for each cluster get spike times during this trial
                for cluster_i = 1 : length(clusters)
                    
                    % which spikes occur on this trial
                    valid_spikes_mask = clusters(cluster_i).spike_times > start_time & ...
                        clusters(cluster_i).spike_times < end_time;
                    
                    % store the spike times, relative to base_t_rep
                    spike_times_rep{type_i}{trial_i}{rep_i}{cluster_i} = ...
                        clusters(cluster_i).spike_times(valid_spikes_mask) - start_time;
                    
                    % get convolved spike rate
                    fr = FiringRate(clusters(cluster_i).spike_times);
                    spike_conv_rep{type_i}{trial_i}{rep_i}{cluster_i} = ...
                        fr.get_convolution(repeats(rep_i).probe_t(offset+(0:n_to_plot-1)));
                    
                end
                
            end
            
        end
        
    end
    
    all_traces = {M_trace, V_trace, T_trace};
    all_traces_rep = {M_trace_rep, V_trace_rep, T_trace_rep};
    axis_labels = {'M', 'V', 'T'};
    
    %% plot each trial
    for type_i = 1 : length(base_t)
        
        
        for cluster_i = 1 : length(clusters)
            
            for trial_i = 1 : length(base_t{type_i})
                
                section_n = 1;
                sp_n = (ceil(section_n/2) - 1)*8 + mod(section_n-1, 2) + 1;
                
                if section_n == 1
                    h_fig = figs.a4figure();
                end
                
                h_ax = [];
                h_ax_rate = [];
                to_scatter = [];
                
                yL = 0;
                yL_rate = 0;
                xL = [inf, -inf];
                
                for trace_i = 1 : length(all_traces)
                    
                    h_ax(end+1) = subplot(8, 2, sp_n + (trace_i-1)*2);
                    
                    trace_plot(h_ax(end), base_t{type_i}{trial_i}, ...
                        all_traces{trace_i}{type_i}{trial_i}, ...
                        [], []);
                    
                    if trace_i == 1
                        title(sprintf('Trial %i, %s', ...
                            trial_id{type_i}{trial_i}, ...
                            fname_suffix{type_i}));
                    end
                    
                    
                    if trace_i == 3
                        
                        n_spikes = length(spike_times{type_i}{trial_i}{cluster_i});
                        h_ax_rate(end+1) = subplot(8, 2, sp_n + trace_i*2);
                        hold on;
                        plot(base_t{type_i}{trial_i}, spike_conv{type_i}{trial_i}{cluster_i}, 'color', [0.5, 0.5, 0.5]);
                        yL_rate = max(yL_rate, max(get(h_ax_rate(end), 'ylim')));
                        to_scatter{end+1} = spike_times{type_i}{trial_i}{cluster_i};
                        set(gca, 'plotboxaspectratio', [3, 1, 1]);
                        box off;
                        
                    end
                    
                    yL = max(yL, max(get(h_ax(end), 'ylim')));
                    xL(1) = min(xL(1), min(get(h_ax(end), 'xlim')));
                    xL(2) = max(xL(2), max(get(h_ax(end), 'xlim')));
                    
                end
                
                for rep_i = 1 : length(base_t_rep{type_i}{trial_i})
                    
                    if strcmp(protocol_type_rep{type_i}{trial_i}{rep_i}, 'StageOnly')
                        
                        section_n = 2;
                        title_str = 'VT';
                        
                    else
                        
                        section_n = 3;
                        title_str = 'V';
                        
                    end
                    
                    sp_n = (ceil(section_n/2) - 1)*8 + mod(section_n-1, 2) + 1;
                    
                    for trace_i = 1 : length(all_traces_rep)
                        
                        h_ax(end+1) = subplot(8, 2, sp_n + (trace_i - 1)*2);
                        
                        trace_plot(h_ax(end), base_t_rep{type_i}{trial_i}{rep_i}, ...
                            all_traces_rep{trace_i}{type_i}{trial_i}{rep_i}, ...
                            [], []);
                        
                        if trace_i == 1
                            title(sprintf('Trial %i, %s (%s)', ...
                                trial_id_rep{type_i}{trial_i}{rep_i}, ...
                                title_str, ...
                                fname_suffix{type_i}));
                        end
                        
                        if trace_i == 3
                            
                            n_spikes = length(spike_times_rep{type_i}{trial_i}{rep_i}{cluster_i});
                            
                            h_ax_rate(end+1) = subplot(8, 2, sp_n + trace_i*2);
                            hold on;
                            plot(base_t_rep{type_i}{trial_i}{rep_i}, spike_conv_rep{type_i}{trial_i}{rep_i}{cluster_i}, 'color', [0.5, 0.5, 0.5]);
                            % collect for scattering
                            to_scatter{end+1} = spike_times_rep{type_i}{trial_i}{rep_i}{cluster_i};
                            yL_rate = max(yL_rate, max(get(h_ax_rate(end), 'ylim')));
                            set(gca, 'plotboxaspectratio', [3, 1, 1]);
                            box off;
                        end
                        
                        yL = max(yL, max(get(h_ax(end), 'ylim')));
                        xL(1) = min(xL(1), min(get(h_ax(end), 'xlim')));
                        xL(2) = max(xL(2), max(get(h_ax(end), 'xlim')));
                        
                    end
                    
                end
                
                for ax_i = 1 : length(h_ax)
                    
                    set(gcf, 'currentaxes', h_ax(ax_i));
                    
                    set(h_ax(ax_i), 'xlim', xL);
                    
                    set(h_ax(ax_i), 'ylim', [-5, yL]);
                    
                    if ax_i == 1
                        ylabel('(cm/s)');
                    end
                    
                    text(xL(1), yL, ...
                        axis_labels{mod(ax_i-1, 3)+1}, ...
                        'fontweight', 'bold', ...
                        'fontsize', 14, ...
                        'horizontalalignment', 'left', ...
                        'verticalalignment', 'top');
                end
                
                for ax_i = 1 : length(h_ax_rate)
                    
                    set(gcf, 'currentaxes', h_ax_rate(ax_i));
                    set(h_ax_rate(ax_i), 'xlim', xL);
                    set(h_ax_rate(ax_i), 'ylim', [-5, yL_rate]);
                    scatter(to_scatter{ax_i}, yL_rate * ones(length(to_scatter{ax_i}), 1), 10, 'k', 'fill')
                    if ax_i == 1
                        ylabel('Convolved rate (Hz)');
                    end
                end
                
                
                FigureTitle(h_fig, sprintf('%s, %s, Cluster %i, %s, Trial %i', ...
                    probe_fnames{probe_i}, ...
                    trial_types{type_i}, ...
                    clusters(cluster_i).id, ...
                    clusters(cluster_i).region_str, ...
                    trial_id{type_i}{trial_i}));
                
                figs.save_fig_to_join();
                
            end  % for trial
            
            figs.join_figs( ...
                sprintf('%s_cluster_%03i_%s_all_trials_with_cluster_raster.pdf', ...
                        probe_fnames{probe_i}, ...
                        clusters(cluster_i).id, ...
                        fname_suffix{type_i}));
            figs.clear_figs();
            
        end  % for cluster
        
    end
    
end






