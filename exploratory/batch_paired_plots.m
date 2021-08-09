experiments = {'darkness', 'visual_flow'};
spiking_types = {'any', 'RS', 'FS'};

for i = 1 : length(experiments)
    for j = 1 : length(spiking_types)
        population_paired(experiments{i}, spiking_types{j});
        population_all_vs_all_paired(experiments{i}, spiking_types{j});
        population_mi_v_depth(experiments{i}, spiking_types{j});
        population_all_vs_all_mi_v_depth(experiments{i}, spiking_types{j});
    end
end
