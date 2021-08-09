% compare mismatch response to MV response

mm_csv_fname = 'C:\Users\Lee\Documents\mvelez\figures\mismatch_nov20\mm_ANOVA\mm_ANOVA_unity_plot.csv';
motion_csv_fname = 'C:\Users\Lee\Documents\mvelez\figures\mismatch_nov20\population_unity_plots\any\population_motion_vs_stationary.csv';

%%
mm_prot_id = 2;   % MVT gain up repsonse
motion_prot_id = 4;  % MV response

mm_table = readtable(mm_csv_fname);
motion_table = readtable(motion_csv_fname);

mm_table = mm_table(mm_table.protocol_id == mm_prot_id, :);
motion_table = motion_table(motion_table.protocol_id == motion_prot_id, :);

n_clusters = size(mm_table, 1);
assert(isequal(mm_table.cluster_id, motion_table.cluster_id));

n_mm_responsive = sum(mm_table.p_val_anova < 0.05);
n_motion_response = sum(motion_table.p_val_signrank < 0.05);
n_mm_and_motion_responsive = sum(mm_table.p_val_anova < 0.05 & motion_table.p_val_signrank < 0.05);
n_mm_not_motion_responsive = n_mm_responsive - n_mm_and_motion_responsive;


fprintf('# mm responsive clusters (gain up): %i/%i (%.1f%%)\n', n_mm_responsive, n_clusters, 100*n_mm_responsive/n_clusters);
fprintf('# motion responsive clusters (protocol %i): %i/%i (%.1f%%)\n', motion_prot_id, n_motion_response, n_clusters, 100*n_motion_response/n_clusters);
fprintf('# mm and motion responsive clusters: %i/%i (%.1f%%)\n', n_mm_and_motion_responsive, n_clusters, 100*n_mm_and_motion_responsive/n_clusters);
fprintf('# mm NOT motion responsive clusters: %i/%i (%.1f%%)\n', n_mm_not_motion_responsive, n_clusters, 100*n_mm_not_motion_responsive/n_clusters);

fprintf('cluster id, mm NOT motion responsive:\n');
mm_not_motion_idx = find(mm_table.p_val_anova < 0.05 & motion_table.p_val_signrank >= 0.05);
for i = 1 : length(mm_not_motion_idx)
    fprintf('    %i: %s, cluster %i\n', i, mm_table.probe_name{mm_not_motion_idx(i)}, mm_table.cluster_id(mm_not_motion_idx(i)));
end
