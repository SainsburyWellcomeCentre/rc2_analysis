function p = mm_do_ANOVA(g1, g2)

group_1 = [ones(numel(g1), 1); 2*ones(numel(g2), 1)];
sg_1 = repmat((1:size(g1, 1))', 1, size(g1, 2));
sg_2 = repmat((1:size(g2, 1))', 1, size(g2, 2));
group_2 = [sg_1(:); sg_2(:)];

p = anovan([g1(:); g2(:)], {group_1, group_2}, 'display', 'off');
