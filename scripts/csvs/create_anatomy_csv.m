% create a csv with anatomy information

experiment_groups   = {'visual_flow', 'darkness', 'mismatch_nov20', 'mismatch_jul21', 'mismatch_darkness_oct21'};

ctl                 = RC2Analysis();
probe_ids           = ctl.get_probe_ids(experiment_groups{:});
tbl                 = table([]);

r                   = 0;

for ii = 1 : length(probe_ids)
    
    data            = ctl.load_formatted_data(probe_ids{ii});
    
    shank_ids       = ctl.get_shank_ids(probe_ids{ii});
    
    for jj = 1 : length(shank_ids)
        
        anatomy         = data.anatomy{jj};
        
        r = r + 1;
        
        tbl.probe_id{r}                     = probe_ids{ii};
        tbl.shank_id(r)                     = shank_ids(jj);
        tbl.VISp1(r)                        = anatomy.VISp_boundaries_from_pia.lower(1);
        tbl.VISp23(r)                       = anatomy.VISp_boundaries_from_pia.lower(2);
        tbl.VISp4(r)                        = anatomy.VISp_boundaries_from_pia.lower(3);
        tbl.VISp5(r)                        = anatomy.VISp_boundaries_from_pia.lower(4);
        tbl.VISp6a(r)                       = anatomy.VISp_boundaries_from_pia.lower(5);
        tbl.VISp6b(r)                       = anatomy.VISp_boundaries_from_pia.lower(6);
        tbl.pia_from_tip(r)                 = anatomy.VISp_boundaries.upper(1);
        tbl.VISp1_from_tip(r)               = anatomy.VISp_boundaries.lower(1);
        tbl.VISp23_from_tip(r)              = anatomy.VISp_boundaries.lower(2);
        tbl.VISp4_from_tip(r)               = anatomy.VISp_boundaries.lower(3);
        tbl.VISp5_from_tip(r)               = anatomy.VISp_boundaries.lower(4);
        tbl.VISp6a_from_tip(r)              = anatomy.VISp_boundaries.lower(5);
        tbl.VISp6b_from_tip(r)              = anatomy.VISp_boundaries.lower(6);
    end
end

tbl(:, 1) = [];
writetable(tbl, 'anatomy_info.csv');