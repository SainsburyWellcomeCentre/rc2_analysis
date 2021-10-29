%% all experiments used
probe_ids = {'CAA-1110262_rec1_rec2_rec3';
             'CAA-1110263_restricted_rec1_rec2_rec3';
             'CAA-1110264_rec1_rec2';
             'CAA-1110265_restricted_rec1_rec2_rec3';
             'CAA-1112222_rec1_rec2_rec3';
             'CAA-1112223_rec1_rec2_rec3';
             'CAA-1112224_rec1_rec2_rec3';
             'CAA-1112416_rec1_rec2_rec3';
             'CAA-1112417_rec1_rec2_rec3';
             'CA_176_1_rec1_rec2_rec3';
             'CA_176_3_rec1_rec2_rec3';
             'CAA-1112872_rec1_rec1b_rec2_rec3';
             'CAA-1112874_rec1_rec2_rec3';
             'CAA-1113219_rec1_rec2_rec3';
             'CAA-1113220_rec1_rec2_rec3';
             'CAA-1113222_rec1_rec2_rec3';
             'CAA-1114977_rec1_rec2_rec3';
             'CAA-1114978_rec1_rec2';
             'CAA-1114979_rec1_rec2';
             'CAA-1114980_rec1_rec2'};

batches_to_use = {1;
                  1:13;
                  1:7;
                  1:12;
                  1:10;
                  1;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10;
                  1:10};
         
         
ctl                         = RC2Analysis();

new_l5                      = nan(length(probe_ids), 1);
new_l5_consistent_batches   = nan(length(probe_ids), 1);
l5_legacy                   = nan(length(probe_ids), 1);
l5_legacy_consistent_batches = nan(length(probe_ids), 1);

for ii = 1 : length(probe_ids)
    
     rec = ctl.load.spikeglx_ap_recording(probe_ids{ii});
    
    hf_power = HighFrequencyPowerProfile(rec, 0);
    
    % new method
    hf_power.compute_power_on_all_channels();
    hf_power.interpolate_across_channels();
    hf_power.smooth_across_channels();
    
    new_l5_consistent_batches(ii) = hf_power.search_for_l5;
    
    hf_power.batches_to_use = batches_to_use{ii};
    
    hf_power.interpolate_across_channels();
    hf_power.smooth_across_channels();
    new_l5(ii) = hf_power.search_for_l5;
    
    
    % switch to legacy mode
%     hf_power = HighFrequencyPower(glx_dir, 0, 'lf');
%     hf_power.compute_legacy();
%     
%     l5_legacy_consistent_batches(ii) = hf_power.search_for_l5;
%     hf_power.batches_to_use = batches_to_use{ii};
%     l5_legacy(ii) = hf_power.search_for_l5;
end