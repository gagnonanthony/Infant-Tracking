#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    GENERATE_MASKS;
    LOCAL_TRACKING
} from '../processes/tracking.nf'

workflow TRACKING {
    take:
        brain_wm_mask_channel
        fodf_channel
        fa_channel
    main:
        masks_channel = brain_wm_mask_channel
            .combine(fa_channel, by: 0)
        
        GENERATE_MASKS(masks_channel)

        tracking_channel = fodf_channel
            .combine(GENERATE_MASKS.out.masks, by: 0)
        
        LOCAL_TRACKING(tracking_channel)
}