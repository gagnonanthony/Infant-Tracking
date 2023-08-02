#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    BET_DWI;
    BET_T2;
    DENOISING;
    TOPUP;
    EDDY_TOPUP;
    N4;
    CROP_DWI;
    CROP_ANAT;
    RESAMPLE_ANAT;
    NORMALIZE;
    RESAMPLE_DWI;
    EXTRACT_B0
} from '../processes/preprocess.nf'

workflow PREPROCESSING {
    take:
        dwi_channel
        rev_channel
        anat_channel
        brain_mask_channel
        wm_mask_channel

    main:

        BET_DWI(dwi_channel)
        DENOISING(BET_DWI.out)
        
        topup_channel = dwi_channel
            .map{[it[0], it[2], it[3]]}
            .combine(DENOISING.out, by: 0)
            .combine(rev_channel, by: 0)
            .map{ sid, bvals, bvecs, dwi, rev -> tuple(sid, dwi, bvals, bvecs, rev)}

        TOPUP(topup_channel)

        eddy_channel = dwi_channel
            .map{[it[0], it[2], it[3]]}
            .combine(DENOISING.out, by: 0)
            .combine(TOPUP.out.topup_result, by: 0)
            .map{ sid, bvals, bvecs, dwi, corrected_b0s, field, movpar -> tuple(sid, dwi, bvals, bvecs, corrected_b0s, field, movpar)}

        EDDY_TOPUP(eddy_channel)
        
        n4_channel = EDDY_TOPUP.out.dwi_bval_bvec
            .combine(EDDY_TOPUP.out.b0_mask, by: 0)

        N4(n4_channel)

        dwi_crop_channel = N4.out
            .join(EDDY_TOPUP.out.b0_mask)
        
        CROP_DWI(dwi_crop_channel)

        anat_crop_channel = anat_channel

        if (params.run_bet_t2w) {
            BET_T2(anat_channel)
            anat_crop_channel = BET_T2.out.bet_t2
        }

        anat_crop_channel = anat_crop_channel
            .combine(brain_mask_channel, by: 0)
            .combine(wm_mask_channel, by:0)        

        CROP_ANAT(anat_crop_channel)
        RESAMPLE_ANAT(CROP_ANAT.out.cropped_t2w_and_mask)

        normalize_channel = CROP_DWI.out.dwi
            .combine(EDDY_TOPUP.out.dwi_bval_bvec.map{[it[0], it[2], it[3]]}, by: 0)
            .combine(CROP_DWI.out.mask, by: 0)

        NORMALIZE(normalize_channel)

        resample_dwi_channel = NORMALIZE.out.dwi_normalized
            .combine(CROP_DWI.out.mask, by: 0)
        
        RESAMPLE_DWI(resample_dwi_channel)

        extract_b0_channel = EDDY_TOPUP.out.dwi_bval_bvec
            .map{[it[0], it[2], it[3]]}
            .combine(RESAMPLE_DWI.out.dwi_resampled, by: 0)
            .map{ sid, bval, bvec, dwi -> tuple(sid, dwi, bval, bvec)}

        EXTRACT_B0(extract_b0_channel)

    emit:
        dwi_bval_bvec = extract_b0_channel
        b0_and_mask = EXTRACT_B0.out.b0_and_mask
        t2w_and_mask = RESAMPLE_ANAT.out.t2w_and_mask
}