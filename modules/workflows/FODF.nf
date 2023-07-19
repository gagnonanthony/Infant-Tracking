#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    FODF_SHELL;
    COMPUTE_FRF;
    MEAN_FRF;
    FODF_METRICS
} from '../processes/FODF.nf'

workflow FODF {
    take:
        dwi_bval_bvec
        b0_mask
        fa_and_md

    main:

        FODF_SHELL(dwi_bval_bvec)

        frf_channel = dwi_bval_bvec
            .combine(b0_mask, by: 0)

        COMPUTE_FRF(frf_channel)

        all_frf = COMPUTE_FRF.out.frf
            .map{[it[1]]}
            .collect()

        frf_for_fodf = COMPUTE_FRF.out.frf

        if (params.mean_frf) {
            MEAN_FRF(all_frf)
            
            frf_for_fodf = COMPUTE_FRF.out.frf
                .merge(MEAN_FRF.out.mean_frf)
                .map{it -> [it[0], it[2]]}
        }

        fodf_channel = FODF_SHELL.out.dwi_fodf
            .combine(b0_mask, by: 0)
            .combine(fa_and_md, by: 0)
            .combine(frf_for_fodf, by: 0)
        
        FODF_METRICS(fodf_channel)

    emit:
        fodf = FODF_METRICS.out.fodf
}