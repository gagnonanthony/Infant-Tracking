#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    EXTRACT_DTI_SHELL;
    DTI_METRICS
} from '../processes/DTI_processes.nf'

workflow DTI {
    take:
        dwi_channel
        b0_and_mask_channel
    main:

        EXTRACT_DTI_SHELL(dwi_channel)

        dti_channel = b0_and_mask_channel
            .map{[it[0], it[2]]}
            .combine(EXTRACT_DTI_SHELL.out.dti_files, by: 0)
            .map{ sid, mask, dwi, bval, bvec -> tuple(sid, dwi, bval, bvec, mask)}
        
        DTI_METRICS(dti_channel)

    emit:
        fa_and_md = DTI_METRICS.out.fa_and_md
}