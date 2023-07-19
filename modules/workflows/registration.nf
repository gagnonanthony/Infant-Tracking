#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    REGISTER_ANAT
} from '../processes/registration_processes.nf'

workflow REGISTRATION {
    take:
        dwi_channel
        t2w_and_mask
        fa_channel
    main:

        register_channel = dwi_channel
            .map{[it[0], it[1], it[2]]}
            .combine(t2w_and_mask, by: 0)
            .combine(fa_channel, by: 0)

        REGISTER_ANAT(register_channel)

    emit:
        warped_anat = REGISTER_ANAT.out.warped_anat
}