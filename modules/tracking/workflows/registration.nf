#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    REGISTER_ANAT
} from '../processes/registration_processes.nf'

workflow REGISTRATION {
    take:
        md_channel
        t2w_and_mask
    main:
        
        register_channel = md_channel
            .combine(t2w_and_mask, by: 0)

        REGISTER_ANAT(register_channel)

    emit:
        warped_anat = REGISTER_ANAT.out.warped_anat
        transfos = REGISTER_ANAT.out.transfos
}