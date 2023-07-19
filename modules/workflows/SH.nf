#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    SH_FITTING_SHELL;
    SH_FITTING
} from '../processes/SH.nf'

workflow SH {
    take:
        dwi_bval_bvec
    main:
        SH_FITTING_SHELL(dwi_bval_bvec)
        SH_FITTING(SH_FITTING_SHELL.out.dwi_sh)
}