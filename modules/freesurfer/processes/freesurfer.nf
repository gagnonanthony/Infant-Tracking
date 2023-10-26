#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FREESURFER {
    cpus params.nb_threads

    input:
        tuple val(sid), path(anat)
    output:
        tuple val(sid), path("$sid/"), emit: folders

    script:
    """
    export SUBJECTS_DIR=.
    recon-all -i $anat -s $sid -all -parallel -openmp $params.nb_threads
    """
}