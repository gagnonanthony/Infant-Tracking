#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FREESURFER {
    cpus params.nb_threads
    if ( ! params.symlink ) {
        publishDir "${params.output_dir}/Freesurfer/Freesurfer/", mode: 'copy'
    } else {
        publishDir "${params.output_dir}/Freesurfer/Freesurfer/", mode: 'symlink'
    }

    input:
        tuple val(sid), path(anat)
    output:
        tuple val(sid), path("$sid/"), emit: folders
        tuple val(sid), path("${sid}__final_t1.nii.gz"), emit: final_t1

    script:
    """
    export SUBJECTS_DIR=.
    recon-all -i $anat -s $sid -all -parallel -openmp $params.nb_threads
    mri_convert $sid/mri/antsdn.brain.mgz ${sid}__final_t1.nii.gz
    """
}