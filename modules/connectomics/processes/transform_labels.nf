#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRANSFORM_LABELS {
    cpus 1
    memory '2 GB'
    label "TRANSFORM_LABELS"

    input:
        tuple val(sid), path(labels), path(t2), path(mat), path(syn), path(masksyn)
    output:
        tuple val(sid), path("${sid}__labels_warped.nii.gz"), emit: labels_warped
    script:
    """
    antsApplyTransforms -d 3 -i $labels -r $t2 -o ${sid}__labels_warped.nii.gz \
        -t $masksyn $syn $mat -n NearestNeighbor
    scil_image_math.py convert ${sid}__labels_warped.nii.gz ${sid}__labels_warped.nii.gz \
        --data_type int16 -f
    """
}