#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GENERATE_MASKS {
    label "GENERATE_MASKS"
    cpus 1

    input:
        tuple val(sid), path(t2w), path(wm_mask), path(fa)
    output:
        tuple val(sid), path("${sid}__wm_mask_final.nii.gz"),
        path("${sid}__brain_mask.nii.gz"), emit: masks
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    mrthreshold $fa ${sid}__fa_mask.nii.gz -abs $params.fa_seeding_mask_thr -nthreads 1 -force
    scil_image_math.py lower_threshold $t2w 0 brain_mask.nii.gz --data_type uint8
    scil_image_math.py erosion brain_mask.nii.gz $params.erosion ${sid}__brain_mask.nii.gz --data_type uint8 -f
    scil_image_math.py union ${sid}__fa_mask.nii.gz $wm_mask\
        wm_mask_temp.nii.gz --data_type uint8 -f
    scil_image_math.py intersection wm_mask_temp.nii.gz brain_mask.nii.gz\
        ${sid}__wm_mask_final.nii.gz --data_type uint8 -f
    """
}

process LOCAL_TRACKING {
    label "TRACKING"
    cpus 2

    input:
        tuple val(sid), path(fodf), path(wm_mask), path(brain_mask)
    output:
        tuple val(sid), path("${sid}__local_tracking.trk"), emit: tractogram
    script:
    """
    scil_compute_local_tracking.py $fodf $wm_mask $wm_mask\
        tmp.trk --algo $params.algo --$params.seeding $params.nb_seeds\
        --seed $params.tracking_seed --step $params.step_size --theta $params.theta\
        --sfthres $params.sfthres --min_length $params.min_len\
        --max_length $params.max_len --sh_basis $params.sh_fitting_basis\
        --compress $params.compress_value
    scil_remove_invalid_streamlines.py tmp.trk\
        ${sid}__local_tracking.trk --remove_single_point
    """
}