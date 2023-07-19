#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process REGISTER_ANAT {
    label "REGISTER_ANAT"
    cpus params.processes_registration

    input:
        tuple val(sid), path(dwi), path(bval), path(t2w), path(brain_mask), path(wm_mask),
        path(fa) 
    output:
        tuple val(sid), path("${sid}__pwd_avg.nii.gz"), emit: pwd_avg
        tuple val(sid), 
        path("${sid}__t2w_warped.nii.gz"),
        path("${sid}__brain_mask_warped.nii.gz"),
        path("${sid}__wm_mask_warped.nii.gz"), emit: warped_anat
        tuple val(sid),
        path("output0GenericAffine.mat"),
        path("outputInverseWarped.nii.gz"),
        path("synoutput0Warp.nii.gz"),
        path("synoutput0InverseWarp.nii.gz"),
        path("maskoutput0Warp.nii.gz"),
        path("maskoutput0InverseWarp.nii.gz"), emit: transfos
    script:
    //** For some reason, mapping of the masks isn't as good as the t2w, performing a final SyN registration **//
    //** to fit the brain mask properly. **//
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    scil_compute_powder_average.py $dwi $bval ${sid}__pwd_avg.nii.gz\
        --b0_thr $params.b0_thr -f
    antsRegistration --dimensionality 3 --float 0\
        --output [output,outputWarped.nii.gz,outputInverseWarped.nii.gz]\
        --interpolation Linear --use-histogram-matching 0\
        --winsorize-image-intensities [0.005,0.995]\
        --initial-moving-transform [${sid}__pwd_avg.nii.gz,$t2w,1]\
        --transform Rigid['0.2']\
        --metric MI[${sid}__pwd_avg.nii.gz,$t2w,1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --transform Affine['0.2']\
        --metric MI[${sid}__pwd_avg.nii.gz,$t2w,1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --verbose 1
    mv outputWarped.nii.gz ${sid}__t2w_affine_warped.nii.gz
    antsRegistration --dimensionality 3 --float 0\
        --output [synoutput,synoutputWarped.nii.gz,synoutputInverseWarped.nii.gz]\
        --interpolation Linear --use-histogram-matching 0\
        --winsorize-image-intensities [0.005,0.995]\
        --transform SyN[0.1,3,0]\
        --metric MI[${sid}__pwd_avg.nii.gz,$t2w,1,32]\
        --metric CC[$fa,$t2w,1,4]\
        --convergence [200x150x200,1e-6,10] --shrink-factors 4x2x1\
        --smoothing-sigmas 3x2x1\
        --verbose 1
    mv synoutputWarped.nii.gz ${sid}__t2w_warped.nii.gz
    antsApplyTransforms -d 3 -i $brain_mask -r ${sid}__pwd_avg.nii.gz\
        -o ${sid}__brain_mask_initial_warp.nii.gz -n NearestNeighbor\
        -t synoutput0Warp.nii.gz output0GenericAffine.mat -u int
    antsRegistration -d 3 --float 0\
        --output [maskoutput,maskoutputWarped.nii.gz,maskoutputInverseWarped.nii.gz]\
        --interpolation NearestNeighbor --use-histogram-matching 0\
        --winsorize-image-intensities [0.005,0.995]\
        --transform SyN[0.1,3,0]\
        --metric MI[${sid}__pwd_avg.nii.gz,${sid}__brain_mask_initial_warp.nii.gz]\
        --convergence [500x250x125,1e-6,10] --shrink-factors 4x2x1\
        --smoothing-sigmas 3x2x1\
        --verbose 1
    mv maskoutputWarped.nii.gz ${sid}__brain_mask_warped.nii.gz
    antsApplyTransforms -d 3 -i $wm_mask -r ${sid}__pwd_avg.nii.gz\
        -o ${sid}__wm_mask_warped.nii.gz -n NearestNeighbor\
        -t maskoutput0Warp.nii.gz synoutput0Warp.nii.gz output0GenericAffine.mat -u int
    """
}