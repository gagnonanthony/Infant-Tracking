#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = false

def fetch_id ( dir, dir_base ) {
    return dir_base.relativize(dir)
        .collect{ it.name }
        .join("_")
}


workflow get_data {
    main:
        if ( !params.input )
            error "You must provide an input folder containing all images using --input"
        
        input = file(params.input)

        // Loading all images.
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { fetch_id(it.parent, input) }
        rev_channel = Channel.fromFilePairs("$input/**/*revb0.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*T2w.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        brain_mask_channel = Channel.fromFilePairs("$input/**/*brain_mask.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        wm_mask_channel = Channel.fromFilePairs("$input/**/*wm_mask.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_channel = dwi_channel.map{sid, bvals, bvecs, dwi -> tuple(sid, dwi, bvals, bvecs)}

    emit:
        dwi = dwi_channel
        rev = rev_channel
        anat = anat_channel
        brain_mask = brain_mask_channel
        wm_mask = wm_mask_channel
}