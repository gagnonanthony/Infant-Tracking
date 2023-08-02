#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = false

def fetch_id ( dir, dir_base ) {
    return dir_base.relativize(dir)
        .collect{ it.name }
        .join("_")
}

// ** Decided to split the data fetching steps for different profiles in different functions ** //
// ** for easier code-reading. ** //
workflow get_data_tracking {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "        --input=/path/to/[input_folder]             Input folder containing multiple subjects for tracking"
            log.info ""
            log.info "                               [Input]"
            log.info "                               ├-- S1"
            log.info "                               |   ├-- *dwi.nii.gz"
            log.info "                               |   ├-- *bval"
            log.info "                               |   ├-- *bvec"
            log.info "                               |   ├-- *revb0.nii.gz"
            log.info "                               |   ├-- *T2w.nii.gz"
            log.info "                               |   ├-- *brain_mask.nii.gz"
            log.info "                               |   └-- *wm_mask.nii.gz"
            log.info "                               └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"
            log.info "                                    ├-- *bval"
            log.info "                                    ├-- *bvec"
            log.info "                                    ├-- *revb0.nii.gz"
            log.info "                                    ├-- *T2w.nii.gz"
            log.info "                                    ├-- *brain_mask.nii.gz"
            log.info "                                    └-- *wm_mask.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }
        
        input = file(params.input)

        // Loading all files.
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { fetch_id(it.parent, input) }
        rev_channel = Channel.fromFilePairs("$input/**/*revb0.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*t2w.nii.gz", size: 1, flat: true)
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

workflow get_data_connectomics {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "     --input=/path/to/[input_folder]             Input folder containing multiple subjects"
            log.info ""
            log.info "                                [Input]"
            log.info "                                ├-- S1"
            log.info "                                |   ├-- *dwi.nii.gz"            
            log.info "                                |   ├-- *bval"            
            log.info "                                |   ├-- *bvec"                
            log.info "                                |   ├-- *t2w.nii.gz"             
            log.info "                                |   ├-- *.trk"                  
            log.info "                                |   ├-- *labels.nii.gz"          
            log.info "                                |   ├-- *peaks.nii.gz"          
            log.info "                                |   ├-- *fodf.nii.gz"            
            log.info "                                |   ├-- OGenericAffine.mat"     
            log.info "                                |   ├-- synoutput0Warp.nii.gz"  
            log.info "                                |   ├-- maskoutput0Warp.nii.gz" 
            log.info "                                |   └-- metrics"
            log.info "                                |       └-- METRIC_NAME.nii.gz  [Optional]"
            log.info "                                └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"         
            log.info "                                    ├-- *bval"                  
            log.info "                                    ├-- *bvec"                  
            log.info "                                    ├-- *t2w.nii.gz"           
            log.info "                                    ├-- *.trk"                   
            log.info "                                    ├-- *labels.nii.gz"         
            log.info "                                    ├-- *peaks.nii.gz"       
            log.info "                                    ├-- *fodf.nii.gz"         
            log.info "                                    ├-- OGenericAffine.mat"   
            log.info "                                    ├-- synoutput0Warp.nii.gz"  
            log.info "                                    ├-- maskoutput0Warp.nii.gz" 
            log.info "                                    └-- metrics"
            log.info "                                        └-- METRIC_NAME.nii.gz  [Optional]"
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)

        // Loading all files.
        tracking_channel = Channel.fromFilePairs("$input/**/*.trk", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        labels_channel = Channel.fromFilePairs("$input/**/*labels.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        dwi_peaks_channel = Channel.fromFilePairs("$input/**/{*dwi.nii.gz,*.bval,*.bvec,*peaks.nii.gz}", size: 4, flat: true)
            { fetch_id(it.parent, input) }
        fodf_channel = Channel.fromFilePairs("$input/**/*fodf.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        metrics_channel = Channel.fromFilePairs("$input/**/metrics/*.nii.gz", size: -1, maxDepth: 2)
            { it.parent.parent.name }
        t2w_channel = Channel.fromFilePairs("$input/**/*t2w_warped.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        transfos_channel = Channel.fromFilePairs("$input/**/{0GenericAffine.mat,synoutput0Warp.nii.gz,maskoutput0Warp.nii.gz}", size: 3, flat: true)
            { fetch_id(it.parent, input) }

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_peaks_channel = dwi_peaks_channel.map{sid, bvals, bvecs, dwi, peaks -> tuple(sid, dwi, bvals, bvecs, peaks)}

        // Setting up transfos channel in this order : sid, affine, syn, masksyn
        transfos_channel = transfos_channel.map{sid, affine, masksyn, syn -> tuple(sid, affine, syn, masksyn)}

        // Flattening metrics channel.
        metrics_channel = metrics_channel.transpose().groupTuple()
                                .flatMap{ sid, metrics -> metrics.collect{ [sid, it] } }

        emit:
            trk = tracking_channel
            labels = labels_channel
            dwi_peaks = dwi_peaks_channel
            fodf = fodf_channel
            metrics = metrics_channel
            t2w = t2w_channel
            transfos = transfos_channel
}