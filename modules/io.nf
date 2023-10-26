#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = false
params.references = false

def fetch_id ( dir, dir_base ) {
    return dir_base.relativize(dir)
        .collect{ it.name }
        .join("_")
}

// ** Getting data for the -profile freesurfer ** //
workflow get_data_freesurfer {
    main:
        if (! params.input ) {
            log.info "You must provide an input folder containing all images required for FreesurferFlow :"
            log.info "        --input=/path/to/[input_folder]           Input folder containing your subjects."
            log.info "                              [input]"
            log.info "                               ├-- S1"
            log.info "                               |   └-- *t1.nii.gz"
            log.info "                               └-- S2"
            log.info "                                   └-- *t1.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)

        // ** Loading files ** //
        anat_channel = Channel.fromFilePairs("$input/**/*t1.nii.gz", size: 1, flat: true)

    emit:
        anat = anat_channel
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
            log.info "                               |   ├-- *dwi.bval"
            log.info "                               |   ├-- *dwi.bvec"
            log.info "                               |   ├-- *revb0.nii.gz"
            log.info "                               |   └-- *t1.nii.gz"
            log.info "                               └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"
            log.info "                                    ├-- *bval"
            log.info "                                    ├-- *bvec"
            log.info "                                    ├-- *revb0.nii.gz"
            log.info "                                    └-- *t1.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }
        
        input = file(params.input)

        // ** Loading all files. ** //
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { fetch_id(it.parent, input) }
        rev_channel = Channel.fromFilePairs("$input/**/*revb0.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*t1.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }

        // ** Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility. ** //
        dwi_channel = dwi_channel.map{sid, bvals, bvecs, dwi -> tuple(sid, dwi, bvals, bvecs)}
        
    emit:
        dwi = dwi_channel
        rev = rev_channel
        anat = anat_channel
}

// ** Getting data for the -profile tracking,infant ** //
workflow get_data_tracking_infant {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "        --input=/path/to/[input_folder]             Input folder containing multiple subjects for tracking"
            log.info ""
            log.info "                               [Input]"
            log.info "                               ├-- S1"
            log.info "                               |   ├-- *dwi.nii.gz"
            log.info "                               |   ├-- *dwi.bval"
            log.info "                               |   ├-- *dwi.bvec"
            log.info "                               |   ├-- *revb0.nii.gz"
            log.info "                               |   ├-- *t2w.nii.gz"
            log.info "                               |   └-- *wm_mask.nii.gz"
            log.info "                               └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"
            log.info "                                    ├-- *bval"
            log.info "                                    ├-- *bvec"
            log.info "                                    ├-- *revb0.nii.gz"
            log.info "                                    ├-- *t2w.nii.gz"
            log.info "                                    └-- *wm_mask.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }
        
        input = file(params.input)

        // ** Loading all files. ** //
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { fetch_id(it.parent, input) }
        rev_channel = Channel.fromFilePairs("$input/**/*revb0.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*t2w.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        wm_mask_channel = Channel.fromFilePairs("$input/**/*wm_mask.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }

        // ** Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility. ** //
        dwi_channel = dwi_channel.map{sid, bvals, bvecs, dwi -> tuple(sid, dwi, bvals, bvecs)}
        
    emit:
        dwi = dwi_channel
        rev = rev_channel
        anat = anat_channel
        wm_mask = wm_mask_channel
}

// ** Fetching data for -profile connectomics ** //
workflow get_data_connectomics {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "     --input=/path/to/[input_folder]             Input folder containing multiple subjects"
            log.info ""
            log.info "                                [Input]"
            log.info "                                ├-- S1"
            log.info "                                |   ├-- *dwi.nii.gz"            
            log.info "                                |   ├-- *dwi.bval"            
            log.info "                                |   ├-- *dwi.bvec"                
            log.info "                                |   ├-- *t1.nii.gz"             
            log.info "                                |   ├-- *.trk"                  
            log.info "                                |   ├-- *labels.nii.gz"          
            log.info "                                |   ├-- *peaks.nii.gz"          
            log.info "                                |   ├-- *fodf.nii.gz"            
            log.info "                                |   ├-- OGenericAffine.mat"     
            log.info "                                |   ├-- output1Warp.nii.gz"  
            log.info "                                |   └-- metrics"
            log.info "                                |       └-- METRIC_NAME.nii.gz  [Optional]"
            log.info "                                └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"         
            log.info "                                    ├-- *bval"                  
            log.info "                                    ├-- *bvec"                  
            log.info "                                    ├-- *t1.nii.gz"           
            log.info "                                    ├-- *.trk"                   
            log.info "                                    ├-- *labels.nii.gz"         
            log.info "                                    ├-- *peaks.nii.gz"       
            log.info "                                    ├-- *fodf.nii.gz"         
            log.info "                                    ├-- OGenericAffine.mat"   
            log.info "                                    ├-- output1Warp.nii.gz"  
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
        t2w_channel = Channel.fromFilePairs("$input/**/*t1.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        transfos_channel = Channel.fromFilePairs("$input/**/{0GenericAffine.mat,output1Warp.nii.gz}", size: 2, flat: true)
            { fetch_id(it.parent, input) }

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_peaks_channel = dwi_peaks_channel.map{sid, bvals, bvecs, dwi, peaks -> tuple(sid, dwi, bvals, bvecs, peaks)}

        emit:
            trk = tracking_channel
            labels = labels_channel
            dwi_peaks = dwi_peaks_channel
            fodf = fodf_channel
            metrics = metrics_channel
            t2w = t2w_channel
            transfos = transfos_channel
}

// ** Fetching data for -profile connectomics,infant ** //
workflow get_data_connectomics_infant {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "     --input=/path/to/[input_folder]             Input folder containing multiple subjects"
            log.info ""
            log.info "                                [Input]"
            log.info "                                ├-- S1"
            log.info "                                |   ├-- *dwi.nii.gz"            
            log.info "                                |   ├-- *dwi.bval"            
            log.info "                                |   ├-- *dwi.bvec"                
            log.info "                                |   ├-- *t2w.nii.gz"             
            log.info "                                |   ├-- *.trk"                  
            log.info "                                |   ├-- *labels.nii.gz"          
            log.info "                                |   ├-- *peaks.nii.gz"          
            log.info "                                |   ├-- *fodf.nii.gz"            
            log.info "                                |   ├-- OGenericAffine.mat"     
            log.info "                                |   ├-- output1Warp.nii.gz"  
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
            log.info "                                    ├-- output1Warp.nii.gz"  
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
        t2w_channel = Channel.fromFilePairs("$input/**/*t2w.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        transfos_channel = Channel.fromFilePairs("$input/**/{0GenericAffine.mat,output1Warp.nii.gz}", size: 2, flat: true)
            { fetch_id(it.parent, input) }

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_peaks_channel = dwi_peaks_channel.map{sid, bvals, bvecs, dwi, peaks -> tuple(sid, dwi, bvals, bvecs, peaks)}

        emit:
            trk = tracking_channel
            labels = labels_channel
            dwi_peaks = dwi_peaks_channel
            fodf = fodf_channel
            metrics = metrics_channel
            t2w = t2w_channel
            transfos = transfos_channel
}

workflow get_data_template {
    main:
    if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "        --input=/path/to/[input_folder]             Input folder containing multiple subjects for tracking"
            log.info ""
            log.info "                               [Input]"
            log.info "                               ├-- S1"
            log.info "                               |   ├-- *dwi.nii.gz"
            log.info "                               |   ├-- *dwi.bvec"
            log.info "                               |   ├-- *fa.nii.gz"
            log.info "                               |   ├-- *t2w.nii.gz"
            log.info "                               └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"
            log.info "                                    ├-- *dwi.bvec"
            log.info "                                    ├-- *fa.nii.gz"
            log.info "                                    └-- *t2w.nii.gz"
            log.info "                               [References]"
            log.info "                               ├-- *fa_ref.nii.gz"
            log.info "                               └-- *t2w_ref.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }
        
        input = file(params.input)
        references = file(params.references)

        // Loading all files.
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bvec}", size: 2, flat: true)
            { fetch_id(it.parent, input) }
        fa_channel = Channel.fromFilePairs("$input/**/*fa.nii.gz", size:1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*t2w.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_ref = Channel.fromPath("$references/*t2w_ref.nii.gz")
        fa_ref = Channel.fromPath("$references/*fa_ref.nii.gz")

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_channel = dwi_channel.map{sid, bvecs, dwi -> tuple(sid, dwi, bvecs)}
        
    emit:
        dwi = dwi_channel
        anat = anat_channel
        fa = fa_channel
        anat_ref = anat_ref
        fa_ref = fa_ref
}