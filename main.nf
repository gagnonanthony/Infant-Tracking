#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false

// Importing modules and processes
include { fetch_id;
          get_data_tracking;
          get_data_connectomics } from "./modules/io.nf"
include { PREPROCESSING } from "./modules/tracking/workflows/preprocessing.nf"
include { DTI } from "./modules/tracking/workflows/DTI.nf"
include { SH } from "./modules/tracking/workflows/SH.nf"
include { REGISTRATION } from "./modules/tracking/workflows/registration.nf"
include { FODF } from "./modules/tracking/workflows/FODF.nf"
include { TRACKING } from "./modules/tracking/workflows/tracking.nf"
include { CONNECTOMICS } from "./modules/connectomics/workflows/connectomics.nf"

workflow {
    if (params.help) display_usage()
    else {
        display_run_info()

        if ( params.run_tracking ) {
            data = get_data_tracking()

            PREPROCESSING(data.dwi, 
                          data.rev, 
                          data.anat, 
                          data.wm_mask)
            
            DTI(PREPROCESSING.out.dwi_bval_bvec,
                PREPROCESSING.out.b0_and_mask)

            if(params.sh_fitting) {
                SH(PREPROCESSING.out.dwi_bval_bvec)
            }

            md_channel = DTI.out.fa_and_md
                .map{ [it[0], it[2]]}

            REGISTRATION(md_channel,
                        PREPROCESSING.out.t2w_and_mask)
            
            b0_mask_channel = PREPROCESSING.out.b0_and_mask
                .map{[it[0], it[2]]}
            
            FODF(PREPROCESSING.out.dwi_bval_bvec,
                b0_mask_channel,
                DTI.out.fa_and_md)

            fa_channel = DTI.out.fa_and_md
                .map{[it[0], it[1]]}

            TRACKING(REGISTRATION.out.warped_anat,
                    FODF.out.fodf,
                    fa_channel)
        }

        if ( params.run_connectomics && params.run_tracking ) {
            tracking = TRACKING.out.trk

            // ** Labels needs to be provided as an input, since they are not computed at some point in the pipeline ** //
            input = file(params.input)
            labels = Channel.fromFilePairs("$input/**/*labels.nii.gz", size: 1, flat: true)
                        { fetch_id(it.parent, input) }

            dwi_peaks = PREPROCESSING.out.dwi_bval_bvec
                .combine(FODF.out.peaks, by: 0)
            fodf = FODF.out.fodf

            // ** Default metrics will be used with combined metrics provided in the input folder ** //
            provided_metrics = Channel.fromFilePairs("$input/**/metrics/*.nii.gz", size: -1, flat: true)
                                    { fetch_id(it.parent, input) }
            def_metrics = DTI.out.fa_and_md
                .combine(DTI.out.ad_and_rd, by: 0)
                .combine(FODF.out.afd_and_nufo, by: 0)
            metrics = provided_metrics
                .combine(def_metrics, by: 0)

            t2w = REGISTRATION.out.warped_anat

            transfos = REGISTRATION.out.transfos

            CONNECTOMICS(tracking,
                         labels,
                         dwi_peaks,
                         fodf,
                         metrics,
                         t2w,
                         transfos)
        }

        if ( params.run_connectomics && !params.run_tracking ) {
            data = get_data_connectomics()

            CONNECTOMICS(data.trk,
                         data.labels,
                         data.dwi_peaks,
                         data.fodf,
                         data.metrics,
                         data.t2w,
                         data.transfos)
        }
    }
}

if (!params.help) {
    workflow.onComplete = {
        log.info "Pipeline completed at : $workflow.complete"
        log.info "Execution status : ${ workflow.success ? 'COMPLETED' : 'FAILED'}"
        log.info "Execution duration : $workflow.duration"
    }
}

def display_usage () {
    usage = file("$projectDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["b0_thr":"$params.b0_thr",
                "initial_bet_f":"$params.initial_bet_f",
                "final_bet_f":"$params.final_bet_f",
                "run_bet_t2w":"$params.run_bet_t2w",
                "bet_t2w_f":"$params.bet_t2w_f",
                "topup_config":"$params.topup_config",
                "encoding_direction":"$params.encoding_direction",
                "readout":"$params.readout",
                "topup_prefix":"$params.topup_prefix",
                "eddy_cmd":"$params.eddy_cmd",
                "topup_bet_f":"$params.topup_bet_f",
                "use_slice_drop_correction":"$params.use_slice_drop_correction",
                "dwi_shell_tolerance":"$params.dwi_shell_tolerance",
                "fa_mask_threshold":"$params.fa_mask_threshold",
                "t2w_resolution":"$params.t2w_resolution",
                "t2w_interpolation":"$params.t2w_interpolation",
                "mask_interpolation":"$params.mask_interpolation",
                "dwi_resolution":"$params.dwi_resolution",
                "dwi_interpolation":"$params.dwi_interpolation",
                "mask_dwi_interpolation":"$params.mask_dwi_interpolation",
                "max_dti_shell_value":"$params.max_dti_shell_value",
                "sh_fitting":"$params.sh_fitting",
                "sh_fitting_order":"$params.sh_fitting_order",
                "sh_fitting_basis":"$params.sh_fitting_basis",
                "min_fodf_shell_value":"$params.min_fodf_shell_value",
                "max_fa_in_ventricle":"$params.max_fa_in_ventricle",
                "min_md_in_ventricle":"$params.min_md_in_ventricle",
                "relative_threshold":"$params.relative_threshold",
                "basis":"$params.basis",
                "sh_order":"$params.sh_order",
                "mean_frf":"$params.mean_frf",
                "fa":"$params.fa",
                "min_fa":"$params.min_fa",
                "min_nvox":"$params.min_nvox",
                "roi_radius":"$params.roi_radius",
                "set_frf":"$params.set_frf",
                "manual_frf":"$params.manual_frf",
                "fa_seeding_mask_thr":"$params.fa_seeding_mask_thr",
                "algo":"$params.algo",
                "seeding":"$params.seeding",
                "nb_seeds":"$params.nb_seeds",
                "tracking_seed":"$params.tracking_seed",
                "step_size":"$params.step_size",
                "theta":"$params.theta",
                "sfthres":"$params.sfthres",
                "min_len":"$params.min_len",
                "max_len":"$params.max_len",
                "erosion":"$params.erosion",
                "compress_value":"$params.compress_value",
                "output_dir":"$params.output_dir",
                "processes_denoise_dwi":"$params.processes_denoise_dwi",
                "processes_eddy":"$params.processes_eddy",
                "processes_registration":"$params.processes_registration",
                "processes_fodf":"$params.processes_fodf",
                "no_pruning":"$params.no_pruning",
                "no_remove_loops":"$params.no_remove_loops",
                "no_remove_outliers":"$params.no_remove_outliers",
                "min_length":"$params.min_length",
                "max_length":"$params.max_length",
                "loop_max_angle":"$params.loop_max_angle",
                "outlier_threshold":"$params.outlier_threshold",
                "nbr_dir":"$params.nbr_dir",
                "para_diff":"$params.para_diff",
                "iso_diff":"$params.iso_diff",
                "processes_commit":"$params.processes_commit",
                "processes_afd_fixel":"$params.processes_afd_fixel",
                "processes_connectivity":"$params.processes_connectivity",
                "run_tracking":"$params.run_tracking",
                "run_connectomics":"$params.run_connectomics"
                ]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
}

def display_run_info () {
    log.info ""
    log.info "Infant-DWI pipeline"
    log.info "========================"
    log.info "Pipeline adapted from the SCIL Tractoflow pipeline " 
    log.info "(https://github.com/scilus/tractoflow.git) and the "
    log.info "Connectoflow Pipeline (https://github.com/scilus/connectoflow.git)."
    log.info "Made for use on newborn diffusion MRI data."
    log.info ""
    log.info "Start time: $workflow.start"
    log.info ""

    log.debug "[Command-line]"
    log.debug "$workflow.commandLine"
    log.debug ""

    log.info "[Git Info]"
    log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
    log.info ""

    log.info "[Inputs]"
    log.info "Input: $params.input"
    log.info "Output Directory: $params.output_dir"
    log.info ""

    if ( params.run_tracking ) {
        log.info "[Tracking Options]"
        log.info ""
        log.info "GLOBAL OPTIONS"
        log.info "Threshold for b0: $params.b0_thr"
        log.info "DWI Shell Tolerance: $params.dwi_shell_tolerance"
        log.info ""
        log.info "BET DWI OPTIONS"
        log.info "Initial fractional value for BET: $params.initial_bet_f"
        log.info "Finale fractional value for BET: $params.final_bet_f"
        log.info ""
        log.info "BET T2W OPTIONS"
        log.info "Run BET on T2W image: $params.run_bet_t2w"
        log.info "Fractional value for T2W BET: $params.bet_t2w_f"
        log.info ""
        log.info "EDDY AND TOPUP OPTIONS"
        log.info "Configuration for topup: $params.topup_config"
        log.info "Encoding direction: $params.encoding_direction"
        log.info "Readout: $params.readout"
        log.info "Topup prefix: $params.topup_prefix"
        log.info "Topup BET fractional value: $params.topup_bet_f"
        log.info "Eddy command: $params.eddy_cmd"
        log.info "Run slice drop correction: $params.use_slice_drop_correction"
        log.info ""
        log.info "NORMALIZE OPTIONS"
        log.info "FA threshold for masking: $params.fa_mask_threshold"
        log.info ""
        log.info "RESAMPLE ANAT OPTIONS"
        log.info "Resampling resolution for T2W: $params.t2w_resolution"
        log.info "Interpolation method for T2W: $params.t2w_interpolation"
        log.info "Interpolation method for masks: $params.mask_interpolation"
        log.info ""
        log.info "RESAMPLE DWI OPTIONS"
        log.info "Resampling resolution for DWI: $params.dwi_resolution"
        log.info "Interpolation method for DWI: $params.dwi_interpolation"
        log.info "Interpolation method for DWI mask: $params.mask_dwi_interpolation"
        log.info ""
        log.info "EXTRACT DWI SHELLS OPTIONS"
        log.info "Maximum DTI shell value: $params.max_dti_shell_value"
        log.info ""
        log.info "SH FITTING OPTIONS"
        log.info "Run SH fitting: $params.sh_fitting"
        log.info "SH fitting order: $params.sh_fitting_order"
        log.info "SH fitting basis: $params.sh_fitting_basis"
        log.info ""
        log.info "FODF OPTIONS"
        log.info "Minimum fODF shell value: $params.min_fodf_shell_value"
        log.info "Maximum FA value in ventricles: $params.max_fa_in_ventricle"
        log.info "Minimum MD value in ventricles: $params.min_md_in_ventricle"
        log.info "Relative threshold (RT): $params.relative_threshold"
        log.info "SH basis: $params.basis"
        log.info "SH order: $params.sh_order"
        log.info ""
        log.info "FRF OPTIONS"
        log.info "Run mean FRF: $params.mean_frf"
        log.info "FA threshold for single fiber voxel: $params.fa"
        log.info "Minimum FA for selecting voxel: $params.min_fa"
        log.info "Minimum number of voxels: $params.min_nvox"
        log.info "ROI radius: $params.roi_radius"
        log.info "Set FRF: $params.set_frf"
        log.info "Manual FRF: $params.manual_frf"
        log.info ""
        log.info "SEEDING AND TRACKING OPTIONS"
        log.info "FA threshold for seeding mask: $params.fa_seeding_mask_thr"
        log.info "Erosion value to apply on brain mask: $params.erosion"
        log.info "Algorithm for tracking: $params.algo"
        log.info "Number of seeds per voxel: $params.nb_seeds"
        log.info "Seeding method: $params.seeding"
        log.info "Step size: $params.step_size"
        log.info "Theta threshold: $params.theta"
        log.info "Spherical function relative threshold: $params.sfthres"
        log.info "Minimum fiber length: $params.min_len"
        log.info "Maximum fiber length: $params.max_len"
        log.info "Random tracking seed: $params.tracking_seed"
        log.info "Compression value: $params.compress_value"
        log.info ""
        log.info "PROCESSES PER TASKS"
        log.info "Processes for denoising DWI: $params.processes_denoise_dwi"
        log.info "Processes for EDDY: $params.processes_eddy"
        log.info "Processes for registration: $params.processes_registration"
        log.info "Processes for FODF: $params.processes_fodf"
        log.info ""
    }

    if ( params.run_connectomics ) {
        log.info "[Connectomics Options]"
        log.info ""
        log.info "DECOMPOSE OPTIONS"
        log.info "No pruning: $params.no_pruning"
        log.info "No remove loops: $params.no_remove_loops"
        log.info "No remove outliers: $params.no_remove_outliers"
        log.info "Minimal outlier length: $params.min_length"
        log.info "Maximal outlier lenght: $params.max_length"
        log.info "Maximum looping angle: $params.loop_max_angle"
        log.info ""
        log.info "COMMIT OPTIONS"
        log.info "Number of directions: $params.nbr_dir"
        log.info "Parallel diffusivity: $params.para_diff"
        log.info "Isotropic diffusivity: $params.iso_diff"
        log.info ""
        log.info "PROCESSES OPTIONS"
        log.info "Number of processes for COMMIT: $params.processes_commit"
        log.info "Number of processes for AFD_FIXEL: $params.processes_afd_fixel"
        log.info "Number of processes for CONNECTIVITY: $params.processes_connectivity"
        log.info "" 
    }
}