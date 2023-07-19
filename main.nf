#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false

// Importing modules and processes
include { get_data } from "./modules/io.nf"
include { PREPROCESSING } from "./modules/workflows/preprocessing.nf"
include { DTI } from "./modules/workflows/DTI.nf"
include { SH } from "./modules/workflows/SH.nf"
include { REGISTRATION } from "./modules/workflows/registration.nf"
include { FODF } from "./modules/workflows/FODF.nf"
include { TRACKING } from "./modules/workflows/tracking.nf"

workflow {
    if (params.help) display_usage()
    else {
        display_run_info()
        data = get_data()

        PREPROCESSING(data.dwi, 
                      data.rev, 
                      data.anat, 
                      data.brain_mask, 
                      data.wm_mask)
        
        DTI(PREPROCESSING.out.dwi_bval_bvec,
            PREPROCESSING.out.b0_and_mask)

        if(params.sh_fitting) {
            SH(PREPROCESSING.out.dwi_bval_bvec)
        }

        fa_channel = DTI.out.fa_and_md
            .map{[it[0], it[1]]}

        REGISTRATION(PREPROCESSING.out.dwi_bval_bvec,
                     PREPROCESSING.out.t2w_and_mask,
                     fa_channel)
        
        b0_mask_channel = PREPROCESSING.out.b0_and_mask
            .map{[it[0], it[2]]}
        
        FODF(PREPROCESSING.out.dwi_bval_bvec,
             b0_mask_channel,
             DTI.out.fa_and_md)
        
        masks_channel = REGISTRATION.out.warped_anat
            .map{[it[0], it[2], it[3]]}

        TRACKING(masks_channel,
                 FODF.out.fodf,
                 fa_channel)

        workflow.onComplete {
            log.info "Pipeline completed at : $workflow.complete"
            log.info "Execution status : ${ workflow.success ? 'COMPLETED' : 'FAILED'}"
            log.info "Execution duration : $workflow.duration"}
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
                "use_brain_mask_as_tracking_mask":"$params.use_brain_mask_as_tracking_mask",
                "processes_denoise_dwi":"$params.processes_denoise_dwi",
                "processes_eddy":"$params.processes_eddy",
                "processes_registration":"$params.processes_registration",
                "processes_fodf":"$params.processes_fodf"
                ]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
}

def display_run_info () {
    log.info ""
    log.info "Infant-Tracking pipeline"
    log.info "========================"
    log.info "Pipeline adapted from the SCIL Tractoflow pipeline " 
    log.info "(https://github.com/scilus/tractoflow.git)"
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

    log.info "Input: $params.input"
}