#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRANSFORM_LABELS } from "../processes/transform_labels.nf"
include { DECOMPOSE_CONNECTIVITY } from "../processes/decompose.nf"
include { COMMIT;
          COMMIT_ON_TRK } from "../processes/commit.nf"
include { COMPUTE_AFD_FIXEL;
          COMPUTE_CONNECTIVITY } from "../processes/compute_metrics.nf"
include { VISUALIZE_CONNECTIVITY } from "../processes/viz.nf"

workflow CONNECTOMICS {
    take:
        tracking_channel
        labels_channel
        dwi_peaks_channel
        fodf_channel
        metrics_channel
        t2w_channel
        transfos_channel
        
    main:

        // ** Transforming labels to diff space ** //
        channel_for_transfo = labels_channel
            .combine(t2w_channel, by: 0)
            .combine(transfos_channel, by: 0)
        TRANSFORM_LABELS(channel_for_transfo)

        // ** If -profile infant is used, first part will be run. COMMIT1 is the only supported ** //
        // ** method as of now, since running commit2 requires a decomposition first, which is not an ** //
        // ** easy task on infant data. This will be improved in the future. ** //
        if ( params.infant_config ) {

            // ** COMMIT1 processing on trk ** //
            commit_channel = tracking_channel
                                .combine(dwi_peaks_channel, by: 0)
            COMMIT_ON_TRK(commit_channel)

            // ** Decomposing tractogram ** //
            decompose_channel = COMMIT_ON_TRK.out.trk_commit
                                .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            DECOMPOSE_CONNECTIVITY(decompose_channel)

            // ** Setting output channel ** //
            afd_fixel_channel = DECOMPOSE_CONNECTIVITY.out.decompose
                                    .combine(fodf_channel, by: 0)
        }
        else {
            // ** Decomposing tractogram ** //
            decompose_channel = tracking_channel
                                    .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            DECOMPOSE_CONNECTIVITY(decompose_channel)

            // ** Running COMMIT1 or COMMIT2 ** //
            commit_channel = DECOMPOSE_CONNECTIVITY.out.decompose
                                .combine(dwi_peaks_channel, by: 0)
            COMMIT(commit_channel)

            // ** Setting output channel ** //
            afd_fixel_channel = COMMIT.out.h5_commit
                                    .combine(fodf_channel, by: 0)
        }

        // ** Computing AFD fixel ** //
        COMPUTE_AFD_FIXEL(afd_fixel_channel)

        // ** Computing Connectivity ** //
        compute_metrics_channel = COMPUTE_AFD_FIXEL.out.decompose_afd
            .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            .combine(metrics_channel, by: 0)
        COMPUTE_CONNECTIVITY(compute_metrics_channel)
        
        // ** Visualizing Connectivity ** //
        VISUALIZE_CONNECTIVITY(COMPUTE_CONNECTIVITY.out.metrics)
}