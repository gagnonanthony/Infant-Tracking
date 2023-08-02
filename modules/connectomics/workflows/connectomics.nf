#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRANSFORM_LABELS } from "../processes/transform_labels.nf"
include { DECOMPOSE_CONNECTIVITY } from "../processes/decompose.nf"
include { COMMIT2 } from "../processes/commit.nf"
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

        channel_for_transfo = labels_channel
            .combine(t2w_channel, by: 0)
            .combine(transfos_channel, by: 0)

        TRANSFORM_LABELS(channel_for_transfo)

        decompose_channel = tracking_channel
            .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)

        DECOMPOSE_CONNECTIVITY(decompose_channel)

        commit_channel = DECOMPOSE_CONNECTIVITY.out.decompose
            .combine(dwi_peaks_channel, by: 0)

        COMMIT2(commit_channel)

        afd_fixel_channel = COMMIT2.out.h5_commit
            .combine(fodf_channel, by: 0)

        COMPUTE_AFD_FIXEL(afd_fixel_channel)

        compute_metrics_channel = COMPUTE_AFD_FIXEL.out.decompose_afd
            .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            .combine(metrics_channel, by: 0)
        
        COMPUTE_CONNECTIVITY(compute_metrics_channel)
        
        VISUALIZE_CONNECTIVITY(COMPUTE_CONNECTIVITY.out.metrics)

}