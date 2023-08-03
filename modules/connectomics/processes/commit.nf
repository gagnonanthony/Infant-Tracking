#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process COMMIT2 {
    cpus params.processes_commit
    memory params.commit_memory_limit
    label "COMMIT2"

    input:
        tuple val(sid), path(h5), path(dwi), path(bval), path(bvec), path(peaks)
    output:
        tuple val(sid), path("${sid}__decompose_commit.h5"), emit: h5_commit
        tuple val(sid), path("${sid}__results_bzs/")
    script:
    """
    scil_run_commit.py $h5 $dwi $bval $bvec ${sid}__results_bzs/ --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b0_thr --nbr_dir $params.nbr_dir \
        --commit2 --ball_stick --para_diff $params.para_diff --iso_diff $params.iso_diff -v
    mv ${sid}__results_bzs/commit_2/decompose_commit.h5 ./${sid}__decompose_commit.h5
    """
}