#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process COMMIT {
    cpus params.processes_commit
    memory params.commit_memory_limit
    label "COMMIT"

    input:
        tuple val(sid), path(trk_h5), path(dwi), path(bval), path(bvec), path(peaks)
    output:
        tuple val(sid), path("${sid}__decompose_commit.h5"), emit: h5_commit
        tuple val(sid), path("${sid}__results_bzs/")
    when:
        params.run_commit
    
    script:
    ball_stick_arg=""
    perp_diff_arg=""
    if ( params.ball_stick ) {
        ball_stick_arg="--ball_stick"
    }
    else {
        perp_diff_arg="--perp_diff $params.perp_diff"
    }
    if ( params.use_commit2 ) {
    """
    scil_run_commit.py $trk_h5 $dwi $bval $bvec "${sid}__results_bzs/" --ball_stick --commit2 --in_peaks $peaks\
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir\
        --para_diff $params.para_diff $perp_diff_arg --iso_diff $params.iso_diff
    mv "${sid}__results_bzs/commit_2/decompose_commit.h5" "./${sid}__decompose_commit.h5"
    """
    }
    else {
    """
    scil_run_commit.py $trk_h5 $dwi $bval $bvec "${sid}__results_bzs/" --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir $ball_stick_arg \
        --para_diff $params.para_diff $perp_diff_arg --iso_diff $params.iso_diff
    mv "${sid}__results_bzs/commit_1/decompose_commit.h5" "./${sid}__decompose_commit.h5"
    """
    }
}

process COMMIT_ON_TRK {
    label "COMMIT"
    cpus params.processes_commit
    memory params.commit_memory_limit

    input:
        tuple val(sid), path(trk_h5), path(dwi), path(bval), path(bvec), path(peaks)
    output:
        tuple val(sid), path("${sid}__essential_tractogram.trk"), emit: trk_commit
        tuple val(sid), path("${sid}__results_bzs/")
    when:
        params.run_commit

    script:
    ball_stick_arg=""
    perp_diff_arg=""
    if ( params.ball_stick ) {
        ball_stick_arg="--ball_stick"
    }
    else {
        perp_diff_arg="--perp_diff $params.perp_diff"
    }
    """
    scil_run_commit.py $trk_h5 $dwi $bval $bvec "${sid}__results_bzs/" --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir $ball_stick_arg \
        --para_diff $params.para_diff $perp_diff_arg --iso_diff $params.iso_diff
    mv "${sid}__results_bzs/commit_1/essential_tractogram.trk" "./${sid}__essential_tractogram.trk"
    """
}