#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process COMMIT {
    cpus params.processes_commit
    memory params.commit_memory_limit

    input:
        tuple val(sid), path(h5), path(dwi), path(bval), path(bvec), path(peaks), path(para_diff), path(iso_diff)
    output:
        tuple val(sid), path("${sid}__decompose_commit.h5"), emit: h5_commit, optional: true
        tuple val(sid), path("${sid}__essential_tractogram.trk"), emit: trk_commit, optional: true
        tuple val(sid), path("${sid}__results_bzs/"), optional: true
        tuple val(sid), path("${sid}__results_bzs_1/"), optional: true
        tuple val(sid), path("${sid}__results_bzs_2/"), optional: true
    when:
        params.run_commit
    
    script:
    def para_diff_arg = para_diff ? "--para_diff \$(cat $para_diff)" : "--para_diff $params.para_diff"
    def iso_diff_arg = iso_diff ? "--iso_diff \$(cat $iso_diff)" : "--iso_diff $params.iso_diff" 
    def perp_diff_arg = ball_stick_arg ? "" : "--perp_diff $params.perp_diff"
    def ball_stick_arg = ball_stick_arg ? "--ball_stick" : ""

    if ( params.use_commit2 ) {
    """
    scil_run_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs/" --ball_stick --commit2 \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir\
        $para_diff_arg $iso_diff_arg
    mv "${sid}__results_bzs/commit_2/decompose_commit.h5" "./${sid}__decompose_commit.h5"
    """
    }
    else if ( params.use_both ) {
    """
    scil_run_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs_1/" --ball_stick --commit2 \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir\
        $para_diff_arg $iso_diff_arg
    scil_run_commit.py ${sid}__results_bzs_1/commit_2/essential_tractogram.trk $dwi $bval $bvec "${sid}__results_bzs_2/"\
        --in_peaks $peaks --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir\
        $para_diff_arg $iso_diff_arg $perp_diff_arg
    mv "${sid}__results_bzs_2/commit_1/essential_tractogram.trk" "./${sid}__essential_tractogram.trk"
    """
    }
    else {
    """
    scil_run_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs/" --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir $ball_stick_arg \
        $para_diff_arg $iso_diff_arg $perp_diff_arg
    mv "${sid}__results_bzs/commit_1/decompose_commit.h5" "./${sid}__decompose_commit.h5"
    """
    }
}

process COMMIT_ON_TRK {
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

process COMPUTE_PRIORS {
    cpus 1

    input:
        tuple val(sid), path(fa), path(ad), path(md)
    output:
        tuple val(sid), path("${sid}__para_diff.txt"), emit: para_diff
        tuple val(sid), path("${sid}__iso_diff.txt"), emit: iso_diff
    when:
        params.run_commit && params.compute_priors

    script:
    """
    scil_compute_NODDI_priors.py $fa $ad $md \
        --out_txt_1fiber ${sid}__para_diff.txt --out_txt_ventricles ${sid}__iso_diff.txt \
        --fa_min $params.fa_min_priors --fa_max $params.fa_max_priors --md_min $params.md_min_priors
    """
}