process SCORE_RANK {
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/04_ranking", mode: 'copy'

    input:
    path hits
    path struct

    output:
    path 'ranked_candidates.tsv', emit: ranked

    script:
    """
    score_rank.py ${hits} ${struct} ranked_candidates.tsv ${params.w_identity} ${params.w_plddt}
    """
}
