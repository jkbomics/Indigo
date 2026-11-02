process HANDOFF {
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/05_handoff", mode: 'copy'

    input:
    path ranked
    path fba

    output:
    path 'handoff.md',      emit: report
    path 'doe_matrix.csv',  emit: doe

    script:
    """
    handoff.py ${ranked} ${fba} ${params.top_n} ${params.doe_runs}
    """
}
