process PATHWAY_FBA {
    conda 'conda-forge::cobra=0.29'
    publishDir "${params.outdir}/03_fba", mode: 'copy'

    input:
    path model

    output:
    path 'fba_summary.tsv', emit: summary
    path 'envelope.tsv',    emit: envelope

    script:
    """
    pathway_fba.py --model ${model} --biomass ${params.biomass_id} --out .
    """
}
