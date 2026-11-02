process STRUCT_QC {
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/02_structure", mode: 'copy'

    input:
    path cand_fasta

    output:
    path 'structure_qc.tsv', emit: table
    path 'models',           emit: models

    script:
    """
    struct_qc.py ${cand_fasta} structure_qc.tsv models
    """
}
