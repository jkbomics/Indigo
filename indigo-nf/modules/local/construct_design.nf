process CONSTRUCT_DESIGN {
    conda "${projectDir}/envs/construct.yml"
    publishDir "${params.outdir}/06_constructs", mode: 'copy'
    errorStrategy 'terminate'   // input problems (bad member id, wrong mutation residue) should stop the run, not retry

    input:
    path fasta
    path ranked
    path constructs

    output:
    path 'shortlist_proteins.fasta',   emit: proteins
    path 'shortlist_codon_opt.fasta',  emit: cds
    path 'order_inserts.fasta',        emit: inserts
    path 'construct_summary.tsv',      emit: summary
    path 'qc_report.tsv',              emit: qc
    path 'construct_design.md',        emit: report
    path 'constructs',                 emit: genbank

    script:
    """
    design_constructs.py \\
        --fasta ${fasta} --ranked ${ranked} --constructs ${constructs} --outdir . \\
        --rbs ${params.rbs} --oh5 ${params.gg_oh5} --oh3 ${params.gg_oh3} \\
        --gc-min ${params.gc_min} --gc-max ${params.gc_max} --seed ${params.seed}
    """
}
