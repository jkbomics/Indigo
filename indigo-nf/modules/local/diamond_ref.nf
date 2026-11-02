process DIAMOND_REF {
    conda 'bioconda::diamond=2.1.9'
    publishDir "${params.outdir}/01_similarity", mode: 'copy'

    input:
    path cand
    path ref

    output:
    path 'cand_vs_ref.tsv', emit: hits

    script:
    """
    diamond makedb --in ${ref} -d ref --threads ${task.cpus}
    diamond blastp -d ref -q ${cand} -o cand_vs_ref.tsv \\
        --evalue 1e-5 --max-target-seqs 1 --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp
    touch cand_vs_ref.tsv
    """
}
