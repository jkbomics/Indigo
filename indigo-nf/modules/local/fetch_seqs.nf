process FETCH_SEQS {
    tag "${csv.baseName}"
    conda 'conda-forge::python=3.11'

    input:
    path csv

    output:
    path "${csv.baseName}.fasta", emit: fasta

    script:
    """
    fetch_uniprot.py ${csv} ${csv.baseName}.fasta
    """
}
