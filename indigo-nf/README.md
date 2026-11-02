# indigo-nf

In silico screening workflow for microbial indigo (E. coli host).

    nextflow run main.nf -profile conda --model iML1515.xml --outdir results

Steps: UniProt fetch -> DIAMOND vs reference enzymes -> AlphaFold DB pLDDT -> FBA (COBRApy) -> ranking -> wet-lab handoff.

Before running:
1. Fill `assets/reference_enzymes.csv` and `assets/candidates_template.csv` (replace every TO_CURATE with a verified UniProt accession).
2. Download iML1515 SBML from BiGG.
3. Pin tool versions in the module `conda` directives.

Scoring weights in `nextflow.config` are placeholders.
