# indigo-nf

In silico screening workflow for microbial indigo (E. coli host).

    nextflow run main.nf -profile conda --model iML1515.xml --outdir results

Steps: UniProt fetch -> DIAMOND vs reference enzymes -> AlphaFold DB pLDDT -> FBA (COBRApy) -> ranking -> wet-lab handoff.

Before running:
1. Fill `assets/reference_enzymes.csv` and `assets/candidates_template.csv` (replace every TO_CURATE with a verified UniProt accession).
2. Download iML1515 SBML from BiGG.
3. Pin tool versions in the module `conda` directives.

Scoring weights in `nextflow.config` are placeholders.

## Finding new candidates

    python3 tools/find_candidates.py            # queries UniProt, writes candidates_review.tsv + candidates_draft.csv
    # check names/organisms in the review table, then copy chosen rows into assets/candidates_template.csv

`ctrl_*` rows in the candidates file are positive controls that are also in the reference set, so they score ~100% identity by construction. New homologs go below them.
StyA/StyB accessions are not filled in: resolve them with the script above and verify before use.

## Construct design (module 6)

Reads `assets/constructs.csv` (which shortlisted candidates go into which construct, optional point mutations),
codon-optimises each CDS for E. coli with DNA Chisel, and writes into `results/06_constructs/`:
`shortlist_proteins.fasta`, `shortlist_codon_opt.fasta`, `order_inserts.fasta` (PET and Golden Gate cassettes),
`constructs/*.gb`, `qc_report.tsv`, `construct_summary.tsv`, `construct_design.md`.
Mutations use UniProt 1-based numbering; the script stops if the wild-type residue does not match.
Defaults for RBS and Golden Gate overhangs are starting points: match them to your vector/toolkit before ordering.
