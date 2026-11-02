#!/usr/bin/env nextflow
include { FETCH_SEQS as FETCH_CAND } from './modules/local/fetch_seqs'
include { FETCH_SEQS as FETCH_REF  } from './modules/local/fetch_seqs'
include { DIAMOND_REF              } from './modules/local/diamond_ref'
include { STRUCT_QC                } from './modules/local/struct_qc'
include { PATHWAY_FBA              } from './modules/local/pathway_fba'
include { SCORE_RANK               } from './modules/local/score_rank'
include { HANDOFF                  } from './modules/local/handoff'

workflow {
    if (!params.model) {
        error "Provide --model (E. coli GEM in SBML, e.g. iML1515 from http://bigg.ucsd.edu/models/iML1515)"
    }

    cand_csv = channel.fromPath(params.candidates, checkIfExists: true)
    ref_csv  = channel.fromPath(params.references, checkIfExists: true)
    model    = channel.fromPath(params.model,      checkIfExists: true)

    // 1. sequences from UniProt
    FETCH_CAND(cand_csv)
    FETCH_REF(ref_csv)

    // 2. similarity of candidates to characterised indigo-forming enzymes
    DIAMOND_REF(FETCH_CAND.out.fasta, FETCH_REF.out.fasta)

    // 3. AlphaFold DB structure confidence
    STRUCT_QC(FETCH_CAND.out.fasta)

    // 4. pathway feasibility / theoretical yield in the host GEM
    PATHWAY_FBA(model)

    // 5. rank + wet-lab handoff package
    SCORE_RANK(DIAMOND_REF.out.hits, STRUCT_QC.out.table)
    HANDOFF(SCORE_RANK.out.ranked, PATHWAY_FBA.out.summary)
}
