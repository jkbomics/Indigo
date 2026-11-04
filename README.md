# indigo-nf: in silico screening for microbial indigo

A Nextflow pipeline that takes a shortlist of candidate indigo-forming enzymes, checks them against known enzymes and predicted structures, tests the pathway in an *E. coli* genome-scale model, ranks the candidates, and writes codon-optimised construct designs for the wet lab.

Status: all seven analysis steps plus the construct-design module ran end to end in a GitHub Codespace on 2026-11-03 (Nextflow 26.10.4). Everything here is a computational prediction or design. Nothing has been tested in the lab yet.

---

## 1. What this project is for

The goal is a standardised, reproducible computational workflow for blue indigo dye made by microbes, ending in a package a wet lab can act on. The scope is:

1. Pick a production host that is easy to scale commercially.
2. Collect public data on the enzymes that make indigo.
3. Rank candidate enzymes and check the pathway is feasible in the host.
4. Design the DNA constructs to order and test.

What this workflow cannot do: it does not predict colourfastness, vat-dyeing behaviour, fibre binding, or whether a given enzyme really oxidises indole. Those need wet-lab data.

Indigo here means the indole dimer (CAS 482-89-3). It is not indigoidine, the blue pigment made by an NRPS in some Streptomyces and Erwinia strains, which has different chemistry.

## 2. Biology in one paragraph

Microbial indigo mostly follows one route: L-tryptophan is converted to indole by tryptophanase (TnaA), an oxygenase converts indole to indoxyl, and two indoxyl molecules dimerise spontaneously (in the presence of oxygen) to indigo. The oxygenase is the choice that matters. Options in this project are a bacterial flavin-containing monooxygenase (FMO), styrene monooxygenase (StyA with its reductase StyB), and naphthalene dioxygenase (NDO). Side products are indirubin (red/purple) and isatin, which reduce colour purity. A second published route makes indole straight from the aromatic pathway through an indole lyase, bypassing tryptophan. That route is not in the FBA model yet (see section 9).

## 3. Microbes and enzymes used

### Production host

***Escherichia coli* K-12 lineage (MG1655 / W3110).** Chosen because:

- it is the best characterised host, with a curated genome-scale model (iML1515)
- it already has native tryptophanase (tnaA), so only an oxygenase has to be added
- tryptophan-overproducing strains exist, so the route can run from glucose
- fast growth, standard tools, and established fermenter scale-up
- indigo is insoluble and accumulates as a precipitate, which helps recovery

Alternatives considered: *Corynebacterium glutamicum* (strong tryptophan producer, fewer genetic tools) and *Pseudomonas putida* (more indole-tolerant, harder to engineer).

### Enzyme sources (UniProt accessions)

Organism names are as listed in UniProt on 2026-11-03.

| Role | UniProt | Source organism | Status |
|---|---|---|---|
| Tryptophanase (TnaA) | P0A853 | *Escherichia coli* K-12 | Reviewed |
| Bacterial FMO (benchmark) | Q83XK4 | *Methylophaga thalassica* (called *M. aminisulfidivorans* in papers and PDB) | Unreviewed; crystal structures exist |
| FMO homolog | A0A1I3XW35 | *Methylophaga sulfidovorans* | Unreviewed |
| FMO homolog | A0A1E3GUV8 | *Methylophaga muralis* | Unreviewed |
| FMO homolog | A0ABN0TY91 | *Methylophaga marina* | Unreviewed; no AlphaFold DB model |
| FMO homolog | A0ABQ5U0Z0 | *Methylophaga thalassica* | Unreviewed; no AlphaFold DB model; 100% identical to the reference |
| Naphthalene dioxygenase alpha (NdoB) | P0A110 | *Pseudomonas putida* | Reviewed; structures exist |
| NdoB homolog | Q51494 | *Pseudomonas aeruginosa* | Reviewed |
| NdoB homolog | O07824 | *Pseudomonas fluorescens* | Reviewed; no AlphaFold DB model |
| NDO beta (NdoC), reference only | P0A112 | *Pseudomonas putida* | Reviewed |
| Styrene monooxygenase A (StyA) | O06834 | *Pseudomonas fluorescens* | Reviewed |
| Styrene monooxygenase B (StyB) | O06835 | *Pseudomonas fluorescens* | Reviewed |

Notes:

- F8U1Q5 and Q93LR7 were returned by UniProt for NdoB but are fragments (192 and 214 residues against 449 for the full-length protein). They are not used.
- NDO is a three-component system. Its ferredoxin and reductase partners are not identified in this project, so NdoB cannot be expected to make indigo on its own. It is kept out of the first wet-lab round.
- F5SYD3 has the same sequence as Q83XK4 and is not used.
- The FMO homologs are annotated in UniProt as trimethylamine monooxygenases. The benchmark carries the same annotation and is known to oxidise indole, but that is not established for the homologs.

## 4. Databases and tools

| Purpose | Resource |
|---|---|
| Sequences and annotation | UniProt (REST API) |
| Predicted structures | AlphaFold DB (API) |
| Experimental structures | PDB (Q83XK4: several entries including 2XVE, 2XVF, 2VQ7, 2VQB; P0A110: 1EG9, 1NDO) |
| Host metabolic model | BiGG Models, iML1515 |
| Sequence similarity | DIAMOND 2.1.9 |
| Flux balance analysis | COBRApy 0.29 |
| Codon optimisation | DNA Chisel 3.2.16 (with Biopython, python-codon-tables) |
| Workflow | Nextflow DSL2 (tested on 26.10.4), conda |

## 5. Pipeline overview

```
candidates_template.csv --> FETCH_CAND --+--> DIAMOND_REF --+
reference_enzymes.csv  --> FETCH_REF ----+                  +--> SCORE_RANK --> HANDOFF
                                         +--> STRUCT_QC ----+         |
iML1515.xml.gz --------> PATHWAY_FBA ----------------------------------+
constructs.csv --+                                                     |
                 +--> CONSTRUCT_DESIGN <-- FETCH_CAND fasta + ranking -+
```

| Step | Module | What it does |
|---|---|---|
| 1 | FETCH_CAND, FETCH_REF | Download protein sequences from UniProt for the candidate and reference CSVs. Rows with `TO_CURATE` are skipped. |
| 2 | DIAMOND_REF | BLASTP of each candidate against the reference set; keeps the best hit (identity, coverage). |
| 3 | STRUCT_QC | Download the AlphaFold DB model per candidate; report mean pLDDT. |
| 4 | PATHWAY_FBA | Add indole to indoxyl to indigo to iML1515; compute theoretical yields and a growth versus production envelope. |
| 5 | SCORE_RANK, HANDOFF | Rank candidates; write the wet-lab handoff (controls, acceptance criteria, a Latin-hypercube DoE matrix). |
| 6 | CONSTRUCT_DESIGN | Codon-optimise the shortlisted CDS for *E. coli* and build pET-style and Golden Gate cassettes with QC. |

### Scoring

```
score = w_identity * (identity/100 * coverage/100) + w_plddt * (mean pLDDT / 100)
```

The default weights are 0.5 and 0.5. They are placeholders and should be recalibrated once wet-lab titres exist. A missing AlphaFold model counts as pLDDT 0, which lowers the score without saying anything about the enzyme.

### FBA model

- Host model: iML1515; biomass reaction `BIOMASS_Ec_iML1515_core_75p37M`.
- Added reactions: `INDOLE_OXY` (indole + O2 + NADPH + H+ to indoxyl + NADP+ + H2O), `INDOXYL_DIMER` (2 indoxyl + O2 to indigo + 2 H2O), and a demand reaction `DM_indigo_c`. Mass and charge balance are checked in the script.
- Native tryptophanase (TRPAS2) is assumed present in the model.
- Scenarios: glucose only (10 mmol/gDW/h) and glucose plus tryptophan (up to 2 mmol/gDW/h), each at O2 limits of 5, 10, 20 and 1000 mmol/gDW/h.
- Envelope: maximum indigo while forcing 0 to 90% of maximum growth.

### Construct design

- Codon optimisation uses the *E. coli* usage table with codon-usage matching (not "always the most frequent codon"), which avoids repetitive sequences.
- Avoided in every CDS: BsaI, BsmBI, NdeI, XhoI sites; A/T runs of 7 or more; G/C runs of 6 or more; internal AGGAGG. GC kept within 0.40 to 0.65 over a 50 bp window (relaxed automatically, with a warning, if impossible).
- Two cassette designs per construct: **PET** (NdeI to XhoI; promoter and RBS come from the vector) and **GG** (BsaI Golden Gate with configurable 4-nt overhangs, default `AATG` and `GCTT`).
- Operons (StyA then StyB) are joined with `TA` plus a default RBS (`AAGGAGATATAAAA`).
- Point mutations (such as Y207S) are made by swapping a single codon in the optimised wild-type CDS, so the control differs from wild type by one codon. The script stops if the stated wild-type residue does not match.
- Final checks: translation identical to the protein, and restriction-site counts as expected in each cassette.

## 6. Results observed (run of 2026-11-03)

### 6.1 Smoke test

All seven analysis steps succeeded in about 3.5 minutes (excluding first-time conda environment builds). CONSTRUCT_DESIGN ran after the candidate file was corrected (see troubleshooting).

### 6.2 FBA (iML1515, theoretical maxima)

| Scenario | O2 limit | Max growth (1/h) | Max indigo (mmol/gDW/h) | Carbon yield |
|---|---|---|---|---|
| Glucose only | 5 | 0.354 | 1.178 | 0.314 |
| Glucose only | 10 | 0.532 | 1.774 | 0.473 |
| Glucose only | 20 | 0.822 | 2.572 | 0.686 |
| Glucose only | 1000 | 0.877 | 2.572 | 0.686 |
| Glucose + tryptophan | 5 | 0.373 | 1.667 | 0.325 |
| Glucose + tryptophan | 10 | 0.551 | 2.470 | 0.482 |
| Glucose + tryptophan | 20 | 0.849 | 3.460 | 0.675 |
| Glucose + tryptophan | 1000 | 0.957 | 3.719 | 0.726 |

Reading it:

- Glucose-only yield tops out at 0.257 mol indigo per mol glucose (69% carbon yield). The hard carbon limit is 0.375 (16 carbons in indigo, 6 in glucose), so the result is consistent.
- Glucose-only indigo rises with the oxygen limit and flattens at 20 mmol/gDW/h. In the tryptophan-fed scenario it keeps rising slightly beyond that. Both oxygenation steps need O2, so this route is oxygen-hungry. At fermenter scale the limit will be oxygen transfer.
- In the tryptophan-fed rows the column `mol_indigo_per_mol_glc` is misleading, because the indigo also draws on tryptophan carbon. Use `carbon_yield` for those rows.
- In the envelope, indigo falls roughly linearly from 2.57 at zero growth to about 0.52 at 80% of maximum growth. The product is not growth-coupled, so expect a growth phase followed by a production phase.
- These are upper bounds. FBA ignores feedback inhibition, indole toxicity, enzyme kinetics and by-products (indirubin, isatin). Use the numbers to compare scenarios, not to forecast titre.
- I did not check which route to indole the model uses in the glucose-only scenario (tryptophanase versus the native indole-3-glycerol phosphate reaction).

### 6.3 Candidate ranking

Rows marked "control" match their own reference at 100% by construction. They are not strong results.

| Candidate | Role | Identity (%) | pLDDT | Score | Flag |
|---|---|---|---|---|---|
| ctrl_bFMO_Q83XK4 | FMO benchmark (control) | 100.0 | 97.7 | 0.989 | - |
| cand_StyA_O06834 | StyA (control) | 100.0 | 97.3 | 0.987 | - |
| ctrl_ndoB_P0A110 | NdoB (control) | 100.0 | 97.1 | 0.986 | - |
| cand_tnaA_ecoli | TnaA (control) | 100.0 | 96.6 | 0.983 | - |
| cand_bFMO_A0A1I3XW35 | FMO, *M. sulfidovorans* | 97.1 | 97.6 | 0.974 | - |
| cand_StyB_O06835 | StyB (control) | 100.0 | 92.9 | 0.965 | - |
| cand_ndoB_Q51494 | NdoB, *P. aeruginosa* | 90.0 | 97.3 | 0.937 | - |
| cand_bFMO_A0A1E3GUV8 | FMO, *M. muralis* | 87.1 | 97.5 | 0.923 | - |
| cand_bFMO_A0ABQ5U0Z0 | FMO, *M. thalassica* | 100.0 | NA | 0.500 | no_af_model |
| cand_bFMO_A0ABN0TY91 | FMO, *M. marina* | 96.9 | NA | 0.485 | no_af_model |
| cand_ndoB_O07824 | NdoB, *P. fluorescens* | 95.8 | NA | 0.479 | no_af_model |

Reading it:

- The score measures similarity to known indigo-forming enzymes plus structure confidence. It does not predict whether a homolog oxidises indole.
- The three `no_af_model` rows score about 0.5 because the structure term is zero, not because they look poor. Their identity is as high as the others.
- The highest-ranked homologs are close variants of the benchmark. The 87% identity FMO is the most informative for diversity.

### 6.4 Construct designs

Five constructs, each in two designs (10 inserts total).

| Construct | Members | PET (bp) | GG (bp) | GC |
|---|---|---|---|---|
| FMO_Q83XK4 | Q83XK4 FMO | 1392 | 1402 | 0.481 |
| FMO_A0A1I3XW35 | FMO homolog (97%) | 1392 | 1402 | 0.483 |
| FMO_A0A1E3GUV8 | FMO homolog (87%) | 1392 | 1402 | 0.485 |
| SMO_StyAB | StyA then StyB operon | 1798 | 1808 | 0.520 |
| NEG_Q83XK4_Y207S | Q83XK4 with Y207S | 1392 | 1402 | 0.483 |

- All six CDS translate back to the input proteins; CAI is 0.69 to 0.73 and GC is 0.485 to 0.532.
- All inserts have the expected site counts: PET has one NdeI and one XhoI; GG has two BsaI; none has BsmBI.
- Residue 207 of the Q83XK4 UniProt sequence is a tyrosine (context `...LVGSSYSAEDI...`), so the Y207S control builds. The source paper's numbering still has to be confirmed (see section 8).

### 6.5 Literature context (from a literature search, verify before citing)

- FMO from *Methylophaga* in *E. coli*: reported up to about 920 mg/L indigo in a 5 L fermenter from tryptophan medium, with scale-up work to 3000 L. Indirubin versus indigo depends strongly on culture conditions.
- StyAB in *E. coli*: about 566 mg/L indigo from 1.2 mg/mL tryptophan.
- De novo from glucose in *E. coli* (Pham et al. 2024): about 965 mg/L in a 3 L fed-batch, using an eight-gene design that includes a plant indole lyase.

## 7. How to reproduce

### Requirements

- A GitHub Codespace (or any Linux machine) with internet access, Java 17 or newer, conda, and about 2 CPUs / 8 GB RAM.
- The iML1515 SBML from BiGG.

### Steps

1. Create a GitHub repository, open a Codespace, and upload `indigo-nf.zip`.

```bash
unzip -o indigo-nf.zip && cd indigo-nf
chmod +x bin/*.py tools/*.py
```

2. Install Nextflow (skip if `nextflow -version` already works).

```bash
java -version
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

3. Configure conda channels (the modules pull from bioconda and conda-forge).

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

4. Download the host model.

```bash
wget http://bigg.ucsd.edu/static/models/iML1515.xml.gz
```

5. Run.

```bash
nextflow run main.nf -profile conda --model iML1515.xml.gz --outdir results
```

The first run is slow because conda builds one environment per tool (DIAMOND, COBRApy, DNA Chisel). Add `-resume` on later runs.

6. Look at the results.

```bash
cat results/03_fba/fba_summary.tsv
cat results/04_ranking/ranked_candidates.tsv
cat results/05_handoff/handoff.md
cat results/06_constructs/construct_design.md
```

### Find new candidates

```bash
python3 tools/find_candidates.py --top 3 --size 15
```

This queries UniProt and writes `candidates_review.tsv` (every hit, with organism, reviewed status and PDB links) and `candidates_draft.csv`. Check protein names and organisms, then copy chosen rows into `assets/candidates_template.csv` and rerun with `-resume`. Extra searches can be added with `--extra "label|role|UniProt query"`.

### Input files

| File | Columns | Purpose |
|---|---|---|
| `assets/candidates_template.csv` | id, uniprot_acc, role | Enzymes to evaluate. Rows marked `ctrl_` are also in the reference set. |
| `assets/reference_enzymes.csv` | id, uniprot_acc, role | Characterised enzymes to compare against. |
| `assets/constructs.csv` | construct_id, members, mutations, notes | Which candidates go into which construct. Members are candidate ids joined by `;`. Mutations use `member_id:Y207S`, UniProt 1-based numbering. |

### Parameters (in `nextflow.config`)

`--model`, `--outdir`, `--biomass_id`, `--top_n`, `--w_identity`, `--w_plddt`, `--doe_runs`, `--constructs`, `--rbs`, `--gg_oh5`, `--gg_oh3`, `--gc_min`, `--gc_max`, `--seed`.

### Project layout

```
indigo-nf/
  main.nf
  nextflow.config
  assets/        candidates_template.csv, reference_enzymes.csv, constructs.csv
  bin/           fetch_uniprot.py, struct_qc.py, pathway_fba.py, score_rank.py,
                 handoff.py, design_constructs.py
  modules/local/ one Nextflow process per step
  envs/          construct.yml (DNA Chisel, pip-installed inside conda)
  tools/         find_candidates.py
results/
  01_similarity  02_structure  03_fba  04_ranking  05_handoff  06_constructs
  pipeline_info  (timeline, report, trace)
```

### Reproducibility notes

- UniProt and AlphaFold DB are queried live, so results can change as the databases change. Record the run date with your results.
- Codon optimisation is seeded (default 1) and was checked to give identical output on repeated runs.
- Pin tool versions in the module `conda` directives if you need a locked environment. DIAMOND 2.1.9, COBRApy 0.29 and DNA Chisel 3.2.16 are pinned now; Python and the other dependencies are not.
- Tested on Nextflow 26.10.4. This version enforces strict syntax, which the scripts here follow.

### Troubleshooting (problems met while building this)

| Symptom | Cause | Fix |
|---|---|---|
| `Statements cannot be mixed with script declarations` | Strict syntax in Nextflow 26.x | Use the shipped `main.nf` (check is inside `workflow {}`). |
| Exit status 126, `Permission denied` on a script | Executable bit lost on upload | `chmod +x bin/*.py tools/*.py` |
| `Properties init: Could not determine current working directory` | Terminal sat in a deleted folder | `cd ~`, then `cd` back into the project. |
| `no_ref_hit` for most candidates | `TO_CURATE` still in the reference CSV | Fill in real accessions. |
| `member '...' not in the candidate FASTA` | IDs in `constructs.csv` do not match `candidates_template.csv` | Make the IDs identical; the error lists the IDs found. |
| Screen of `~` and `(END)` | Output opened in the `less` pager | Press `q`. |

## 8. Before sending anything to the wet lab

Computational designs only. Check the following first:

- **Y207S numbering.** The script uses UniProt numbering. Confirm that position 207 is the intended residue in the source paper (Cho et al. 2011) and its Y207S structure entries in the PDB.
- **RBS and overhangs** are defaults. Match them to your vector or toolkit.
- **Tags.** The pET design ends in a stop codon (no C-terminal tag). N-terminal tags depend on the backbone. Tags can change FMO activity.
- **Promoter, backbone and RBS strength** are not designed here.
- **StyA and StyB** are both from *P. fluorescens* and must be co-expressed.
- **tnaA** is assumed native to the host. Add a row to `constructs.csv` to overexpress it.
- **Freedom to operate.** Patents on microbial indigo exist. Check before commercial use. For work in Germany, contained use of genetically modified organisms is regulated (Gentechnikgesetz).

### Suggested first round

1. Q83XK4 FMO as the benchmark.
2. A0A1I3XW35 FMO (97% identity) and A0A1E3GUV8 FMO (87% identity) as variants.
3. StyA + StyB as a second system.
4. Controls: empty vector, TnaA only, Q83XK4 Y207S mutant (reported to lack indole oxygenation), tryptophan-fed versus unfed.
5. Leave NdoB out until its partner proteins are identified.

Readouts: indigo by A600 in DMSO after extraction, confirmed by HPLC-DAD; LC-MS [M+H]+ m/z 263.0815 (C16H10N2O2); indirubin (about 540 nm) to indigo ratio; at least three biological replicates. Minimum titre and purity thresholds are still to be defined by the lab.

The DoE matrix in `results/05_handoff/doe_matrix.csv` is a Latin hypercube over temperature (25-37 C), pH (6.5-7.5), tryptophan (0-2 g/L) and IPTG (0.05-0.5 mM). These are starting ranges, not validated conditions.

## 9. Known gaps and next steps

- Add the indole-lyase route to the FBA model, and report flux through each indole-supplying reaction.
- RBS strength and 5' mRNA structure prediction for each construct.
- Docking of indole into Q83XK4 and homologs (experimental structures are available in the PDB).
- Identify the ferredoxin and reductase for NDO, or drop NDO.
- Predict structures with ColabFold for the three accessions without an AlphaFold DB model.
- Growth-coupled knockout design (OptKnock or minimal cut sets).
- Calibrate the scoring weights with wet-lab titres, then feed results back into the model.

## 10. References

- Cho H.J. et al. (2011). Structural and functional analysis of bacterial flavin-containing monooxygenase reveals its ping-pong-type reaction mechanism. J Struct Biol 175:39.
- Alfieri A. et al. (2008). Revealing the moonlighting role of NADP in the structure of a flavin-containing monooxygenase. PNAS 105:6572.
- Carredano E. et al. (2000). Substrate binding site of naphthalene 1,2-dioxygenase: functional implications of indole binding. J Mol Biol 296:701.
- Han G.H. et al. Bio-indigo production with recombinant E. coli harbouring a flavin-containing monooxygenase gene (fmo). Chosun University; see J Biotechnol 150 (2010) and Process Biochem.
- Du L., Yue J., Zhu Y., Yin S. (2022). Production of indigo by recombinant Escherichia coli with expression of monooxygenase, tryptophanase, and molecular chaperone. Foods 11:2117.
- Pham N.N. et al. (2024). Auto-inducible synthetic pathway in E. coli enhanced sustainable indigo production from glucose. Metab Eng 85:14-25.
- Du J., Yang D., Luo Z.W., Lee S.Y. (2018). Metabolic engineering of Escherichia coli for the production of indirubin from glucose. J Biotechnol 267:19-28.
- Murdock D. et al. (1993). Construction of metabolic operons catalyzing the de novo biosynthesis of indigo in Escherichia coli. Bio/Technology 11:381.
- Monk J.M. et al. (2017). iML1515, a knowledgebase that computes Escherichia coli traits. Nat Biotechnol 35:904.
- Zulkower V., Rosser S. (2020). DNA Chisel, a versatile sequence optimizer. Bioinformatics 36:4508.
- Buchfink B., Reuter K., Drost H.G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat Methods 18:366.
- Ebrahim A. et al. (2013). COBRApy: constraints-based reconstruction and analysis for Python. BMC Syst Biol 7:74.

Literature figures in section 6.5 come from search snippets and abstracts; check the original papers before quoting them. Reference details were compiled from search results and memory, so verify them before citing.
