# Construct design

Computational designs only. Nothing here has been synthesised or tested.

## Parameters

- Codon optimisation: DNA Chisel, E. coli usage table, match_codon_usage, seed 1
- GC window: 0.4-0.65 (50 bp window)
- Avoided: BsaI, BsmBI, NdeI, XhoI sites; A/T runs of 7, G/C runs of 6; internal AGGAGG
- Operon spacer: TA + RBS `AAGGAGATATAAAA` (default; tune with an RBS calculator)
- Golden Gate overhangs: 5' `AATG`, 3' `GCTT` (defaults; match your toolkit)
- PET design assumes the vector supplies promoter and RBS upstream of NdeI

## Constructs

| construct | design | members | length (bp) | GC | sites as expected |
|---|---|---|---|---|---|
| FMO_Q83XK4 | PET | ctrl_bFMO_Q83XK4 | 1392 | 0.481 | True |
| FMO_Q83XK4 | GG | ctrl_bFMO_Q83XK4 | 1402 | 0.481 | True |
| FMO_A0A1I3XW35 | PET | cand_bFMO_A0A1I3XW35 | 1392 | 0.483 | True |
| FMO_A0A1I3XW35 | GG | cand_bFMO_A0A1I3XW35 | 1402 | 0.483 | True |
| FMO_A0A1E3GUV8 | PET | cand_bFMO_A0A1E3GUV8 | 1392 | 0.485 | True |
| FMO_A0A1E3GUV8 | GG | cand_bFMO_A0A1E3GUV8 | 1402 | 0.485 | True |
| SMO_StyAB | PET | cand_StyA_O06834;cand_StyB_O06835 | 1798 | 0.520 | True |
| SMO_StyAB | GG | cand_StyA_O06834;cand_StyB_O06835 | 1808 | 0.520 | True |
| NEG_Q83XK4_Y207S | PET | ctrl_bFMO_Q83XK4_Y207S | 1392 | 0.483 | True |
| NEG_Q83XK4_Y207S | GG | ctrl_bFMO_Q83XK4_Y207S | 1402 | 0.483 | True |

## Warnings and things to check

- ctrl_bFMO_Q83XK4_Y207S: numbering is UniProt 1-based; context ...LVGSSYSAEDI... - confirm it matches the position intended in the source paper.

## Not covered yet

- RBS strength / 5' mRNA structure, tags, plasmid backbone, promoter choice
- tnaA is assumed to be native to the host; add a construct row if you want it overexpressed
- Order synthesis only after checking the sequences against your vector and toolkit
