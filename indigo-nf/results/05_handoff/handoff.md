# Indigo - wet-lab handoff

## Top candidates (in silico ranking)

| rank | id | UniProt | score | pLDDT | best ref | flags |
|---|---|---|---|---|---|---|
| 1 | ctrl_bFMO_Q83XK4 | Q83XK4 | 0.9885 | 97.7 | ref_bFMO_maminisulfidivorans | - |
| 2 | cand_StyA_O06834 | O06834 | 0.9865 | 97.3 | ref_styA_pfluorescens | - |
| 3 | ctrl_ndoB_P0A110 | P0A110 | 0.9855 | 97.1 | ref_ndoB_pseudomonas_9816 | - |
| 4 | cand_tnaA_ecoli | P0A853 | 0.983 | 96.6 | ref_tnaA_ecoli | - |
| 5 | cand_bFMO_A0A1I3XW35 | A0A1I3XW35 | 0.9735 | 97.6 | ref_bFMO_maminisulfidivorans | - |

## Pathway feasibility (FBA, host GEM)

```
scenario	o2_limit	max_growth	max_indigo	mol_indigo_per_mol_glc	carbon_yield
glucose_only	5	0.3535	1.1780	0.1178	0.3141
glucose_only	10	0.5315	1.7741	0.1774	0.4731
glucose_only	20	0.8218	2.5721	0.2572	0.6859
glucose_only	1000	0.8770	2.5721	0.2572	0.6859
glucose_plus_trp	5	0.3727	1.6667	0.1667	0.3252
glucose_plus_trp	10	0.5513	2.4696	0.2470	0.4819
glucose_plus_trp	20	0.8487	3.4601	0.3460	0.6751
glucose_plus_trp	1000	0.9571	3.7186	0.3719	0.7256
```

## Controls

- Empty-vector host
- Host + tnaA only (no oxygenase)
- No-tryptophan vs tryptophan-fed
- Indigo standard for calibration

## Acceptance criteria (set thresholds before testing)

- Indigo: A600 in DMSO after extraction, confirmed by HPLC-DAD
- LC-MS: indigo [M+H]+ m/z 263.0815 (C16H10N2O2)
- Indirubin (about 540 nm) and isatin as by-products; report indirubin/indigo ratio
- Minimum titre / purity: TO BE DEFINED
- Biological replicates: n>=3

## DoE

See doe_matrix.csv (starting ranges; adjust to what your lab can run).
