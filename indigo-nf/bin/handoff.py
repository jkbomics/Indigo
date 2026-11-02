#!/usr/bin/env python3
"""Build the wet-lab handoff: handoff.md (candidates, in silico summary, acceptance criteria)
and doe_matrix.csv (Latin hypercube over starting ranges - adjust to your lab)."""
import csv, random, sys

ranked_f, fba_f, top_n, runs = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
ranked = list(csv.DictReader(open(ranked_f), delimiter="\t"))[:top_n]
fba = open(fba_f).read().strip().splitlines()

random.seed(1)
factors = {"temp_C": (25, 37), "pH": (6.5, 7.5), "trp_g_per_L": (0, 2), "IPTG_mM": (0.05, 0.5)}


def lhs(n, lo, hi):
    pts = [(i + random.random()) / n for i in range(n)]
    random.shuffle(pts)
    return [lo + p * (hi - lo) for p in pts]


cols = {k: lhs(runs, *v) for k, v in factors.items()}
with open("doe_matrix.csv", "w", newline="") as o:
    w = csv.writer(o)
    w.writerow(["run"] + list(factors))
    for i in range(runs):
        w.writerow([i + 1] + [f"{cols[k][i]:.2f}" for k in factors])

with open("handoff.md", "w") as o:
    o.write("# Indigo - wet-lab handoff\n\n## Top candidates (in silico ranking)\n\n")
    o.write("| rank | id | UniProt | score | pLDDT | best ref | flags |\n|---|---|---|---|---|---|---|\n")
    for i, r in enumerate(ranked, 1):
        o.write(f"| {i} | {r['id']} | {r['acc']} | {r['score']} | {r['mean_plddt']} | {r['best_ref']} | {r['flags']} |\n")
    o.write("\n## Pathway feasibility (FBA, host GEM)\n\n```\n" + "\n".join(fba) + "\n```\n")
    o.write("\n## Controls\n\n- Empty-vector host\n- Host + tnaA only (no oxygenase)\n- No-tryptophan vs tryptophan-fed\n- Indigo standard for calibration\n")
    o.write("\n## Acceptance criteria (set thresholds before testing)\n\n")
    o.write("- Indigo: A600 in DMSO after extraction, confirmed by HPLC-DAD\n")
    o.write("- LC-MS: indigo [M+H]+ m/z 263.0815 (C16H10N2O2)\n")
    o.write("- Indirubin (about 540 nm) and isatin as by-products; report indirubin/indigo ratio\n")
    o.write("- Minimum titre / purity: TO BE DEFINED\n- Biological replicates: n>=3\n")
    o.write("\n## DoE\n\nSee doe_matrix.csv (starting ranges; adjust to what your lab can run).\n")
