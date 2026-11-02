#!/usr/bin/env python3
"""Rank candidate enzymes: score = w_id * (identity * query coverage to best reference hit)
+ w_pl * mean pLDDT/100. Weights are placeholders until wet-lab data calibrates them."""
import csv, sys

hits_f, struct_f, out_f, w_id, w_pl = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])

best = {}
for line in open(hits_f):
    if not line.strip():
        continue
    q, s, pid, length, ev, bs, qcov = line.rstrip("\n").split("\t")
    if q not in best or float(bs) > best[q][4]:
        best[q] = (s, float(pid), float(qcov), float(ev), float(bs))

rows = []
for r in csv.DictReader(open(struct_f), delimiter="\t"):
    key = "|".join([r["id"], r["acc"], r["role"]])
    h = best.get(key)
    id_score = (h[1] / 100) * (h[2] / 100) if h else 0.0
    pl = float(r["mean_plddt"]) / 100 if r["mean_plddt"] != "NA" else 0.0
    flags = []
    if not h:
        flags.append("no_ref_hit")
    if r["status"] != "ok":
        flags.append("no_af_model")
    elif pl < 0.70:
        flags.append("low_plddt")
    rows.append({"id": r["id"], "acc": r["acc"], "role": r["role"],
                 "best_ref": h[0].split("|")[0] if h else "NA",
                 "pident": f"{h[1]:.1f}" if h else "NA", "qcov": f"{h[2]:.1f}" if h else "NA",
                 "mean_plddt": r["mean_plddt"], "score": round(w_id * id_score + w_pl * pl, 4),
                 "flags": ",".join(flags) or "-"})

rows.sort(key=lambda x: -x["score"])
with open(out_f, "w", newline="") as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0].keys()) if rows else ["id"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)
