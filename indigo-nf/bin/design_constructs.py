#!/usr/bin/env python3
"""Shortlist -> codon-optimised E. coli CDS -> construct designs.

Inputs
  --fasta        protein FASTA from the pipeline (headers: >id|acc|role)
  --ranked       ranked_candidates.tsv (used only to annotate warnings)
  --constructs   CSV: construct_id,members,mutations,notes
                   members   : candidate ids joined by ';' (operon order)
                   mutations : 'member_id:Y207S' joined by ';' (UniProt numbering, 1-based)

Outputs (in --outdir)
  shortlist_proteins.fasta, shortlist_codon_opt.fasta, order_inserts.fasta,
  construct_summary.tsv, qc_report.tsv, construct_design.md, constructs/*.gb

Two cassette designs are produced per construct:
  PET : NdeI ... XhoI (pET-style; RBS and promoter come from the vector)
  GG  : BsaI Golden Gate with configurable 4-nt overhangs
Defaults are starting points; match overhangs and RBS to your own toolkit/vector.
"""
import argparse, csv, math, os, random, sys
import numpy as np
import python_codon_tables as pct
from Bio.Restriction import BsaI, BsmBI, NdeI, XhoI
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from dnachisel import (AvoidPattern, CodonOptimize, DnaOptimizationProblem,
                       EnforceGCContent, EnforceTranslation, reverse_translate)

TABLE = pct.get_codons_table("e_coli")
W = {c: f / max(d.values()) for aa, d in TABLE.items() for c, f in d.items()}
AA_CODONS = {aa: sorted(d, key=lambda c: -d[c]) for aa, d in TABLE.items()}
ENZ = {"NdeI": NdeI, "XhoI": XhoI, "BsaI": BsaI, "BsmBI": BsmBI}
AVOID_PATTERNS = ["AAAAAAA", "TTTTTTT", "GGGGGG", "CCCCCC", "AGGAGG"]  # homopolymers, internal SD
AVOID_ENZ = ["BsaI_site", "BsmBI_site", "NdeI_site", "XhoI_site"]
PAD5, PAD3 = "TTAATT", "AATTAA"  # leading/trailing bases so enzymes cut efficiently


def rc(s):
    return str(Seq(s).reverse_complement())


def gc(s):
    return (s.count("G") + s.count("C")) / len(s)


def cai(cds):
    cod = [cds[i:i + 3] for i in range(0, len(cds) - 3, 3)]
    ws = [W[c] for c in cod if c not in ("ATG", "TGG") and c in W]
    return math.exp(sum(math.log(w) for w in ws) / len(ws))


def max_run(s):
    best = run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def site_counts(s):
    return {n: len(e.search(Seq(s))) for n, e in ENZ.items()}


def forbidden_hits(cds):
    hits = [n for n, k in site_counts(cds).items() if k]
    hits += [p for p in AVOID_PATTERNS if p in cds or p in rc(cds)]
    return hits


def optimise(protein, seed, gc_min, gc_max):
    """Return (cds, note). Relaxes the GC window step by step if it cannot be met."""
    from dnachisel.DnaOptimizationProblem.NoSolutionError import NoSolutionError
    windows = [(gc_min, gc_max), (gc_min - 0.05, gc_max + 0.05), (0.30, 0.75)]
    for i, (lo, hi) in enumerate(windows):
        random.seed(seed)
        np.random.seed(seed)
        start = reverse_translate(protein + "*")
        cons = [EnforceTranslation(), EnforceGCContent(mini=lo, maxi=hi, window=50)]
        cons += [AvoidPattern(p) for p in AVOID_ENZ + AVOID_PATTERNS]
        prob = DnaOptimizationProblem(sequence=start, constraints=cons,
                                      objectives=[CodonOptimize(species="e_coli", method="match_codon_usage")],
                                      logger=None)
        prob.max_random_iters = 5000
        try:
            prob.resolve_constraints()
        except NoSolutionError:
            continue
        prob.optimize()
        cds = prob.sequence
        if str(Seq(cds).translate()) != protein + "*":
            sys.exit("internal error: translation mismatch after optimisation")
        note = "" if i == 0 else f"GC window relaxed to {lo:.2f}-{hi:.2f}"
        return cds, note
    sys.exit("could not satisfy constraints even with the widest GC window")


def mutate_codon(wt_cds, pos, wt_aa, new_aa):
    """Swap one codon in an optimised CDS (minimal change vs the wild-type construct)."""
    cur = str(Seq(wt_cds[3 * (pos - 1): 3 * pos]).translate())
    if cur != wt_aa:
        sys.exit(f"mutation {wt_aa}{pos}{new_aa}: residue {pos} is {cur}, not {wt_aa}")
    for c in AA_CODONS[new_aa]:
        cand = wt_cds[:3 * (pos - 1)] + c + wt_cds[3 * pos:]
        if not forbidden_hits(cand):
            return cand
    sys.exit(f"mutation {wt_aa}{pos}{new_aa}: no codon avoids forbidden sites")


def read_fasta(path):
    out = {}
    for rec in SeqIO.parse(path, "fasta"):
        cid, acc, role = (rec.id.split("|") + ["", ""])[:3]
        out[cid] = (acc, role, str(rec.seq))
    return out


def build(mode, cdss, rbs, oh5, oh3):
    """Return (sequence, features[(start,end,type,label)])."""
    feats, seq = [], ""

    def add(s, typ=None, label=None):
        nonlocal seq
        if typ:
            feats.append((len(seq), len(seq) + len(s), typ, label))
        seq += s

    if mode == "PET":
        add(PAD5, "misc_feature", "5' pad")
        add("CAT", "misc_feature", "NdeI (CATATG; ATG is the start codon)")
        body_first = cdss[0][1]
    else:
        add(PAD5, "misc_feature", "5' pad")
        add("GGTCTC", "misc_feature", "BsaI")
        add("A", "misc_feature", "N1")
        add(oh5, "misc_feature", "5' overhang")
        body_first = cdss[0][1][3:] if oh5.endswith("ATG") else cdss[0][1]
    for i, (label, cds) in enumerate(cdss):
        if i == 0:
            add(body_first, "CDS", label)
        else:
            add("TA", "misc_feature", "spacer")
            add(rbs, "regulatory", "RBS")
            add(cds, "CDS", label)
    if mode == "PET":
        add("CTCGAG", "misc_feature", "XhoI")
    else:
        add(oh3, "misc_feature", "3' overhang")
        add("T", "misc_feature", "N1")
        add("GAGACC", "misc_feature", "BsaI (reverse)")
    add(PAD3, "misc_feature", "3' pad")
    return seq, feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--ranked", required=True)
    ap.add_argument("--constructs", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--rbs", default="AAGGAGATATAAAA")
    ap.add_argument("--oh5", default="AATG")
    ap.add_argument("--oh3", default="GCTT")
    ap.add_argument("--gc-min", type=float, default=0.40)
    ap.add_argument("--gc-max", type=float, default=0.65)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(f"{a.outdir}/constructs", exist_ok=True)

    prots = read_fasta(a.fasta)
    ranked = {r["id"]: r for r in csv.DictReader(open(a.ranked), delimiter="\t")}
    rows = list(csv.DictReader(open(a.constructs)))
    warnings = []

    # -- parse constructs ------------------------------------------------
    needed, muts = [], {}
    for r in rows:
        members = [m.strip() for m in r["members"].split(";") if m.strip()]
        for m in members:
            if m not in prots:
                sys.exit(f"construct {r['construct_id']}: member '{m}' not in the candidate FASTA")
            if m not in needed:
                needed.append(m)
            fl = ranked.get(m, {}).get("flags", "-")
            if fl not in ("-", ""):
                warnings.append(f"{m}: pipeline flags = {fl}")
        r["_members"] = members
        r["_muts"] = {}
        for item in filter(None, (r.get("mutations") or "").split(";")):
            mid, mut = item.split(":")
            if mid not in members:
                sys.exit(f"construct {r['construct_id']}: mutation target {mid} not a member")
            r["_muts"][mid] = mut.strip()

    # -- optimise wild-type CDS once per member --------------------------
    cds_wt, prot_out, cds_out, qc = {}, {}, {}, []
    for m in needed:
        p = prots[m][2]
        if not p.startswith("M"):
            warnings.append(f"{m}: protein does not start with Met")
        cds_wt[m], note = optimise(p, a.seed, a.gc_min, a.gc_max)
        if note:
            warnings.append(f"{m}: {note}")
        prot_out[m], cds_out[m] = p, cds_wt[m]

    # -- variants ---------------------------------------------------------
    var_cds, var_prot = {}, {}
    for r in rows:
        for m, mut in r["_muts"].items():
            wt_aa, pos, new_aa = mut[0], int(mut[1:-1]), mut[-1]
            p = prots[m][2]
            if p[pos - 1] != wt_aa:
                ctx = p[max(0, pos - 6): pos + 5]
                sys.exit(f"{r['construct_id']}: {mut}: residue {pos} of {m} is {p[pos - 1]} "
                         f"(context ...{ctx}...). Check numbering against the paper.")
            key = f"{m}_{mut}"
            var_cds[key] = mutate_codon(cds_wt[m], pos, wt_aa, new_aa)
            var_prot[key] = p[:pos - 1] + new_aa + p[pos:]
            prot_out[key], cds_out[key] = var_prot[key], var_cds[key]
            warnings.append(f"{key}: numbering is UniProt 1-based; context ...{p[max(0, pos - 6): pos + 5]}... "
                            f"- confirm it matches the position intended in the source paper.")

    # -- QC per CDS --------------------------------------------------------
    for k, cds in cds_out.items():
        prot = prot_out[k]
        qc.append({"item": k, "kind": "CDS", "aa": len(prot), "nt": len(cds), "gc": f"{gc(cds):.3f}",
                   "cai": f"{cai(cds):.3f}", "max_homopolymer": max_run(cds),
                   "forbidden": ",".join(forbidden_hits(cds)) or "-",
                   "translation_ok": str(str(Seq(cds).translate()) == prot + "*")})

    # -- constructs --------------------------------------------------------
    inserts, summary = [], []
    for r in rows:
        cid = r["construct_id"]
        cdss = []
        for m in r["_members"]:
            k = f"{m}_{r['_muts'][m]}" if m in r["_muts"] else m
            cdss.append((k, cds_out[k]))
        for mode in ("PET", "GG"):
            seq, feats = build(mode, cdss, a.rbs, a.oh5, a.oh3)
            cnt = site_counts(seq)
            exp = {"PET": {"NdeI": 1, "XhoI": 1, "BsaI": 0, "BsmBI": 0},
                   "GG": {"NdeI": 0, "XhoI": 0, "BsaI": 2, "BsmBI": 0}}[mode]
            ok = cnt == exp
            if not ok:
                warnings.append(f"{cid} {mode}: site counts {cnt} differ from expected {exp}")
            rec = SeqRecord(Seq(seq), id=f"{cid}_{mode}", name=f"{cid}_{mode}"[:16],
                            description=f"{cid} {mode} cassette (design, not verified in the lab)",
                            annotations={"molecule_type": "DNA", "topology": "linear"})
            for s, e, typ, label in feats:
                q = {"label": [label]}
                if typ == "CDS":
                    q.update({"codon_start": [1], "transl_table": [11],
                              "translation": [prot_out[label]]})
                rec.features.append(SeqFeature(FeatureLocation(s, e, strand=1), type=typ, qualifiers=q))
            SeqIO.write(rec, f"{a.outdir}/constructs/{cid}.{mode}.gb", "genbank")
            inserts.append((f"{cid}|{mode}|{len(seq)}bp", seq))
            qc.append({"item": f"{cid}.{mode}", "kind": "insert", "aa": "-", "nt": len(seq),
                       "gc": f"{gc(seq):.3f}", "cai": "-", "max_homopolymer": max_run(seq),
                       "forbidden": "-", "translation_ok": f"sites_ok={ok} {cnt}"})
            summary.append({"construct": cid, "design": mode, "members": ";".join(k for k, _ in cdss),
                            "length_bp": len(seq), "gc": f"{gc(seq):.3f}", "sites_as_expected": ok,
                            "notes": r.get("notes", "")})

    # -- write outputs -----------------------------------------------------
    with open(f"{a.outdir}/shortlist_proteins.fasta", "w") as o:
        for k, p in prot_out.items():
            o.write(f">{k}\n{p}\n")
    with open(f"{a.outdir}/shortlist_codon_opt.fasta", "w") as o:
        for k, c in cds_out.items():
            o.write(f">{k}|E_coli_optimised\n{c}\n")
    with open(f"{a.outdir}/order_inserts.fasta", "w") as o:
        for h, s in inserts:
            o.write(f">{h}\n{s}\n")
    for name, data, cols in (("construct_summary.tsv", summary, list(summary[0])),
                             ("qc_report.tsv", qc, list(qc[0]))):
        with open(f"{a.outdir}/{name}", "w", newline="") as o:
            w = csv.DictWriter(o, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(data)

    with open(f"{a.outdir}/construct_design.md", "w") as o:
        o.write("# Construct design\n\n")
        o.write("Computational designs only. Nothing here has been synthesised or tested.\n\n")
        o.write("## Parameters\n\n")
        o.write(f"- Codon optimisation: DNA Chisel, E. coli usage table, match_codon_usage, seed {a.seed}\n")
        o.write(f"- GC window: {a.gc_min}-{a.gc_max} (50 bp window)\n")
        o.write("- Avoided: BsaI, BsmBI, NdeI, XhoI sites; A/T runs of 7, G/C runs of 6; internal AGGAGG\n")
        o.write(f"- Operon spacer: TA + RBS `{a.rbs}` (default; tune with an RBS calculator)\n")
        o.write(f"- Golden Gate overhangs: 5' `{a.oh5}`, 3' `{a.oh3}` (defaults; match your toolkit)\n")
        o.write("- PET design assumes the vector supplies promoter and RBS upstream of NdeI\n\n")
        o.write("## Constructs\n\n| construct | design | members | length (bp) | GC | sites as expected |\n|---|---|---|---|---|---|\n")
        for s in summary:
            o.write(f"| {s['construct']} | {s['design']} | {s['members']} | {s['length_bp']} | {s['gc']} | {s['sites_as_expected']} |\n")
        o.write("\n## Warnings and things to check\n\n")
        for w in sorted(set(warnings)) or ["none"]:
            o.write(f"- {w}\n")
        o.write("\n## Not covered yet\n\n- RBS strength / 5' mRNA structure, tags, plasmid backbone, promoter choice\n")
        o.write("- tnaA is assumed to be native to the host; add a construct row if you want it overexpressed\n")
        o.write("- Order synthesis only after checking the sequences against your vector and toolkit\n")
    print(f"{len(cds_out)} CDS, {len(inserts)} inserts written")


if __name__ == "__main__":
    main()
