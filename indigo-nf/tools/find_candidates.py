#!/usr/bin/env python3
"""Query UniProt for candidate indigo-forming enzymes and write
  candidates_review.tsv : every hit, with label, reviewed status, EC, PDB links (read this)
  candidates_draft.csv  : top hits per label in the id,uniprot_acc,role format the pipeline expects

Nothing here is auto-trusted: open the review table, check protein name and organism,
then copy the rows you want into assets/candidates_template.csv.

Usage:
  python3 tools/find_candidates.py
  python3 tools/find_candidates.py --top 3 --size 15
  python3 tools/find_candidates.py --extra "P450|p450_oxygenase|ec:1.14.14.1 AND organism_name:Bacillus"
"""
import argparse, csv, io, sys, urllib.parse, urllib.request

BASE = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,protein_name,gene_names,organism_name,length,reviewed,ec,xref_pdb"

# label -> (role, UniProt query). Edit freely; syntax: https://www.uniprot.org/help/query-fields
DEFAULT = {
    "bFMO":   ("flavin_monooxygenase",            "ec:1.14.13.148 AND organism_name:Methylophaga"),
    "NdoB":   ("naphthalene_dioxygenase_alpha",   "gene:ndoB AND organism_name:Pseudomonas"),
    "StyA":   ("styrene_monooxygenase_A",         "gene:styA AND organism_name:Pseudomonas"),
    "StyB":   ("styrene_monooxygenase_B",         "gene:styB AND organism_name:Pseudomonas"),
    "TnaA":   ("tryptophanase",                   "gene:tnaA AND organism_id:83333"),
}


def search(query, size):
    params = urllib.parse.urlencode({"query": query, "format": "tsv", "fields": FIELDS, "size": size})
    with urllib.request.urlopen(f"{BASE}?{params}", timeout=60) as r:
        return r.read().decode()


def parse(tsv_text):
    return list(csv.DictReader(io.StringIO(tsv_text), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", type=int, default=10, help="max hits per query")
    ap.add_argument("--top", type=int, default=2, help="rows per label in the draft CSV")
    ap.add_argument("--extra", action="append", default=[], help='"label|role|uniprot query" (repeatable)')
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    queries = dict(DEFAULT)
    for e in a.extra:
        label, role, q = e.split("|", 2)
        queries[label] = (role, q)

    review, draft = [], []
    for label, (role, q) in queries.items():
        try:
            rows = parse(search(q, a.size))
        except Exception as ex:
            print(f"WARN {label}: query failed ({ex})", file=sys.stderr)
            continue
        if not rows:
            print(f"WARN {label}: no hits for: {q}", file=sys.stderr)
            continue
        # reviewed (Swiss-Prot) entries first, then as UniProt returned them
        rows.sort(key=lambda r: r.get("Reviewed", "") != "reviewed")
        for r in rows:
            review.append({"label": label, "role": role, "acc": r["Entry"], "name": r["Protein names"],
                           "gene": r["Gene Names"], "organism": r["Organism"], "length": r["Length"],
                           "reviewed": r["Reviewed"], "ec": r["EC number"], "pdb": r["PDB"], "query": q})
        for r in rows[:a.top]:
            draft.append({"id": f"cand_{label}_{r['Entry']}", "uniprot_acc": r["Entry"], "role": role})

    cols = ["label", "role", "acc", "name", "gene", "organism", "length", "reviewed", "ec", "pdb", "query"]
    with open(f"{a.out}/candidates_review.tsv", "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(review)
    with open(f"{a.out}/candidates_draft.csv", "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=["id", "uniprot_acc", "role"])
        w.writeheader()
        w.writerows(draft)
    print(f"{len(review)} hits -> candidates_review.tsv; {len(draft)} draft rows -> candidates_draft.csv")


if __name__ == "__main__":
    main()
