#!/usr/bin/env python3
"""For each >id|acc|role entry in a FASTA, pull the AlphaFold DB model and report
mean pLDDT (CA atoms) and the fraction of residues with pLDDT >= 70."""
import json, os, sys, urllib.request


def get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def main(fasta, out_tsv, model_dir):
    os.makedirs(model_dir, exist_ok=True)
    rows = []
    for line in open(fasta):
        if not line.startswith(">"):
            continue
        cid, acc, role = line[1:].strip().split("|")[:3]
        try:
            meta = json.loads(get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"))[0]
            pdb = get(meta["pdbUrl"]).decode()
        except Exception:
            rows.append((cid, acc, role, "NA", "NA", "no_model"))
            continue
        with open(f"{model_dir}/{cid}_{acc}.pdb", "w") as fh:
            fh.write(pdb)
        pl = [float(l[60:66]) for l in pdb.splitlines()
              if l.startswith("ATOM") and l[12:16].strip() == "CA"]
        mean = sum(pl) / len(pl)
        frac = sum(p >= 70 for p in pl) / len(pl)
        rows.append((cid, acc, role, f"{mean:.1f}", f"{frac:.2f}", "ok"))
    with open(out_tsv, "w") as o:
        o.write("id\tacc\trole\tmean_plddt\tfrac_plddt_ge70\tstatus\n")
        for r in rows:
            o.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
