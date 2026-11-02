#!/usr/bin/env python3
"""Fetch UniProt FASTA for each row of a CSV (id,uniprot_acc,role).
Output header: >id|acc|role. Rows with empty or TO_CURATE accessions are skipped."""
import csv, sys, time, urllib.request, urllib.error


def fetch(acc, tries=3):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 * (i + 1))
        except urllib.error.URLError:
            time.sleep(2 * (i + 1))
    return None


def main(csv_in, fasta_out):
    n = 0
    with open(csv_in) as f, open(fasta_out, "w") as out:
        for row in csv.DictReader(f):
            acc = row["uniprot_acc"].strip()
            if not acc or acc.startswith("TO_CURATE"):
                print(f"skip {row['id']}: no accession", file=sys.stderr)
                continue
            txt = fetch(acc)
            if not txt:
                print(f"WARN: could not fetch {acc}", file=sys.stderr)
                continue
            seq = "".join(txt.splitlines()[1:])
            out.write(f">{row['id']}|{acc}|{row['role']}\n{seq}\n")
            n += 1
    if n == 0:
        sys.exit("no sequences fetched - check the CSV")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
