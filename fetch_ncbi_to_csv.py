#!/usr/bin/env python3
"""
fetch_ncbi_to_csv.py
Given an ID file (one ID per line), fetch GenBank records in chunks using Entrez.epost+efetch,
parse records (Biopython SeqIO) and write a CSV/Excel summarizing fields.

Output columns include:
ID, DEFINITION, marker (inferred), RNA Sequence, Length (RNA), Host, Geo Loc Name, Collection Date, Organism, Mol Type, Product
"""
import argparse
from pathlib import Path
import time
import io
import logging
import requests
import pandas as pd
from Bio import Entrez, SeqIO

# Adjust Entrez settings via environment or modify here
Entrez.email = "nowshadiarash@gmail.com"  # replace or use env var
Entrez.api_key = ""  # optional: set via env var if available

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

FIELDNAMES = [
    "ID","DEFINITION","marker","RNA Sequence","Length (RNA)",
    "Host","Geo Loc Name","Collection Date","Organism","Mol Type","Product"
]

def parse_record(rec):
    desc = rec.description or ""
    seq_str = str(rec.seq) if rec.seq is not None else ""
    length = len(seq_str)
    entry = dict.fromkeys(FIELDNAMES, "N/A")
    entry.update({
        "ID": rec.id,
        "DEFINITION": desc,
        "marker": "N/A",
        "RNA Sequence": seq_str,
        "Length (RNA)": length
    })
    for feat in rec.features:
        if feat.type == "source":
            q = feat.qualifiers
            entry.update({
                "Host":           q.get("host", ["N/A"])[0],
                "Geo Loc Name":   q.get("geo_loc_name", ["N/A"])[0],
                "Collection Date":q.get("collection_date", ["N/A"])[0],
                "Organism":       q.get("organism", ["N/A"])[0],
                "Mol Type":       q.get("mol_type", ["N/A"])[0]
            })
        elif feat.type == "CDS":
            q = feat.qualifiers
            entry["Product"] = q.get("product", ["N/A"])[0]
    return entry

def epost_efetch_chunks(id_list, bulk=50, retries=3, delay=0.1):
    results = []
    for i in range(0, len(id_list), bulk):
        chunk = id_list[i:i+bulk]
        ids_str = ",".join(chunk)
        # epost
        try:
            post = Entrez.epost(db="nucleotide", id=ids_str)
            data = Entrez.read(post); post.close()
            we, qk = data["WebEnv"], data["QueryKey"]
        except Exception as e:
            logger.error("EPost failed for chunk: %s", e)
            continue
        # efetch
        for attempt in range(1, retries+1):
            try:
                fetch = Entrez.efetch(db="nucleotide", WebEnv=we, query_key=qk, rettype="gb", retmode="text")
                text = fetch.read(); fetch.close()
                recs = list(SeqIO.parse(io.StringIO(text), "gb"))
                results.extend(recs)
                break
            except Exception as e:
                logger.warning("EFetch attempt %d failed: %s", attempt, e)
                time.sleep(delay*(2**(attempt-1)))
    return results

def ids_from_file(path):
    p = Path(path)
    lines = [l.strip() for l in p.read_text(encoding='utf-8').splitlines() if l.strip()]
    return lines

def main():
    p = argparse.ArgumentParser(description="Fetch GenBank records and write CSV")
    p.add_argument("-i","--input", required=True)
    p.add_argument("-o","--output", required=True)
    p.add_argument("--bulk-size", type=int, default=50)
    args = p.parse_args()

    ids = ids_from_file(args.input)
    if not ids:
        logger.error("No IDs found in input.")
        return
    logger.info("Fetching %d IDs...", len(ids))
    recs = epost_efetch_chunks(ids, bulk=args.bulk_size)
    parsed = [parse_record(r) for r in recs]
    df = pd.DataFrame.from_records(parsed, columns=FIELDNAMES)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    logger.info("Saved CSV: %s (%d records)", out, len(df))

if __name__ == "__main__":
    main()
