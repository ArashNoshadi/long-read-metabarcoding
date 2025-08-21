#!/usr/bin/env python3
"""
efetch_fasta_fetcher.py
Batch-download FASTA and GenBank records via NCBI EFetch for a list of IDs.

Outputs:
- combined FASTA file
- combined GenBank file (text)
"""
import argparse
import requests
from pathlib import Path
from time import sleep
from tqdm import tqdm

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def fetch_batch(ids, rettype, retmode="text"):
    params = {"db": "nuccore", "id": ",".join(ids), "rettype": rettype, "retmode": retmode}
    r = requests.get(EUTILS_URL, params=params, timeout=60)
    r.raise_for_status()
    return r.text

def fetch_ids_to_files(id_file, out_fasta, out_gb, batch_size=500, delay=0.34):
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(out_gb).parent.mkdir(parents=True, exist_ok=True)

    with open(id_file, "r", encoding="utf-8") as fh:
        ids = [line.strip() for line in fh if line.strip()]

    with open(out_fasta, "w", encoding="utf-8") as fasta_fh, open(out_gb, "w", encoding="utf-8") as gb_fh:
        for i in tqdm(range(0, len(ids), batch_size), desc="Fetching batches"):
            batch = ids[i:i+batch_size]
            # FASTA
            try:
                fasta_text = fetch_batch(batch, rettype="fasta", retmode="text")
                fasta_fh.write(fasta_text)
            except Exception as e:
                print(f"FASTA batch error at {i}: {e}")
            sleep(delay)
            # GenBank
            try:
                gb_text = fetch_batch(batch, rettype="gb", retmode="text")
                gb_fh.write(gb_text)
            except Exception as e:
                print(f"GenBank batch error at {i}: {e}")
            sleep(delay)
    print(f"Saved FASTA -> {out_fasta} and GenBank -> {out_gb}")

def main():
    parser = argparse.ArgumentParser(description="Batch efetch FASTA/GenBank for IDs")
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-gb", required=True)
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--delay", type=float, default=0.34)
    args = parser.parse_args()
    fetch_ids_to_files(args.id_file, args.out_fasta, args.out_gb, batch_size=args.batch_size, delay=args.delay)

if __name__ == "__main__":
    main()
