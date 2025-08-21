
---

# `esearch_download.py`
```python
#!/usr/bin/env python3
"""
esearch_download.py
Download ESearch XML from NCBI Entrez using a given query term.
Outputs: XML file with ESearch results.
"""
import argparse
import requests
from pathlib import Path
from time import sleep
from tqdm import tqdm

NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

def download_esearch_xml(term, db="nuccore", retmax=1000000, out_path="esearch.xml", retmode="xml", delay=0.34):
    params = {"db": db, "term": term, "retmax": retmax, "retmode": retmode}
    print(f"ESearch: term={term[:120]}... retmax={retmax}")
    with requests.get(NCBI_ESEARCH_URL, params=params, stream=True, timeout=60) as r:
        r.raise_for_status()
        out = Path(out_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "wb") as fh:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    fh.write(chunk)
    print(f"Saved ESearch XML -> {out_path}")
    # polite delay
    sleep(delay)

def main():
    p = argparse.ArgumentParser(description="Run ESearch and save XML.")
    p.add_argument("--term", required=True, help="Entrez ESearch term (wrap in quotes)")
    p.add_argument("--db", default="nuccore")
    p.add_argument("--retmax", type=int, default=5000000)
    p.add_argument("--out", default="raw_esearch.xml")
    args = p.parse_args()
    download_esearch_xml(args.term, db=args.db, retmax=args.retmax, out_path=args.out)

if __name__ == "__main__":
    main()
