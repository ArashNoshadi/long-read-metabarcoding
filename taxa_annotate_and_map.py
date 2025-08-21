#!/usr/bin/env python3
"""
taxa_annotate_and_map.py
Annotate parsed GenBank CSV/Excel with Main_Organism, Family mapping, cp and f-h values.

Inputs:
- parsed CSV/Excel (columns should include 'Organism' or similar)
- an optional family lookup (CSV) with columns: Family, cp, f-h

Outputs:
- annotated Excel with added columns: Main_Organism, Family, Species, cp, f-h
"""
import argparse
import pandas as pd
import re
from pathlib import Path
from gbif_check import gbif_match  # local import (ensure gbif_check.py is in same folder)

MAIN_NAME_RE = re.compile(r'^\s*([A-Za-z][A-Za-z\-]+)')

BAD_VALUES = {"uncultured","environmental","metagenome","unknown","nematode","sample","sp."}

def extract_main_name(organism):
    if pd.isna(organism): return ""
    m = MAIN_NAME_RE.match(str(organism))
    return m.group(1).lower() if m else str(organism).split()[0].lower()

def extract_species(organism):
    if pd.isna(organism): return ""
    parts = str(organism).split()
    if len(parts) >= 2:
        return " ".join(parts[1:])
    return ""

def detect_family_by_lookup(main_name, family_lookup):
    # family_lookup: dict genus->family or family->genera mapping depending on file
    if not main_name:
        return "unknown"
    # direct mapping genus->family
    fam = family_lookup.get(main_name.lower())
    if fam:
        return fam
    return "unknown"

def load_family_lookup(path):
    if not Path(path).exists():
        return {}
    df = pd.read_csv(path)
    # expect columns: Genus, Family
    lookup = {}
    if "Genus" in df.columns and "Family" in df.columns:
        for _, r in df.iterrows():
            lookup[str(r["Genus"]).strip().lower()] = str(r["Family"]).strip().lower()
    return lookup

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--family-lookup", default="data/family_lookup.csv")
    p.add_argument("--gbif-cache", default="cache/gbif_cache.json")
    args = p.parse_args()

    path = Path(args.input)
    if path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)

    family_lookup = load_family_lookup(args.family_lookup)
    df["Organism"] = df["Organism"].astype(str)
    df["Main_Organism"] = df["Organism"].apply(extract_main_name)
    df["Species"] = df["Organism"].apply(extract_species)
    df["Family"] = df["Main_Organism"].apply(lambda x: detect_family_by_lookup(x, family_lookup))

    # GBIF check for names that remain unknown or suspect
    for idx, r in df.iterrows():
        if df.at[idx, "Family"] == "unknown" and df.at[idx, "Main_Organism"]:
            nm = df.at[idx, "Main_Organism"]
            res = gbif_match(nm, cache_path=args.gbif_cache)
            if res and isinstance(res, dict) and res.get("phylum", "").lower() == "nematoda":
                if res.get("family"):
                    df.at[idx, "Family"] = res.get("family").lower()

    # load cp/f-h mapping if exist
    try:
        lookup2 = pd.read_csv("data/family_ecology_lookup.csv")
        lookup_map = {str(r["Family"]).strip().lower(): {"cp": r.get("cp",""), "f-h": r.get("f-h","")} for _, r in lookup2.iterrows()}
    except Exception:
        lookup_map = {}

    df["cp"], df["f-h"] = zip(*df["Family"].apply(lambda x: (lookup_map.get(x, {}).get("cp","Unknown"), lookup_map.get(x, {}).get("f-h","Unknown"))))

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(args.output, index=False)
    print("Saved annotated file:", args.output)

if __name__ == "__main__":
    main()
