#!/usr/bin/env python3
"""
gbif_check.py
Simple GBIF name-match helper with local caching.

Usage:
python gbif_check.py --name "Rhabditis sp." --cache cache/gbif.json
"""
import requests
import json
from pathlib import Path
import argparse
import time

GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"

def load_cache(path):
    p = Path(path)
    if p.exists():
        return json.loads(p.read_text(encoding='utf-8'))
    return {}

def save_cache(cache, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(cache, indent=2), encoding='utf-8')

def gbif_match(name, cache_path="cache/gbif_cache.json", pause=0.2):
    cache = load_cache(cache_path)
    if name in cache:
        return cache[name]
    params = {"name": name}
    r = requests.get(GBIF_MATCH_URL, params=params, timeout=15)
    if r.status_code == 200:
        data = r.json()
    else:
        data = {"error": f"status_{r.status_code}"}
    cache[name] = data
    save_cache(cache, cache_path)
    time.sleep(pause)
    return data

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", required=True)
    p.add_argument("--cache", default="cache/gbif_cache.json")
    args = p.parse_args()
    res = gbif_match(args.name, cache_path=args.cache)
    print(json.dumps(res, indent=2))

if __name__ == "__main__":
    main()
