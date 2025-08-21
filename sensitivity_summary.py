#!/usr/bin/env python3
"""
sensitivity_summary.py
Compute counts of sequences meeting length thresholds (60,70,80,90,100%) relative to expected lengths.

Usage:
python sensitivity_summary.py --input-files 18s.xlsx 28s.xlsx cox1.xlsx --expected-lengths '{"18S":1700,"28S":3500,"cox1":700}' --output results/sensitivity_summary.xlsx
"""
import argparse
import pandas as pd
import json
from pathlib import Path

def run(files, expected_lengths, output):
    rows = []
    for f in files:
        marker = Path(f).stem.lower()
        # map typical names
        key = None
        for k in expected_lengths.keys():
            if k.lower() in marker:
                key = k
                break
        if key is None:
            # try filename directly as key
            key = Path(f).stem
        expected = expected_lengths[key]
        df = pd.read_excel(f)
        total = len(df)
        if "Length (RNA)" not in df.columns and "Length (nt)" not in df.columns:
            print(f"Warning: no length column found in {f}; skipping")
            continue
        length_col = "Length (RNA)" if "Length (RNA)" in df.columns else "Length (nt)"
        for p in [60,70,80,90,100]:
            thr = expected * (p/100.0)
            cnt = df[length_col].dropna().ge(thr).sum()
            rows.append({'Marker': key, 'Threshold_%': p, 'Length_threshold': thr, 'Count': int(cnt), 'Percent_of_total': round(100*cnt/total,2)})
    out_df = pd.DataFrame(rows)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_excel(output, index=False)
    print("Saved sensitivity summary:", output)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-files", nargs="+", required=True)
    parser.add_argument("--expected-lengths", required=True, help='JSON e.g. {"18S":1700,"28S":3500,"cox1":700}')
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    expected = json.loads(args.expected_lengths)
    run(args.input_files, expected, args.output)

if __name__ == "__main__":
    main()
