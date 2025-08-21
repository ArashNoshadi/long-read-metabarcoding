#!/usr/bin/env python3
# generate_sensitivity_plot.py
"""
Generate sensitivity plot (60-100%) based on length columns present in marker files.
Usage:
  python generate_sensitivity_plot.py --files data/18s.xlsx data/28s.xlsx data/cox1.xlsx \
    --expected '{"18S":1700,"28S":3500,"cox1":700}' --out results/sensitivity_plot.png --out-xlsx results/sensitivity_summary.xlsx
"""
import argparse, json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def get_marker_key(fname, expected_keys):
    name = Path(fname).stem.lower()
    for k in expected_keys:
        if k.lower() in name:
            return k
    # fallback: return first key
    return list(expected_keys)[0]

def run(files, expected, out_png, out_xlsx):
    thresholds = [60,70,80,90,100]
    rows = []
    for f in files:
        df = pd.read_excel(f)
        key = get_marker_key(f, expected.keys())
        expected_len = expected[key]
        if "Length (RNA)" in df.columns:
            length_col = "Length (RNA)"
        elif "Length (nt)" in df.columns:
            length_col = "Length (nt)"
        else:
            # try to compute from sequence column
            if "RNA Sequence" in df.columns:
                df["Length_tmp"] = df["RNA Sequence"].astype(str).str.replace('-','').str.len()
                length_col = "Length_tmp"
            else:
                continue
        total = len(df)
        for p in thresholds:
            thr = expected_len * (p/100.0)
            cnt = int(df[length_col].dropna().ge(thr).sum())
            rows.append({"Marker": key, "Threshold_%": p, "Count": cnt, "Percent": round(100*cnt/total if total>0 else 0,2)})
    out_df = pd.DataFrame(rows)
    out_df.to_excel(out_xlsx, index=False)
    # plot
    plt.figure(figsize=(8,5))
    for key, grp in out_df.groupby("Marker"):
        plt.plot(grp['Threshold_%'], grp['Count'], marker='o', label=key)
    plt.xlabel("Length threshold (%) of expected length")
    plt.ylabel("Number of sequences retained")
    plt.title("Length-threshold sensitivity (60-100%)")
    plt.grid(True)
    plt.legend()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Saved plot {out_png} and table {out_xlsx}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--files", nargs="+", required=True)
    p.add_argument("--expected", required=True, help='JSON string e.g. {"18S":1700,"28S":3500,"cox1":700}')
    p.add_argument("--out", default="results/sensitivity_plot.png")
    p.add_argument("--out-xlsx", default="results/sensitivity_summary.xlsx")
    args = p.parse_args()
    expected = json.loads(args.expected)
    run(args.files, expected, args.out, args.out_xlsx)

if __name__ == "__main__":
    main()
