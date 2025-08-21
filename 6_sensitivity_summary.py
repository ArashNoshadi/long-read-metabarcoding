import argparse
import os
import pandas as pd

DEFAULT_PATHS = {
    '18S': r'G:\Paper\nema-Nanopore-Sequencing\by all gene\combined_results-6\18S_updated.xlsx',
    '28S': r'G:\Paper\nema-Nanopore-Sequencing\by all gene\combined_results-6\28S_updated.xlsx',
    'cox1': r'G:\Paper\nema-Nanopore-Sequencing\by all gene\combined_results-6\cox1_updated.xlsx',
}
EXPECTED = {'18S': 1700, '28S': 3500, 'cox1': 700}
THRESHOLDS = [60, 70, 80, 90, 100]

def analyze_file(path, expected_len):
    out_rows = []
    if not os.path.exists(path):
        print(f"Warning: input not found: {path}")
        return out_rows
    try:
        df = pd.read_excel(path, engine='openpyxl')
    except Exception as e:
        print(f"Error reading {path}: {e}")
        return out_rows

    # Accept potential column names for length
    possible_cols = ["Length (RNA)", "Length", "Length_RNA", "length", "length_rna"]
    length_col = None
    for c in possible_cols:
        if c in df.columns:
            length_col = c
            break

    if length_col is None:
        print(f"No length column found in {path}. Skipping length counts for this file.")
        return out_rows

    lengths = pd.to_numeric(df[length_col], errors='coerce').dropna()
    total = int(len(lengths))
    for p in THRESHOLDS:
        thr = expected_len * (p / 100.0)
        cnt = int((lengths >= thr).sum())
        out_rows.append({
            "Marker": os.path.basename(path),
            "Marker_short": None,  # fill later
            "Threshold_pct": p,
            "Length_threshold": thr,
            "Count_ge_threshold": cnt,
            "Total_sequences_with_length": total,
            "Percent_ge_threshold": (cnt / total * 100) if total > 0 else 0.0
        })
    return out_rows

def main(args):
    # build marker->path mapping
    inputs = {}
    if args.inputs:
        # if user passed 3 paths, map in order; otherwise expect 'key=path' pairs
        if len(args.inputs) == 3:
            inputs = {'18S': args.inputs[0], '28S': args.inputs[1], 'cox1': args.inputs[2]}
        else:
            # parse key=path
            for item in args.inputs:
                if '=' in item:
                    k,v = item.split('=',1)
                    inputs[k] = v
    else:
        inputs = DEFAULT_PATHS.copy()

    rows = []
    for marker, path in inputs.items():
        marker_key = marker
        expected = EXPECTED.get(marker_key, None)
        if expected is None:
            print(f"No expected length for marker '{marker_key}', skipping.")
            continue
        res = analyze_file(path, expected)
        for r in res:
            r["Marker_short"] = marker_key
            rows.append(r)

    if not rows:
        print("No results to write.")
        return

    out_df = pd.DataFrame(rows)
    out_excel = args.out
    
    # --- خطوط اضافه شده برای حل مشکل ---
    # نام پوشه را از مسیر کامل فایل خروجی استخراج می‌کند
    output_dir = os.path.dirname(out_excel)
    # پوشه را در صورتی که وجود نداشته باشد، ایجاد می‌کند
    os.makedirs(output_dir, exist_ok=True)
    # ------------------------------------

    out_csv = os.path.splitext(out_excel)[0] + ".csv"
    out_df.to_excel(out_excel, index=False)
    out_df.to_csv(out_csv, index=False)
    print(f"Saved sensitivity summary to: {out_excel} and {out_csv}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Sensitivity summary for length thresholds.")
    p.add_argument('--inputs', nargs='+', help="Either three paths (18S 28S cox1) or key=path pairs")
    p.add_argument('--out', default=r'G:\Paper\nema-Nanopore-Sequencing\by all gene\sensitivity_summary\sensitivity_summary.xlsx',
                   help="Output Excel file (default under by all gene)")
    args = p.parse_args()
    main(args)