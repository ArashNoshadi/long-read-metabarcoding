#!/usr/bin/env python3
# generate_filter_summary.py
"""
Generate filter_summary.xlsx by aggregating per-step CSVs and parsed files.
Usage:
  python generate_filter_summary.py --parsed parsed/18s_parsed.csv parsed/28s_parsed.csv parsed/cox1_parsed.csv \
    --ambig ambiguous_removed_18S.csv ambiguous_removed_28S.csv ambiguous_removed_cox1.csv \
    --dedup dedup_removed_18S.csv dedup_removed_28S.csv dedup_removed_cox1.csv \
    --numt COI_numt_report.csv --out results/filter_summary.xlsx
"""
import argparse
from pathlib import Path
import pandas as pd

def load_or_empty(p):
    p = Path(p)
    if not p.exists():
        return pd.DataFrame()
    try:
        if p.suffix.lower() in ['.xlsx','.xls']:
            return pd.read_excel(p)
        else:
            return pd.read_csv(p)
    except Exception:
        return pd.DataFrame()

def summarize_counts(df, name):
    if df is None or df.empty:
        return pd.DataFrame([{"Step": name, "Count": 0}])
    return pd.DataFrame([{"Step": name, "Count": len(df)}])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--parsed", nargs="+", help="Parsed files (CSV/XLSX) for markers")
    p.add_argument("--ambig", nargs="*", help="ambiguous_removed CSVs")
    p.add_argument("--dedup", nargs="*", help="dedup_removed CSVs")
    p.add_argument("--numt", help="COI_numt_report.csv")
    p.add_argument("--out", default="results/filter_summary.xlsx")
    args = p.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    writer = pd.ExcelWriter(out_path, engine='openpyxl')

    # 1) Summary of parsed files
    summary_rows = []
    if args.parsed:
        for pf in args.parsed:
            df = load_or_empty(pf)
            summary_rows.append({"File": Path(pf).name, "Parsed_Count": len(df)})
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name="Parsed_counts", index=False)

    # 2) ambiguous removed details
    amb_sheets = []
    if args.ambig:
        for a in args.ambig:
            df = load_or_empty(a)
            sheet = Path(a).stem
            if df.empty:
                pd.DataFrame().to_excel(writer, sheet_name=sheet)
            else:
                # Put a limited preview and full count
                df.head(100).to_excel(writer, sheet_name=sheet, index=False)
                # also write count sheet
                amb_sheets.append({"File": Path(a).name, "Removed": len(df)})
    if amb_sheets:
        pd.DataFrame(amb_sheets).to_excel(writer, sheet_name="ambig_summary", index=False)

    # 3) dedup removed
    dedup_sheets = []
    if args.dedup:
        for d in args.dedup:
            df = load_or_empty(d)
            sheet = Path(d).stem
            if df.empty:
                pd.DataFrame().to_excel(writer, sheet_name=sheet)
            else:
                df.head(100).to_excel(writer, sheet_name=sheet, index=False)
                dedup_sheets.append({"File": Path(d).name, "Removed": len(df)})
    if dedup_sheets:
        pd.DataFrame(dedup_sheets).to_excel(writer, sheet_name="dedup_summary", index=False)

    # 4) NUMT report summary
    if args.numt:
        df = load_or_empty(args.numt)
        if not df.empty:
            # main summary
            total = len(df)
            retained = int(df['Retained'].sum()) if 'Retained' in df.columns else df[df['Reason']=='ok'].shape[0]
            flagged = total - retained
            summary_df = pd.DataFrame([{"Total": total, "Retained": retained, "Flagged_as_NUMT": flagged}])
            summary_df.to_excel(writer, sheet_name="numt_summary", index=False)
            df.head(200).to_excel(writer, sheet_name="numt_report_preview", index=False)
        else:
            pd.DataFrame().to_excel(writer, sheet_name="numt_report_preview")
            pd.DataFrame([{"Total":0}]).to_excel(writer, sheet_name="numt_summary", index=False)

    writer.save()
    print(f"filter_summary written to {out_path}")

if __name__ == "__main__":
    main()
