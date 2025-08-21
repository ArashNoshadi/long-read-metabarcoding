#!/usr/bin/env python3
"""
numt_filter.py

Six-frame translation NUMT detection and cleaning for COI sequences.

Input: Excel/CSV file with at least columns: ID, RNA Sequence
Outputs:
- <out_prefix>_numt_report.csv  : per-sequence flags and reasons
- <out_prefix>_cleaned.fasta    : retained sequences as FASTA
- prints summary to stdout

Defaults follow manuscript:
- translation table: 5 (invertebrate mitochondrial)
- expected protein length for full COI ~218 aa
- full-record retention threshold: min_protein_len = max(150, ceil(0.60*expected))
- fragment-aware: min_aa_abs = 30, min_aa_pct = 0.80
- GC outlier: mean ± 2*SD (flagged)
"""
import argparse
import math
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path

def translate_six_frames(nuc_seq: str, table: int=5):
    s = str(nuc_seq).replace('-', '').replace('U','T').upper()
    frames = []
    for off in (0,1,2):
        sub = s[off:]
        trim = (len(sub)//3)*3
        prot = str(Seq(sub[:trim]).translate(table=table, to_stop=False))
        frames.append(prot)
    rc = str(Seq(s).reverse_complement())
    for off in (0,1,2):
        sub = rc[off:]
        trim = (len(sub)//3)*3
        prot = str(Seq(sub[:trim]).translate(table=table, to_stop=False))
        frames.append(prot)
    return frames

def has_internal_stop(prot: str):
    return ('*' in prot[:-1]) if len(prot) > 1 else True

def gc_content(seq: str):
    s = str(seq).replace('-', '').upper()
    if len(s) == 0:
        return np.nan
    g = s.count('G')
    c = s.count('C')
    return 100.0 * (g + c) / len(s)

def evaluate_row(seq, expected_protein=218, table=5, min_protein_len_full=None, min_aa_abs=30, min_aa_pct=0.8, gc_mean=None, gc_std=None, gc_sd_threshold=2.0):
    s = str(seq).replace('U','T').replace(' ','').upper()
    length_nt = len(s.replace('-', ''))
    result = {"Length_nt": length_nt, "GC": gc_content(s)}
    if length_nt < 90:
        result.update({"Retained": False, "Reason": "too_short"})
        return result
    # compute min_protein_len_full if not provided
    if min_protein_len_full is None:
        min_protein_len_full = max(150, math.ceil(0.6 * expected_protein))
    frames = translate_six_frames(s, table=table)
    frame_info = []
    any_ok = False
    for prot in frames:
        prot_len = len(prot.replace('*',''))
        internal = has_internal_stop(prot)
        ok = (not internal) and (prot_len >= min_protein_len_full)
        frame_info.append({"protein_len": prot_len, "has_internal_stop": internal, "ok_full": ok})
        if ok:
            any_ok = True
    # fragment-aware test:
    # compute threshold based on sequence length
    fragment_threshold = max(min_aa_abs, math.ceil(min_aa_pct * (length_nt / 3.0)))
    any_ok_fragment = False
    for prot in frames:
        prot_len = len(prot.replace('*',''))
        internal = has_internal_stop(prot)
        ok_frag = (not internal) and (prot_len >= fragment_threshold)
        if ok_frag:
            any_ok_fragment = True
            break

    # decide retention: prefer full-record rule if length ~ expected; otherwise fragment rule
    retained = False
    reason_list = []
    if length_nt >= int(0.9 * expected_protein * 3):  # nearly full-length
        retained = any_ok
        if not any_ok:
            reason_list.append("no_frame_passed_full")
    else:
        retained = any_ok_fragment
        if not any_ok_fragment:
            reason_list.append("no_frame_passed_fragment")

    # GC outlier
    gc_val = result["GC"]
    gc_outlier = False
    if not math.isnan(gc_val) and gc_mean is not None and gc_std is not None and gc_std > 0:
        if (gc_val < (gc_mean - gc_sd_threshold * gc_std)) or (gc_val > (gc_mean + gc_sd_threshold * gc_std)):
            gc_outlier = True
            reason_list.append("gc_outlier")

    result.update({"Retained": bool(retained), "Reason": ";".join(reason_list) if reason_list else "ok", "GC_outlier": gc_outlier, "Frame_Info": frame_info})
    return result

def run_filter(df, seq_col="RNA Sequence", id_col="ID", out_prefix="results/cox1", expected_protein=218, table=5, min_protein_len_full=None, min_aa_abs=30, min_aa_pct=0.8, gc_sd_threshold=2.0):
    df = df.copy()
    df[seq_col] = df[seq_col].astype(str)
    # compute GC mean/std
    df["_gc_temp"] = df[seq_col].apply(lambda s: gc_content(s))
    gc_mean = df["_gc_temp"].mean(skipna=True)
    gc_std = df["_gc_temp"].std(skipna=True)
    results = []
    for idx, row in df.iterrows():
        seq = row[seq_col]
        idv = row.get(id_col, f"row{idx}")
        res = evaluate_row(seq, expected_protein=expected_protein, table=table, min_protein_len_full=min_protein_len_full, min_aa_abs=min_aa_abs, min_aa_pct=min_aa_pct, gc_mean=gc_mean, gc_std=gc_std, gc_sd_threshold=gc_sd_threshold)
        res_row = {"ID": idv, "Length_nt": res["Length_nt"], "GC": res["GC"], "Retained": res["Retained"], "Reason": res["Reason"], "GC_outlier": res["GC_outlier"]}
        results.append(res_row)
    rpt = pd.DataFrame(results)
    # save outputs
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    rpt_csv = f"{out_prefix}_numt_report.csv"
    rpt.to_csv(rpt_csv, index=False)
    # cleaned FASTA
    kept = rpt[rpt["Retained"]]["ID"].values.tolist()
    id_to_seq = {str(row[id_col]): str(row[seq_col]) for _, row in df.iterrows()}
    records = []
    for k in kept:
        seqstr = id_to_seq.get(str(k), "")
        if seqstr:
            records.append(SeqRecord(Seq(seqstr.replace('U','T')), id=str(k), description="retained"))
    fasta_out = f"{out_prefix}_cleaned.fasta"
    SeqIO.write(records, fasta_out, "fasta")
    return rpt_csv, fasta_out

def main():
    parser = argparse.ArgumentParser(description="NUMT filter for COI")
    parser.add_argument("--input", required=True, help="Input csv/xlsx with ID and RNA Sequence")
    parser.add_argument("--seq-col", default="RNA Sequence")
    parser.add_argument("--id-col", default="ID")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--table", type=int, default=5)
    parser.add_argument("--min-protein-full", type=int, default=None)
    parser.add_argument("--min-aa-abs", type=int, default=30)
    parser.add_argument("--min-aa-pct", type=float, default=0.8)
    parser.add_argument("--gc-sd", type=float, default=2.0)
    args = parser.parse_args()
    path = Path(args.input)
    if path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    rpt_csv, fasta_out = run_filter(df, seq_col=args.seq_col, id_col=args.id_col, out_prefix=args.out_prefix, expected_protein=218, table=args.table, min_protein_len_full=args.min_protein_full, min_aa_abs=args.min_aa_abs, min_aa_pct=args.min_aa_pct, gc_sd_threshold=args.gc_sd)
    print("Report:", rpt_csv)
    print("Cleaned FASTA:", fasta_out)

if __name__ == "__main__":
    import argparse
    main()
