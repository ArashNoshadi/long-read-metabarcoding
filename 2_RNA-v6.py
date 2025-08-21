#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fetch_ncbi_to_csv.py

ورودی: فایل ID (هر خط یک ID) یا یک XML شامل تگ <Id>...</Id>
خروجی: CSV حاوی رکوردهای کامل NCBI (ستون‌ها قابل توسعه‌اند)
"""

import os
import re
import io
import time
import ssl
import certifi
import socket
import argparse
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import pandas as pd
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
from Bio import Entrez

# ---------- تنظیمات پیش‌فرض ----------
os.environ['SSL_CERT_FILE'] = certifi.where()
socket.setdefaulttimeout(10)

# مقداردهی Entrez (ایمیل را حتماً تغییر دهید یا از متغیر محیطی استفاده کنید)
Entrez.email = os.environ.get("NCBI_EMAIL", "nowshadiarash@gmail.com")
Entrez.tool  = "nema_itss_fetcher"
Entrez.api_key = os.environ.get("NCBI_API_KEY", "81b8f961fd2e2bab9a24e2bbe90fb6d17609")  # اگر ندارید رشته خالی بذارید

# ---------- لاگینگ ----------
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

# ---------- نقشه‌ی markerها (قابل توسعه) ----------
marker_map = [
    # ── 5S ریبوزومی ─────────────────────────────────────────────────────────────
    (r"\b5[\s\-]?s[\s\-]?rRNA\b",                          "5S"),

    # ── External Transcribed Spacer (ETS) ──────────────────────────────────────
    (r"\bexternal[\s\-]?transcribed[\s\-]?spacer\b",       "ETS"),
    (r"\bets\b",                                           "ETS"),

    # ── Internal Transcribed Spacer (ITS) ──────────────────────────────────────
    (r"\bITS[\s\-]?region\b",                              "ITS"),
    (r"\bits[\s\-]?1\b",                                   "ITS1"),
    (r"\bits[\s\-]?2\b",                                   "ITS2"),
    (r"\binternal[\s\-]?transcribed[\s\-]?spacer[\s\-]?1\b","ITS1"),
    (r"\binternal[\s\-]?transcribed[\s\-]?spacer[\s\-]?2\b","ITS2"),

    # ── 18S (SSU) ────────────────────────────────────────────────────────────────
    (r"\b18[\s\-]?s\b",                                     "18S"),
    (r"\bssu[\s\-]?rrna\b",                                 "18S"),
    (r"\brRNA[\s\-]?small[\s\-]?subunit\b",                 "18S"),
    (r"\bsmall[\s\-]?subunit[\s\-]?ribosomal\b",            "18S"),
    (r"\b18[\s\-]?s[-\s]?rDNA\b",                           "18S"),

    # ── 28S (LSU) ───────────────────────────────────────────────────────────────
    (r"\b28[\s\-]?s\b",                                     "28S"),
    (r"\blsu[\s\-]?rrna\b",                                 "28S"),
    (r"\brRNA[\s\-]?large[\s\-]?subunit\b",                 "28S"),
    (r"\blarge[\s\-]?subunit[\s\-]?ribosomal\b",            "28S"),
    (r"\b28[\s\-]?s[-\s]?rDNA\b",                           "28S"),

    # ── 12S (mito SSU) ──────────────────────────────────────────────────────────
    (r"\b12[\s\-]?s\b",                                     "12S"),
    (r"\brrn[sS]\b",                                        "12S"),

    # ── 16S (mito LSU) ──────────────────────────────────────────────────────────
    (r"\b16[\s\-]?s\b",                                     "16S"),
    (r"\brrn[lL]\b",                                        "16S"),

    # ── 5.8S ─────────────────────────────────────────────────────────────────────
    (r"\b5[.\s\-]?8[\s\-]?s\b",                             "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?rRNA\b",                  "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?subunit\b",               "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?ribosomal\b",             "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?ribosomal[\s\-]?RNA\b",   "5.8S"),

    # ── 23S ─────────────────────────────────────────────────────────────────────
    (r"\b23[\s\-]?s\b",                                     "23S"),
    (r"\b23[\s\-]?s[\s\-]?rRNA\b",                          "23S"),

    # ── ژن‌های میتوکندریایی: COX1, COX2, COX3, CYTB ────────────────────────────
    #   COX1
    (r"\bcox[\s\-]?1\b",                                    "COX1"),
    (r"\bco[iI]\b",                                         "COX1"),
    (r"\bcytochrome[\s\-]?c[\s\-]?oxidase[\s\-]?subunit[\s]?1\b","COX1"),
    (r"\bcytochrome[\s\-]?oxidase[\s\-]?subunit\s?I\b",     "COX1"),
    #   COX2
    (r"\bcox[\s\-]?2\b",                                    "COX2"),
    (r"\bcytochrome[\s\-]?c[\s\-]?oxidase[\s\-]?subunit[\s]?2\b","COX2"),
    (r"\bcytochrome[\s\-]?oxidase[\s\-]?subunit\s?II\b",    "COX2"),
    #   COX3
    (r"\bcox[\s\-]?3\b",                                    "COX3"),
    (r"\bcytochrome[\s\-]?c[\s\-]?oxidase[\s\-]?subunit[\s]?3\b","COX3"),
    (r"\bcytochrome[\s\-]?oxidase[\s\-]?subunit\s?III\b",   "COX3"),
    #   CYTB
    (r"\bcytb\b",                                           "CYTB"),
    (r"\bcytochrome[\s\-]?b\b",                             "CYTB"),

    # ── NADH Dehydrogenase (ND1–ND6) ────────────────────────────────────────────
    (r"\bnd([1-6])\b",                                      lambda m: f"ND{m.group(1)}"),
    (r"\bnad[hH]?[\s\-]?dehydrogenase[\s\-]?subunit[\s]?([1-6])\b",
                                                           lambda m: f"ND{m.group(1)}"),

    # ── ژن‌های ATP ───────────────────────────────────────────────────────────────
    (r"\batp[fF][\s\-]?atp[hH]\b",                          "atpF-atpH"),
    (r"\batp[\s\-]?ase[\s\-]?subunit\s?6\b",                "ATP6"),
    (r"\batp6\b",                                           "ATP6"),
    (r"\batp[\s\-]?synthetase\s?subunit\s?alpha\b",         "ATPα"),

    # ── گیاهی (rbcL, matK) ───────────────────────────────────────────────────────
    (r"\brbc[\s\-]?l\b",                                    "rbcL"),
    (r"\bmat[\s\-]?k\b",                                    "matK"),

    # ── ترو-psb و psb-k/psb-i ────────────────────────────────────────────────────
    (r"\btrn[hH][\s\-]?psb[aA]\b",                          "trnH-psbA"),
    (r"\bpsb[kK][\s\-]?psb[iI]\b",                          "psbK-psbI"),

    # ── متغیرهای منطقه‌ای 18S/28S ───────────────────────────────────────────────
    (r"\bssu[-\s]?v4\b",                                    "18S-V4"),
    (r"\bssu[-\s]?v9\b",                                    "18S-V9"),
    (r"\blsu[-\s]?d2\b",                                    "28S-D2"),
    (r"\b28[\s\-]?s[-\s]?d2\b",                             "28S-D2"),

    # ── اضافه: تشخیص 'mitochondrion' و 'mitochondrial' ─────────────────────────
    (r"\bmitochondrion\b",                                  "mitochondrion"),
    (r"\bmitochondrial\b",                                  "mitochondrion"),
]


compiled_markers = [(re.compile(p, re.IGNORECASE), v) for p, v in marker_map]

# ---------- ستون‌های خروجی ----------
FIELDNAMES = [
    "ID","DEFINITION","marker","RNA Sequence","Length (RNA)",
    "Host","Geo Loc Name","Collection Date","Organism","Mol Type","Product"
]

# ---------- توابع کمکی ----------
def ensure_parent(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

def load_ids(path: Path):
    text = path.read_text(encoding='utf-8', errors='ignore')
    # اگر XML باشد: تگ <Id>12345</Id>
    ids = re.findall(r"<Id>(\d+)</Id>", text)
    if ids:
        logger.info("Found %d IDs from XML-style file.", len(ids))
        return ids
    # در غیر این صورت هر خط را بررسی کن
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    # اگر خطوط شامل چیزهایی غیر عددی باشند، سعی کن فقط اعداد بیرون بکشی
    if all(re.fullmatch(r"\d+", s) for s in lines):
        logger.info("Found %d IDs from plain-line file.", len(lines))
        return lines
    # ممکن است فایل شامل ترکیبی از معیارها باشد؛ جستجوی طولانی‌تر برای ارقام
    ids_from_any = re.findall(r"\b(\d{5,})\b", text)  # IDs طولانی‌تر از 5 رقم
    logger.info("Extracted %d numeric IDs (fallback).", len(ids_from_any))
    return ids_from_any

def check_internet(host="8.8.8.8", port=53, timeout=3):
    try:
        socket.setdefaulttimeout(timeout)
        s = socket.socket()
        s.connect((host, port))
        s.close()
        return True
    except OSError:
        return False

def parse_record(rec):
    desc = rec.description or ""
    found = []
    for regex, canon in compiled_markers:
        for m in regex.finditer(desc):
            val = canon(m) if callable(canon) else canon
            found.append(val)
    # حذف تکراری‌ها و حفظ ترتیب
    seen = set()
    markers = [x for x in found if not (x in seen or seen.add(x))]
    marker_str = "|".join(markers) if markers else "N/A"

    try:
        seq_str = str(rec.seq)
        length = len(rec.seq)
    except (UndefinedSequenceError, TypeError):
        seq_str = "N/A"
        length = 0

    entry = dict.fromkeys(FIELDNAMES, "N/A")
    entry.update({
        "ID": rec.id,
        "DEFINITION": desc,
        "marker": marker_str,
        "RNA Sequence": seq_str,
        "Length (RNA)": length
    })

    for feat in rec.features:
        q = feat.qualifiers
        if feat.type == "source":
            entry.update({
                "Host":           q.get("host", ["N/A"])[0],
                "Geo Loc Name":   q.get("geo_loc_name", ["N/A"])[0],
                "Collection Date":q.get("collection_date", ["N/A"])[0],
                "Organism":       q.get("organism", ["N/A"])[0],
                "Mol Type":       q.get("mol_type", ["N/A"])[0]
            })
        elif feat.type == "CDS":
            entry["Product"] = q.get("product", ["N/A"])[0]

    return entry

# ---------- تابع fetch برای یک bulk از IDها ----------
def fetch_bulk(id_list, efetch_retries=3, delay=0.1, email=None, api_key=None):
    """
    ورودی: id_list => لیست ID به صورت رشته
    خروجی: لیست دیکشنری‌های parse_record
    """
    time.sleep(delay)
    ids_str = ",".join(id_list)
    try:
        post = Entrez.epost(db="nucleotide", id=ids_str)
        data = Entrez.read(post); post.close()
        we, qk = data["WebEnv"], data["QueryKey"]
    except Exception as e:
        logger.error("EPost FAILED: %s", e)
        return []

    recs = []
    for attempt in range(1, efetch_retries + 1):
        try:
            params = {"db":"nucleotide","WebEnv":we,"query_key":qk,
                      "rettype":"gb","retmode":"text"}
            # در صورتی که api_key موجود باشد آن را اضافه کن
            if api_key:
                params["api_key"] = api_key
            resp = requests.get(Entrez.base + "efetch.fcgi", params=params, timeout=30, verify=True)
            resp.raise_for_status()
            recs = list(SeqIO.parse(io.StringIO(resp.text), "gb"))
            break
        except Exception as e:
            logger.warning("EFetch attempt %d failed: %s", attempt, e)
            time.sleep(delay * (2 ** (attempt-1)))
    if not recs:
        logger.error("All EFetch retries failed for chunk starting with %s", id_list[:1])
        return []
    return [parse_record(r) for r in recs]

# ---------- تابع برای نوشتن تدریجی CSV ----------
def append_df_to_csv(df: pd.DataFrame, path: Path):
    # در صورتی که فایل وجود ندارد هدر بنویسد، وگرنه append بدون هدر
    header = not path.exists()
    df.to_csv(path, mode='a', index=False, header=header, encoding='utf-8')

# ---------- تابع اصلی پردازش ----------
def process_ids_to_csv(input_path, output_csv,
                       batch_size=200, bulk_fetch_size=50, max_workers=3,
                       efetch_retries=3, delay=0.1, check_interval=3):
    input_path = Path(input_path)
    output_csv = Path(output_csv)
    ensure_parent(output_csv)

    ids = load_ids(input_path)
    if not ids:
        logger.error("No IDs found in %s", input_path)
        return

    # تقسیم به batch های بزرگ‌تر (برای کنترل جریان کاری)
    batches = [ids[i:i+batch_size] for i in range(0, len(ids), batch_size)]
    logger.info("Total IDs: %d, Batches: %d", len(ids), len(batches))

    total_written = 0

    for bi, batch in enumerate(batches, 1):
        logger.info("Processing batch %d/%d (size=%d)", bi, len(batches), len(batch))
        # اطمینان از وصل بودن اینترنت
        while not check_internet():
            logger.warning("No internet. Sleeping %d seconds...", check_interval)
            time.sleep(check_interval)

        # شکستن batch به قطعات کوچکتر برای epost+efetch
        chunks = [batch[i:i+bulk_fetch_size] for i in range(0, len(batch), bulk_fetch_size)]
        results = []
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [pool.submit(fetch_bulk, ch, efetch_retries, delay, Entrez.email, Entrez.api_key) for ch in chunks]
            for fut in as_completed(futures):
                try:
                    res = fut.result() or []
                    results.extend(res)
                except Exception as e:
                    logger.error("Worker error: %s", e)

        if results:
            df = pd.DataFrame.from_records(results, columns=FIELDNAMES)
            append_df_to_csv(df, output_csv)
            total_written += len(df)
            logger.info("Appended %d rows to %s (total written: %d)", len(df), output_csv, total_written)
        else:
            logger.warning("No records fetched for this batch")

        # کوتاه مهلت بین batch ها (قابل تنظیم)
        time.sleep(0.5)

    logger.info("Done. Total rows written: %d", total_written)

# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(description="Fetch NCBI records for a list of IDs and save to CSV.")
    p.add_argument("-i", "--input", required=True, help="Path to input file (one ID per line or XML with <Id>)")
    p.add_argument("-o", "--output", required=True, help="Path to output CSV file")
    p.add_argument("--batch-size", type=int, default=200, help="Batch size for outer loop")
    p.add_argument("--bulk-size", type=int, default=50, help="Chunk size for epost/efetch")
    p.add_argument("--workers", type=int, default=3, help="Max worker threads")
    p.add_argument("--retries", type=int, default=3, help="EFetch retries per chunk")
    p.add_argument("--delay", type=float, default=0.1, help="Base delay between requests")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    logger.info("Starting fetch: input=%s output=%s", args.input, args.output)
    process_ids_to_csv(
        input_path=args.input,
        output_csv=args.output,
        batch_size=args.batch_size,
        bulk_fetch_size=args.bulk_size,
        max_workers=args.workers,
        efetch_retries=args.retries,
        delay=args.delay
    )
