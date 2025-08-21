# Nematode long-read metabarcoding pipeline

This repository reproduces the retrieval, curation, filtering and analysis pipeline described in the manuscript:
> He et al., *Paving the way for deeper insights into nematode community composition with long-read metabarcoding*

Date of primary retrieval: **January 2025**

## Contents
- `esearch_download.py` – retrieve ESearch XML results from NCBI Entrez (raw XML)
- `extract_ids.py` – parse ESearch XML and extract IDs
- `efetch_fasta_fetcher.py` – batch EFetch to download FASTA/GenBank
- `fetch_ncbi_to_csv.py` – full GenBank -> CSV parser (including features)
- `numt_filter.py` – six-frame NUMT detection and cleaning for COI
- `sensitivity_summary.py` – length-threshold sensitivity analysis
- `gbif_check.py` – GBIF name-matching helper (caching)
- `taxa_annotate_and_map.py` – map organism→Main_Organism/Family and add cp / f-h
- `cdhit_dedup.sh` – optional CD-HIT dedupe script for FASTA files
- `merged_analysis_optimized_v2.py` – core analysis (alignment, metrics, reports) *(provided separately)*

## Requirements
- Python 3.9+
- Packages: `biopython`, `pandas`, `numpy`, `requests`, `tqdm`, `numba`, `psutil`, `openpyxl`, `matplotlib`
- Optional: MUSCLE v5, CD-HIT, OpenCL/PyOpenCL if GPU acceleration desired.

Install Python dependencies:
```bash
python -m pip install biopython pandas numpy requests tqdm numba psutil openpyxl matplotlib
