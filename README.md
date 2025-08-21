# Nematode long-read metabarcoding pipeline

This repository reproduces the retrieval, curation, filtering and analysis pipeline described in the manuscript:

> **He et al.**, *Paving the way for deeper insights into nematode community composition with long-read metabarcoding*

**Date of primary retrieval:** **January 2025**

---

## Overview (reader-facing)

This README is a concise, reader-oriented summary of the pipeline and the role of each script. The goal is to make it straightforward for a reviewer or new user to understand what each component does, what inputs are required, and how to run the pipeline end-to-end on small test data or scale it up for full retrieval.

Key points:

* The pipeline retrieves nematode sequences from NCBI Entrez (ESearch → EFetch), parses GenBank records into tabular format, applies quality control and deduplication, performs COI NUMT filtering, checks taxonomic names against GBIF, and runs downstream analyses (alignment, diversity, pairwise metrics).
* All steps store provenance (raw XML, timestamps, logs) so results can be audited and reproduced.
* The codebase is modular: retrieval, parsing, filtering, and analysis are separate scripts so reviewers can run or inspect exact parts independently.

---

## Contents (scripts & files)

* `esearch_download.py`
  Retrieve ESearch XML results from NCBI Entrez (stores raw XML per query / marker).

* `extract_ids.py`
  Parse ESearch XML files and extract accession/gi/UID lists (text or CSV output). Adds simple provenance (source file + timestamp).

* `efetch_fasta_fetcher.py`
  Batch EFetch utility to download GenBank/FASTA records using lists of UIDs/accessions. Handles batching, retries, optional compression, and logs outputs (`data/genbank_raw/<marker>/batch_*.xml`).

* `fetch_ncbi_to_csv.py`
  Full GenBank → CSV/Excel parser. Extracts accession, definition, full sequence, and annotated features (organism, host, collection\_date, country / geo\_loc\_name, feature products and locations). Adds `source_file`, `fetch_timestamp`, and raw XML path to each parsed row.

* `numt_filter.py`
  COI-specific six-frame NUMT detection and cleaning. Uses invertebrate mitochondrial translation table (NCBI table 5). Flags near-full-length ORFs, internal stops, GC outliers and produces a cleaned FASTA plus per-record report.

* `sensitivity_summary.py`
  Length-threshold sensitivity analysis (e.g., test 60/70/80/90/100% thresholds) and NUMT parameter sweeps; emits tables used in supplementary material.

* `gbif_check.py`
  GBIF name-matching helper with local caching. Queries GBIF species match API and caches responses to avoid repeated API calls. Produces GBIF match reports and a filtered set keeping records where GBIF phylum = `Nematoda`.

* `taxa_annotate_and_map.py`
  Map raw organism strings → curated `Main_Organism`, family assignment, and add ecological/trophic metadata (c–p groups, functional guilds). Uses curated mapping tables for consistency.

* `cdhit_dedup.sh`
  Optional helper script to run CD-HIT for exact or near-identical sequence deduplication (useful on very large FASTA collections). Provided as a shell wrapper with recommended parameters.

* `merged_analysis_optimized_v2.py` *(provided separately)*
  Core analysis: alignment wrapper interface (MUSCLE/MAFFT), pairwise similarity/overlap calculations, nucleotide diversity (π) and Tajima’s D estimates. Expects aligned FASTA as input; produces Excel workbooks with pairwise tables and diversity stats.

---

## Requirements

* **Python:** 3.9+

* **Python packages:** `biopython`, `pandas`, `numpy`, `requests`, `tqdm`, `numba`, `psutil`, `openpyxl`, `matplotlib`
  Install with:

  ```bash
  python -m pip install biopython pandas numpy requests tqdm numba psutil openpyxl matplotlib
  ```

* **Optional / recommended tools:**

  * **MUSCLE v5** (recommended for large alignments)
  * **CD-HIT v4.x** (fast clustering / deduplication)
  * **OpenCL / PyOpenCL** — optional GPU acceleration for computationally heavy routines (the code falls back to CPU if unavailable)

---

## Environment variables & configuration

* `ENTREZ_EMAIL` — **required** by NCBI Entrez; set to the user's contact email.
  Example (bash):

  ```bash
  export ENTREZ_EMAIL="your.email@example.com"
  ```

* `ENTREZ_API_KEY` — optional, strongly recommended for large retrievals to increase rate limits and throughput. Set similarly:

  ```bash
  export ENTREZ_API_KEY="your_ncbi_api_key"
  ```

* Per-script configuration is handled by CLI flags (see usage examples below) or a small YAML/JSON config file (if the script supports it). Typical tunables: `--batch-size`, `--sleep` (seconds between requests), `--min-aa`, `--min-aa-pct`, `--expected-length`.

---

## Quick examples (how to run the main pieces)

> **Note:** The examples assume you have set `ENTREZ_EMAIL` and optionally `ENTREZ_API_KEY`.

1. **ESearch (retrieve UIDs for a marker):**

   ```bash
   python esearch_download.py \
     --term "txid6231[Organism:exp] AND (COI OR cox1 OR COX1)" \
     --retmax 100000 \
     --out data/raw_esearch/COI_esearch.xml
   ```

2. **Extract IDs from ESearch XML:**

   ```bash
   python extract_ids.py --xml data/raw_esearch/COI_esearch.xml --id-out data/parsed/COI_ids.txt
   ```

3. **Batch EFetch (download GenBank XML / FASTA in batches):**

   ```bash
   python efetch_fasta_fetcher.py \
     --id-file data/parsed/COI_ids.txt \
     --marker COI \
     --batch-size 500 \
     --output-dir data/genbank_raw/COI \
     --compress
   ```

4. **Parse GenBank XML → CSV (extract sequences + features):**

   ```bash
   python fetch_ncbi_to_csv.py \
     --input-dir data/genbank_raw/COI \
     --out data/parsed/COI_parsed.csv
   ```

5. **NUMT filtering for COI (six-frame translation):**

   ```bash
   python numt_filter.py \
     --input data/parsed/COI_parsed.csv \
     --acc-col accession \
     --seq-col sequence \
     --fasta-out data/cleaned_fastas/COI_cleaned.fasta \
     --report results/COI_numt_report.csv \
     --min-aa 150 \
     --min-aa-pct 0.8
   ```

6. **Length sensitivity (example for 18S expected length):**

   ```bash
   python sensitivity_summary.py \
     --input data/parsed/18S_parsed.csv \
     --out SUPPLEMENTARY_TABLES/length_sensitivity_18S.csv \
     --expected-len 1800 \
     --thresholds 0.6 0.7 0.8 0.9 1.0
   ```

7. **GBIF name matching & caching:**

   ```bash
   python gbif_check.py \
     --names data/parsed/COI_parsed.csv --name-col organism \
     --cache-dir data/gbif_cache --out SUPPLEMENTARY_TABLES/gbif_match_report.csv
   ```

8. **Optional: deduplicate FASTA with CD-HIT:**

   ```bash
   bash cdhit_dedup.sh input.fasta output_dedup.fasta 1.00
   ```

9. **Core analysis (alignment must be produced first; merged\_analysis expects aligned FASTA):**

   ```bash
   # align with MUSCLE v5 (or via provided wrapper)
   bash alignment_wrapper.sh data/cleaned_fastas/COI_cleaned.fasta data/cleaned_fastas/COI_aligned.fasta
   # then:
   python merged_analysis_optimized_v2.py \
     --fasta data/cleaned_fastas/COI_aligned.fasta \
     --taxmap taxonomic_mapping.tsv \
     --out results/genetic_summary_optimized.xlsx
   ```

---

## Pipeline steps — brief explanation for each stage

1. **Retrieval (ESearch → EFetch)**

   * Query NCBI Entrez with taxon filter (`txid6231[Organism:exp]`) plus gene markers. Save raw ESearch XML and all batch EFetch outputs (raw XML or compressed XML) in `data/` for provenance.

2. **Parsing**

   * Convert GenBank XML to tabular CSV/Excel with per-record fields: accession, sequence, definition, organism, host, collection\_date, country/geo\_loc\_name, lat/long (if present), features and product annotations. Add `source_file` and `fetch_timestamp`.

3. **Quality control & normalization**

   * Nucleotide normalization (U→T, uppercase, strip non-ATCGN), ambiguous label filtering (e.g., `uncultured`, `environmental sample`), deduplication (collapse identical sequences or accessions), and length filtering (marker-specific expected lengths and sensitivity analysis).

4. **Taxonomic cross-checking (GBIF)**

   * Use GBIF species match API to validate organism names and ensure phylum = `Nematoda` where possible. Cache results in `data/gbif_cache/` to avoid repeated API calls.

5. **NUMT filtering (COI)**

   * Six-frame translation with invertebrate mt code (table 5), ORF finding, internal stop checks, GC outlier detection. Produce cleaned FASTA and a CSV report with flags (`pass`, `numt_candidate`, `internal_stops`, `gc_outlier`, `too_short`).

6. **Alignment & pairwise similarity**

   * Produce multiple sequence alignment (MUSCLE v5 recommended for large datasets). Compute pairwise % identity and overlap (non-gap shared positions).

7. **Diversity & neutrality tests**

   * Compute nucleotide diversity (π) and simple Tajima’s D on groups with adequate sample size (script reports stats and warns when sample sizes are too small for reliable inference).

8. **Geographic mapping and ecological annotation**

   * Where collection location metadata exists, export interactive Plotly maps and static PNGs; annotate taxa with curated family and ecological/trophic assignments (trophic c–p groups, functional guilds).

9. **Outputs & reports**

   * Cleaned FASTA files, per-record CSVs, NUMT reports, sensitivity tables, GBIF match tables, Excel workbooks with pairwise and diversity summaries, interactive maps (HTML), and a `logs/` folder with run logs.

---

## Provenance & reproducibility

* Raw XML files, parsed CSVs, and run logs are retained in `data/` and `logs/` to allow exact re-runs and auditing.
* Each parsed record contains `source_file`, `fetch_timestamp`, and `raw_xml_path`.
* Provide `environment.yml` or `requirements.txt` with exact tool versions used. For archival, export `conda list --explicit` or `pip freeze > requirements_frozen.txt`.

---

## Tips & troubleshooting (common issues)

* **NCBI rate limits / HTTP 429**: set `ENTREZ_API_KEY` and increase sleep between requests. Use smaller batch sizes if problems persist.
* **Missing geographic metadata**: many GenBank records lack coordinates. Where possible, contact original submitters or supplement with sample metadata.
* **NUMT edge cases**: adjust `--min-aa` and `--min-aa-pct` and run the sensitivity script. Review `COI_numt_report.csv` for borderline cases.
* **Large alignments memory**: use MUSCLE v5 `-super5` or run alignments on HPC nodes; consider chunked alignment/merging strategies if datasets are enormous.
* **GBIF API failures**: the GBIF helper caches responses; if the API is down, operate from cache and re-run checks when service resumes.

---

## Suggested outputs for reviewers

Provide the following with the manuscript submission:

* `data/raw_esearch/` — raw ESearch XML
* `data/genbank_raw/` — raw EFetch XML batches
* `data/parsed/` — parsed CSV/XLSX records
* `data/cleaned_fastas/` — cleaned FASTA by marker
* `SUPPLEMENTARY_TABLES/` — `filter_summary.xlsx`, `length_sensitivity.xlsx`, `numt_sensitivity.xlsx`, `gbif_match_report.xlsx`, `trophic_cp_assignments.xlsx`
* `results/` — `genetic_summary_optimized.xlsx`, `diversity_stats.xlsx`, `maps/*.html`
* `logs/` — run logs for each script

A `QuickStart` minimal dataset and a single-command run (on small example data) is recommended so reviewers can validate the pipeline quickly without re-downloading large datasets.

---

## Citation, license & contact

If you use or adapt this pipeline, please cite the manuscript:

> He et al. (2025). *Paving the way for deeper insights into nematode community composition with long-read metabarcoding.*

License: include a `LICENSE` file in the repository (MIT is a common, permissive choice). Replace example contact below with a real address in the repository:

**Contact:** Arash Noshadi — `Nowshadiarash@gmail.com`

---

*This README is written for readers and reviewers to understand the structure, purpose, and usage of the pipeline. If you want, I can also expand any of the per-script usage examples into full CLI help text or produce a short `QuickStart.md` with a single reproducible test run (no large downloads).*
