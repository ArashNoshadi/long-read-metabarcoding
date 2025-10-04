# Nematode Long-Read Metabarcoding Pipeline

This repository provides a reproducible implementation of the sequence retrieval, curation, filtering, and analysis pipeline outlined in the manuscript:

> **Noshadi et al. (2025)**. *Paving the way for deeper insights into nematode community composition with long-read metabarcoding*. Soil Biology and Biochemistry. https://doi.org/10.1016/j.soilbio.2025.110001

**Date of primary retrieval:** January 2025

---

## Overview

This README offers a clear, user-friendly summary of the pipeline's structure and components. It is designed to help reviewers, researchers, or new users quickly grasp the purpose of each script, required inputs, and how to execute the pipeline end-to-end—whether on small test datasets for validation or at full scale for comprehensive analyses.

Key highlights:
- The pipeline fetches nematode sequences from NCBI Entrez (via ESearch and EFetch), parses GenBank records into structured tabular formats, applies rigorous quality controls and deduplication, filters potential COI NUMTs (nuclear mitochondrial DNA segments), validates taxonomic names against GBIF, and conducts downstream analyses including sequence alignment, diversity metrics, and pairwise comparisons.
- Provenance is maintained throughout: raw XML files, timestamps, and logs are stored for full auditability and reproducibility.
- Modularity is a core design principle—retrieval, parsing, filtering, and analysis scripts are independent, allowing users to inspect, modify, or run specific stages as needed.

---

## Contents

- **`esearch_download.py`**: Performs ESearch queries on NCBI Entrez to retrieve XML results containing sequence identifiers (e.g., accessions or UIDs) for specified markers; stores raw XML outputs per query/marker for traceability.

- **`extract_ids.py`**: Parses ESearch XML files to extract lists of accessions, GIs, or UIDs; outputs in text or CSV format with added provenance metadata (e.g., source file and timestamp).

- **`efetch_fasta_fetcher.py`**: A batch-oriented EFetch tool that downloads GenBank or FASTA records using ID lists; supports batching, automatic retries, optional compression, and detailed logging (outputs to `data/genbank_raw/<marker>/batch_*.xml`).

- **`fetch_ncbi_to_csv.py`**: Comprehensive parser converting GenBank XML to CSV/Excel; extracts key fields like accession, sequence, definition, organism, host, collection date, location (country/geo_loc_name), and annotated features (e.g., product names and positions); includes provenance columns such as `source_file`, `fetch_timestamp`, and raw XML path.

- **`numt_filter.py`**: Specialized filter for COI sequences to detect and remove potential NUMTs using six-frame translation (invertebrate mitochondrial code, NCBI table 5); identifies near-full-length ORFs, internal stop codons, and GC content outliers; generates a cleaned FASTA file and a detailed per-record report.

- **`sensitivity_summary.py`**: Conducts sensitivity analyses for length thresholds (e.g., 60–100% of expected marker length) and NUMT detection parameters; produces tabular summaries for supplementary materials.

- **`gbif_check.py`**: Validates taxonomic names against the GBIF API with local caching to minimize repeated queries; generates match reports and filters records to retain only those with GBIF phylum assignment as `Nematoda`.

- **`taxa_annotate_and_map.py`**: Standardizes raw organism names to a curated `Main_Organism` field, assigns family-level taxonomy, and annotates with ecological metadata (e.g., c–p values, functional guilds) using predefined mapping tables.

- **`cdhit_dedup.sh`**: Shell wrapper for CD-HIT to perform exact or near-identical sequence deduplication on large FASTA files; includes recommended parameters for efficiency.

- **`merged_analysis_optimized_v2.py`** *(provided separately)*: Central analysis script handling sequence alignment (via MUSCLE/MAFFT interfaces), pairwise similarity/overlap calculations, nucleotide diversity (π), and Tajima’s D; requires aligned FASTA input and outputs Excel workbooks with metrics.

---

## Requirements

- **Python**: 3.9 or higher.

- **Python packages**: `biopython`, `pandas`, `numpy`, `requests`, `tqdm`, `numba`, `psutil`, `openpyxl`, `matplotlib`.
  Install via:
  ```bash:disable-run
  python -m pip install biopython pandas numpy requests tqdm numba psutil openpyxl matplotlib
  ```

- **Optional/recommended tools**:
  - **MUSCLE v5**: For efficient large-scale alignments.
  - **CD-HIT v4.x**: For rapid sequence clustering and deduplication.
  - **OpenCL/PyOpenCL**: For GPU-accelerated computations (code gracefully falls back to CPU if unavailable).

---

## Environment Variables and Configuration

- **`ENTREZ_EMAIL`**: Required for NCBI Entrez compliance; set to your email address.
  Example:
  ```bash
  export ENTREZ_EMAIL="your.email@example.com"
  ```

- **`ENTREZ_API_KEY`**: Optional but highly recommended for high-volume queries to boost rate limits.
  Example:
  ```bash
  export ENTREZ_API_KEY="your_ncbi_api_key"
  ```

- Script configurations are primarily managed via command-line flags (e.g., `--batch-size`, `--sleep` for request delays, `--min-aa`, `--min-aa-pct`, `--expected-length`). Some scripts also support YAML/JSON config files for batch runs.

---

## Quick Examples

**Note**: Ensure `ENTREZ_EMAIL` (and optionally `ENTREZ_API_KEY`) is set before running.

1. **ESearch (Retrieve UIDs for a Marker)**:
   ```bash
   python esearch_download.py \
     --term "txid6231[Organism:exp] AND (COI OR cox1 OR COX1)" \
     --retmax 100000 \
     --out data/raw_esearch/COI_esearch.xml
   ```

2. **Extract IDs from ESearch XML**:
   ```bash
   python extract_ids.py --xml data/raw_esearch/COI_esearch.xml --id-out data/parsed/COI_ids.txt
   ```

3. **Batch EFetch (Download GenBank Records)**:
   ```bash
   python efetch_fasta_fetcher.py \
     --id-file data/parsed/COI_ids.txt \
     --marker COI \
     --batch-size 500 \
     --output-dir data/genbank_raw/COI \
     --compress
   ```

4. **Parse GenBank XML to CSV**:
   ```bash
   python fetch_ncbi_to_csv.py \
     --input-dir data/genbank_raw/COI \
     --out data/parsed/COI_parsed.csv
   ```

5. **NUMT Filtering for COI**:
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

6. **Length Sensitivity Analysis (e.g., for 18S)**:
   ```bash
   python sensitivity_summary.py \
     --input data/parsed/18S_parsed.csv \
     --out SUPPLEMENTARY_TABLES/length_sensitivity_18S.csv \
     --expected-len 1800 \
     --thresholds 0.6 0.7 0.8 0.9 1.0
   ```

7. **GBIF Name Matching**:
   ```bash
   python gbif_check.py \
     --names data/parsed/COI_parsed.csv --name-col organism \
     --cache-dir data/gbif_cache --out SUPPLEMENTARY_TABLES/gbif_match_report.csv
   ```

8. **Deduplicate FASTA with CD-HIT**:
   ```bash
   bash cdhit_dedup.sh input.fasta output_dedup.fasta 1.00
   ```

9. **Core Analysis (After Alignment)**:
   ```bash
   # Align (using provided wrapper or directly with MUSCLE)
   bash alignment_wrapper.sh data/cleaned_fastas/COI_cleaned.fasta data/cleaned_fastas/COI_aligned.fasta
   # Analyze
   python merged_analysis_optimized_v2.py \
     --fasta data/cleaned_fastas/COI_aligned.fasta \
     --taxmap taxonomic_mapping.tsv \
     --out results/genetic_summary_optimized.xlsx
   ```

---

## Pipeline Steps

1. **Retrieval (ESearch → EFetch)**: Queries NCBI Entrez for nematode sequences (taxon ID 6231) filtered by markers; stores raw ESearch XML and batched EFetch outputs in `data/` for complete traceability.

2. **Parsing**: Transforms GenBank XML into tabular formats (CSV/Excel), extracting essential metadata (e.g., accession, sequence, organism, location) and adding provenance details.

3. **Quality Control & Normalization**: Normalizes sequences (e.g., U→T, uppercase, remove non-standard bases), filters ambiguous or low-quality entries (e.g., "uncultured"), deduplicates, and applies marker-specific length thresholds with sensitivity testing.

4. **Taxonomic Validation (GBIF)**: Matches organism names to GBIF database entries (with caching); retains only nematode-confirmed records and generates detailed reports.

5. **NUMT Filtering (COI-Specific)**: Employs six-frame translation to flag NUMT artifacts via ORF integrity, stop codons, and GC anomalies; outputs cleaned FASTA and flagged records.

6. **Alignment & Pairwise Metrics**: Generates alignments (MUSCLE v5 preferred) and computes pairwise identity/overlap on non-gapped positions.

7. **Diversity & Neutrality Analyses**: Calculates nucleotide diversity (π) and Tajima’s D for sufficiently sampled groups, with warnings for low-sample scenarios.

8. **Geographic & Ecological Annotation**: Creates interactive maps (Plotly/HTML) from location metadata; assigns curated taxonomic families and ecological traits (e.g., trophic guilds).

9. **Outputs & Reports**: Produces cleaned FASTAs, CSVs, reports (e.g., NUMT, GBIF), Excel summaries, maps, and logs in structured directories.

---

## Provenance & Reproducibility

- All raw data (XML, CSVs) and logs are preserved in `data/` and `logs/` for re-runs and verification.
- Parsed records include provenance fields like `source_file` and `fetch_timestamp`.
- Include an `environment.yml` or `requirements.txt` with pinned versions; for exact replication, use `conda list --explicit` or `pip freeze > requirements_frozen.txt`.

---

## Tips & Troubleshooting

- **NCBI Rate Limits (HTTP 429)**: Use `ENTREZ_API_KEY`, increase `--sleep`, or reduce `--batch-size`.
- **Geographic Data Gaps**: Supplement with external sources if needed, as many records lack coordinates.
- **NUMT Edge Cases**: Tune parameters and review reports; use sensitivity analyses for optimization.
- **Memory for Large Alignments**: Opt for MUSCLE's `-super5` mode or HPC resources; consider partitioning datasets.
- **GBIF API Issues**: Rely on cached results during outages and retry later.

---

## Suggested Outputs for Reviewers

Include these in manuscript submissions:
- `data/raw_esearch/`: Raw ESearch XML.
- `data/genbank_raw/`: Batched EFetch XML.
- `data/parsed/`: Parsed CSVs/XLSX.
- `data/cleaned_fastas/`: Filtered FASTAs.
- `SUPPLEMENTARY_TABLES/`: Summaries (e.g., filter, length/NUMT sensitivity, GBIF reports, trophic assignments).
- `results/`: Metrics (e.g., genetic summaries, diversity stats) and maps (`*.html`).
- `logs/`: Script execution logs.

For quick validation, provide a `QuickStart.md` with a minimal dataset and one-command execution (avoiding large downloads).

---

## Citation, License & Contact

When using or adapting this pipeline, cite:

> Noshadi et al. (2025). *Paving the way for deeper insights into nematode community composition with long-read metabarcoding*. Soil Biology and Biochemistry. https://doi.org/10.1016/j.soilbio.2025.110001

License: MIT (see `LICENSE` file).

**Contact**: Arash Noshadi — `Nowshadiarash@gmail.com`

---

*This README focuses on clarity and usability for readers and reviewers. For detailed per-script help, run scripts with `--help`. If needed, a `QuickStart.md` can be added for a streamlined test run.*
```
