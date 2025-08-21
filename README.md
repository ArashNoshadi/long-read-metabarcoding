Nematode Long-Read Metabarcoding Pipeline
This repository contains the complete pipeline for retrieval, curation, filtering, and analysis of nematode metabarcoding data, as described in the manuscript:

He et al., Paving the way for deeper insights into nematode community composition with long-read metabarcoding

Date of primary retrieval: January 2025
Overview
This pipeline processes long-read metabarcoding data for nematode communities, targeting three marker genes (18S, 28S, COI) from NCBI Entrez, performing quality control, taxonomic cross-checking, NUMT filtering for COI, and downstream analyses including alignment, nucleotide diversity, and geographic mapping. All steps are fully reproducible, with detailed logs, provenance tracking, and supplementary outputs for transparency.
Repository Structure
├── README.md
├── SUPPLEMENTARY_METHODS.md
├── environment.yml
├── requirements.txt
├── data/
│   ├── raw_esearch/                # Raw ESearch XML per marker
│   ├── genbank_raw/                # Raw EFetch returns (per accession)
│   ├── parsed/                     # Parsed CSV / Excel outputs
│   ├── cleaned_fastas/             # Cleaned FASTA per marker
│   ├── data_examples/              # Small example dataset for testing
├── scripts/
│   ├── fetch_records_from_ids.py   # Batch EFetch for FASTA/GenBank
│   ├── parse_esearch_xml.py        # Parse ESearch XML and extract IDs
│   ├── numt_filter.py              # Six-frame NUMT detection for COI
│   ├── sensitivity_summary.py      # Length-threshold sensitivity analysis
│   ├── merged_analysis_optimized_v2.py # Core analysis (alignment, metrics)
│   ├── alignment_wrapper.sh        # Wrappers for MUSCLE / CD-HIT
├── logs/                           # Run logs
├── results/                        # Analysis outputs (e.g., genetic_summary_optimized.xlsx)
├── SUPPLEMENTARY_TABLES/           # Supplementary tables (filter summaries, sensitivity)

Requirements

Python: 3.9+
Packages: biopython, pandas, numpy, requests, tqdm, numba, psutil, openpyxl, matplotlib, plotly
Optional: MUSCLE v5, CD-HIT v4.8.1, OpenCL/PyOpenCL for GPU acceleration

Install dependencies:
python -m pip install biopython pandas numpy requests tqdm numba psutil openpyxl matplotlib plotly

Conda environment:
conda env create -f environment.yml
conda activate nema_pipeline

Pipeline Steps
1. Data Retrieval

Time window: January 2025 (exact timestamps in SUPPLEMENTARY_METHODS.md)
Target: Nematoda (txid6231) sequences for 18S, 28S, and COI
Method: Entrez ESearch/EFetch via Biopython
Query template:from Bio import Entrez
Entrez.email = os.environ.get("ENTREZ_EMAIL")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")  # Optional
term_template = "txid6231[Organism:exp] AND ({gene_query})"
handle = Entrez.esearch(db="nuccore", term=term_template.format(gene_query="18S OR SSU OR small subunit"),
                        retmax=100000, retmode="xml")
esearch_xml = handle.read()


Batch EFetch: Batch size = 500 (configurable, logged in SUPPLEMENTARY_METHODS.md)
Outputs:
Raw ESearch XML: data/raw_esearch/<marker>_esearch.xml
Raw EFetch XML: data/genbank_raw/<marker>_efetch_batch_*.xml
Parsed CSV/Excel: data/parsed/<marker>_parsed.csv, Supplementary_Material_1.xlsx



2. Parsing & Provenance

Script: scripts/parse_esearch_xml.py
Parsed fields: Accession, definition, sequence, features (organism, host, collection_date, country/geo_loc_name, product), feature locations
Outputs: data/parsed/<marker>_parsed.csv, data/parsed/<marker>_records.xlsx
Provenance: Each row includes source_file, fetch_timestamp, raw_xml_path

3. Quality Control

Logs: Saved in logs/
Summaries: SUPPLEMENTARY_TABLES/filter_summary.xlsx

3.1 Nucleotide Normalization

Replace U→T, convert to uppercase, remove non-ATCGN characters (ambiguous N kept if ≤ 5%)
Outputs: data/cleaned_fastas/<marker>_cleaned.fasta

3.2 Ambiguous-Label Filtering

Filter terms: scripts/bad_value_list.txt (e.g., "uncultured", "environmental sample")
Outputs: SUPPLEMENTARY_TABLES/ambiguous_removed_<marker>.csv

3.3 Deduplication

Collapse duplicate accession IDs (keep first, map in dedup_accession_map.csv)
Collapse identical sequences (pandas/MD5 or CD-HIT with -c 1.0)
Outputs: SUPPLEMENTARY_TABLES/filter_summary.xlsx

3.4 Length Filtering

Threshold: ≥ 80% of expected length per marker (defined in SUPPLEMENTARY_METHODS.md)
Sensitivity analysis: 60%, 70%, 80%, 90%, 100% thresholds
Outputs: SUPPLEMENTARY_TABLES/length_sensitivity.xlsx

4. Taxonomic Cross-Checking (GBIF)

Query GBIF API: https://api.gbif.org/v1/species/match?name=<urlencoded_name>
Cache: data/gbif_cache/
Keep records with GBIF phylum = "Nematoda"
Outputs: SUPPLEMENTARY_TABLES/ambiguous_taxonomic_matches.csv, gbif_match_report.xlsx

5. NUMT Detection (COI-specific)

Script: scripts/numt_filter.py
Steps:
Normalize: U→T, strip non-ATCG, discard len < 90 nt
Six-frame translation (invertebrate mitochondrial code, NCBI table 5)
Rules:
Near-full-length COI: ≥ 150 aa, no internal stops
Fragments: ORF ≥ 30 aa and ≥ 80% of expected length
GC outliers: Exclude if outside mean ± 2*SD and fails ORF




Outputs: COI_numt_report.csv, data/cleaned_fastas/COI_cleaned.fasta
Sensitivity: Tested min_aa_abs = {120,150,180,200}, min_aa_pct = {0.6,0.8,1.0}
Results: SUPPLEMENTARY_TABLES/numt_sensitivity.xlsx



6. Alignment & Pairwise Similarity

Tools: MUSCLE v5 (-super5 for large datasets), scripts/alignment_wrapper.sh
Grouping: By Family/Genus/Species using taxonomic_mapping.tsv
Metrics:
Similarity: % identical positions (Nei & Li, 1979)
Overlap: % shared non-gap positions
Changes: % mismatches or gaps


Implementation: NumPy/Numba, optional OpenCL in merged_analysis_optimized_v2.py
Outputs: results/genetic_summary_optimized.xlsx

7. Nucleotide Diversity & Tajima’s D

Target: Species with ≥ 2 sequences (preferably ≥ 200)
Metrics:
Nucleotide diversity (π): Vectorized NumPy
Tajima’s D: Compared to Watterson’s θ (Tajima, 1989)


Outputs: results/diversity_stats.xlsx

8. Geographic Mapping

Source: country, geo_loc_name, or explicit lat/long
Visualization: Plotly Express scatter_geo (size ~ log(sample_count))
Outputs: results/maps/<marker>_geo_map.html, static PNGs

9. Taxonomic & Ecological Assignments

Family: Based on Andrássy (2005–2009), Abebe et al. (2006), Ahmad & Jairajpuri (2010), Pena-Santiago (2021)
Trophic/c–p groups: Ferris (2010)
Outputs: SUPPLEMENTARY_TABLES/trophic_cp_assignments.xlsx

10. Reproducibility

Environment: environment.yml or requirements.txtname: nema_pipeline
dependencies:
  - python=3.9
  - biopython=1.79
  - pandas=1.3.5
  - numpy=1.21
  - numba=0.53
  - matplotlib=3.5
  - plotly=5.8
  - muscle=5.1
  - cd-hit=4.8.1


Provenance: Raw XMLs, logs, commit history, SUPPLEMENTARY_METHODS.md
QuickStart: data_examples/ and QuickStart.md for end-to-end testing

Outputs for Reviewers

Tables: SUPPLEMENTARY_TABLES/{filter_summary, numt_sensitivity, length_sensitivity, trophic_cp_assignments, gbif_match_report}.xlsx
Data: data/parsed/*, data/cleaned_fastas/*, data/raw_esearch/*, data/genbank_raw/*
Logs: logs/
Results: results/{genetic_summary_optimized, diversity_stats}.xlsx, results/maps/*.html

Example Script: fetch_records_from_ids.py
Below is a complete script for batch EFetching GenBank records.

import os
import time
from Bio import Entrez
import xml.etree.ElementTree as ET
from tqdm import tqdm
import logging

Setup logging
logging.basicConfig(    filename="logs/fetch_records.log",    level=logging.INFO,    format="%(asctime)s - %(levelname)s - %(message)s")
Entrez configuration
Entrez.email = os.environ.get("ENTREZ_EMAIL")Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
def batch_fetch_records(id_list, db="nuccore", batch_size=500, output_dir="data/genbank_raw"):    """Fetch GenBank records in batches and save as XML."""    os.makedirs(output_dir, exist_ok=True)    total_ids = len(id_list)
for start in tqdm(range(0, total_ids, batch_size), desc="Fetching batches"):
    end = min(start + batch_size, total_ids)
    batch_ids = id_list[start:end]
    
    try:
        handle = Entrez.efetch(
            db=db,
            id=",".join(batch_ids),
            retmode="xml",
            rettype="gb"
        )
        xml_data = handle.read()
        handle.close()
        
        # Save batch XML
        batch_file = f"{output_dir}/batch_{start//batch_size + 1}.xml"
        with open(batch_file, "wb") as f:
            f.write(xml_data)
        
        logging.info(f"Saved batch {start//batch_size + 1} to {batch_file}")
        
        # Respect NCBI rate limits
        time.sleep(0.34)  # ~3 requests per second with API key
        
    except Exception as e:
        logging.error(f"Error fetching batch {start//batch_size + 1}: {str(e)}")
        continue

def read_id_list(file_path):    """Read list of accession IDs from file."""    with open(file_path, "r") as f:        return [line.strip() for line in f if line.strip()]
def main():    # Example: Read IDs from parsed ESearch results    id_file = "data/parsed/id_list.txt"  # Assume this exists from parse_esearch_xml.py    id_list = read_id_list(id_file)
# Fetch records for each marker
markers = ["18S", "28S", "COI"]
for marker in markers:
    logging.info(f"Starting fetch for {marker}")
    output_dir = f"data/genbank_raw/{marker}"
    batch_fetch_records(id_list, output_dir=output_dir)
    logging.info(f"Completed fetch for {marker}")

if name == "main":    main()
