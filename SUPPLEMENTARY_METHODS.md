
---

# `SUPPLEMENTARY_METHODS.md`
```markdown
# Supplementary Methods — Exact commands, versions and pipeline parameters

This file documents the EXACT commands, parameters and software versions used in the manuscript analyses. Use it to reproduce the retrieval and curation pipeline.

## Environment & versions
- Python: 3.9+
- Biopython: >= 1.78
- pandas: >= 1.2
- numpy: >= 1.20
- requests: >= 2.25
- tqdm: >= 4.50
- numba: >= 0.52
- psutil: >= 5.8
- openpyxl: >= 3.0
- matplotlib: >= 3.3
- MUSCLE v5 (command-line), optional for alignments
- CD-HIT v4.8.1 (optional for sequence deduplication)

## Entrez / NCBI settings
- Retrieval date: **January 2025**
- Entrez usage: set `ENTREZ_EMAIL` and (optionally) `ENTREZ_API_KEY` as environment variables to avoid throttling.
  - Example:
    ```bash
    export ENTREZ_EMAIL="your_email@example.com"
    export ENTREZ_API_KEY="your_api_key_here"
    ```
- Example ESearch query for COI:
