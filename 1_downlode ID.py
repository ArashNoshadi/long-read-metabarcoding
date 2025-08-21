import requests
from tqdm import tqdm

def download_with_progress(url: str, output_path: str) -> None:
    """
    Download the content at `url` and write it to `output_path` as text,
    while displaying a progress bar in the terminal.
    """
    # Single GET (streaming) — no HEAD
    with requests.get(url, stream=True) as resp:
        resp.raise_for_status()
        total_size = int(resp.headers.get('Content-Length', 0))

        # Open the output file in text mode
        with open(output_path, "w", encoding="utf-8") as f, \
             tqdm(
                 total=total_size or None,   # None => indeterminate bar if size unavailable
                 unit='B',
                 unit_scale=True,
                 unit_divisor=1024,
                 desc="Downloading"
             ) as bar:
            for chunk in resp.iter_content(chunk_size=8192):
                if not chunk:
                    continue
                text = chunk.decode('utf-8', errors='replace')
                f.write(text)
                bar.update(len(chunk))

    print(f"\nDownload complete: {output_path}")

if __name__ == "__main__":
    URL = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        "esearch.fcgi?db=nuccore&term=txid6231%5BOrganism%3Aexp%5D+cox"
        "&retmode=xml&retmax=5635000"
    )
    OUTPUT_FILE = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\data\nema-cox.txt"
    download_with_progress(URL, OUTPUT_FILE)
