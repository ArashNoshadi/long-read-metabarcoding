#!/usr/bin/env python3
"""
extract_ids.py
Parse an ESearch XML file and extract <Id> entries into a plain text file (one ID per line).
"""
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

def extract_ids(xml_path, out_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    ids = [elem.text for elem in root.findall(".//Id")]
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        for id_ in ids:
            fh.write(id_.strip() + "\n")
    print(f"Extracted {len(ids)} IDs to {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Extract IDs from ESearch XML")
    parser.add_argument("--input", required=True, help="ESearch XML file")
    parser.add_argument("--output", required=True, help="Output ID file (one ID per line)")
    args = parser.parse_args()
    extract_ids(args.input, args.output)

if __name__ == "__main__":
    main()
