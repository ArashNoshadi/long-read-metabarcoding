#!/usr/bin/env python3
# plot_geo_map.py
"""
Plot geographic distribution of sequences by country using Plotly.
Input: CSV/XLSX with either:
 - aggregated rows: columns ['country','count']
 - raw rows: includes 'Geo Loc Name' or 'country' for each record; script will aggregate counts by country.
Outputs:
 - interactive HTML (plotly), PNG (static)
Usage:
  python plot_geo_map.py --input data/cox1_parsed.csv --country-col "Geo Loc Name" --out-prefix results/cox1_geo
"""
import argparse
from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.io as pio
import pycountry

def normalize_country_name(name):
    if pd.isna(name): return None
    s = str(name).strip()
    # try pycountry
    try:
        return pycountry.countries.lookup(s).name
    except Exception:
        # fallback
        return s

def run(input_file, country_col, out_prefix):
    df = pd.read_excel(input_file) if str(input_file).lower().endswith(('.xlsx','.xls')) else pd.read_csv(input_file)
    if country_col not in df.columns:
        # attempt common columns
        for c in ['country','Country','geo_loc_name','Geo Loc Name']:
            if c in df.columns:
                country_col = c
                break
    countries = df[country_col].astype(str).fillna('').str.strip()
    countries = countries.replace({'': None})
    df['country_norm'] = countries.apply(lambda x: normalize_country_name(x) if x and x.lower()!='nan' else None)
    agg = df.groupby('country_norm').size().reset_index(name='count').dropna(subset=['country_norm'])
    if agg.empty:
        print("No country data found.")
        return
    # map country names to ISO alpha-3 codes required by plotly
    def to_iso3(name):
        try:
            return pycountry.countries.lookup(name).alpha_3
        except Exception:
            return None
    agg['iso3'] = agg['country_norm'].apply(to_iso3)
    agg_valid = agg.dropna(subset=['iso3'])
    if agg_valid.empty:
        print("No valid ISO3 country codes; check country names.")
    fig = px.scatter_geo(agg_valid, locations="iso3", color="count", size="count",
                         hover_name="country_norm", projection="natural earth",
                         title="Sequence counts by country")
    html_out = f"{out_prefix}_geo.html"
    png_out = f"{out_prefix}_geo.png"
    Path(html_out).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(html_out)
    # write static PNG (requires kaleido)
    try:
        pio.write_image(fig, png_out, scale=2)
    except Exception as e:
        print("Warning: could not write PNG (kaleido missing). Error:", e)
    print(f"Saved interactive map: {html_out} and static PNG: {png_out}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--country-col", default="Geo Loc Name")
    p.add_argument("--out-prefix", default="results/geo")
    args = p.parse_args()
    run(args.input, args.country_col, args.out_prefix)

if __name__ == "__main__":
    main()
