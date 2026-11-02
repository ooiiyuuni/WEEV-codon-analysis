#!/usr/bin/env python3
# ============================================================
# 01_extract_metadata.py
# Description: Extract metadata for 52 WEEV sequences (Table S1)
# Input:  data/WEEV/52.gb
# Output: results/tables/metadata.csv
# Columns: Accession, Version, Definition, Organism, Host, Country, Isolate
# Dependencies: biopython (tested with v1.81+)
# Run: python scripts/01_extract_metadata.py
# ============================================================

from pathlib import Path
from Bio import SeqIO
import csv

ROOT     = Path(__file__).resolve().parent.parent
GB_FILE  = ROOT / "data" / "WEEV" / "52.gb"
OUT_FILE = ROOT / "results" / "tables" / "metadata.csv"

# Ensure output directory exists
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

rows = []
for record in SeqIO.parse(str(GB_FILE), "genbank"):
    host, country, isolate, strain = "", "", "", ""
    for feature in record.features:
        if feature.type == "source":
            q = feature.qualifiers
            host    = q.get("host",    [""])[0]
            country = q.get("country", q.get("geo_loc_name", [""]))[0]
            isolate = q.get("isolate", [""])[0]
            strain  = q.get("strain",  [""])[0]
            break

    acc         = record.id.split(".")[0]   # remove version
    version_str = record.id                 # keep version

    rows.append({
        "Accession":  acc,
        "Version":    version_str,
        "Definition": record.description,
        "Organism":   record.annotations.get("organism", ""),
        "Host":       host,
        "Country":    country,
        "Isolate":    isolate if isolate else strain,
    })

with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"\nExtracted {len(rows)} sequences.")
print(f"Results saved to {OUT_FILE}\n")

for r in rows[:5]:
    print(r)