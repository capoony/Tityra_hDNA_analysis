#!/bin/bash

###############################################################################
# Aggregate Mitochondrial Gene Metadata
# Description: Parse the headers of the concatenated genes fasta, which list
#              every gene/sample used in the analysis (COI, CYTB, ND2, CO2),
#              extract the species and accession ID per gene, look up the
#              matching individual metadata in data/Mito/{gene}_{genus}.tsv,
#              and write a single consolidated metadata table.
###############################################################################

set -euo pipefail

# Set working directory
WD=/media/inter/mkapun/projects/Tityra

# Directories
CONCAT_DIR="${WD}/results/concatenated_genes_assembly"
MITO_DIR="${WD}/data/Mito"
OUTPUT_DIR="${CONCAT_DIR}"

# Input concatenated fasta
CONCAT_FASTA="${CONCAT_DIR}/concatenated_all_genes.fasta"

# Log file
LOG="${WD}/logs/aggregate_mito_metadata.log"
mkdir -p "${WD}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "==================================================================="
echo "Aggregating mitochondrial gene metadata"
echo "Started: $(date)"
echo "==================================================================="

mkdir -p "${OUTPUT_DIR}"

OUT_TABLE="${OUTPUT_DIR}/mito_gene_metadata.csv"

python3 << PYEOF
import os
import re
import csv

concat_fasta = "${CONCAT_FASTA}"
mito_dir = "${MITO_DIR}"
out_table = "${OUT_TABLE}"

# Metadata columns to carry over from the per-individual TSVs (Seq intentionally excluded)
META_COLS = [
    "GenBank_ID", "TaxName", "GenBankTitle", "Journal",
    "SubmissionTitle", "Consortium", "Lat_Lon", "CollCountry",
    "Tissue", "Identifier", "Collector"
]

GENES = ["COI", "CYTB", "ND2", "CO2"]

def load_metadata_tsv(path):
    """Load a metadata TSV keyed by Accession. Returns dict or None if missing."""
    if not os.path.isfile(path):
        return None
    table = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = (row.get("Accession") or "").strip()
            if acc:
                table[acc] = row
    return table

rows_out = []

if not os.path.isfile(concat_fasta):
    print(f"Error: concatenated fasta not found: {concat_fasta}")
    raise SystemExit(1)

with open(concat_fasta) as f:
    for line in f:
        line = line.strip()
        if not line.startswith(">"):
            continue
        # Header format:
        # >Species genes=N/4 COI:Species|Accession CYTB:Species|Accession ...
        body = line[1:].strip()
        # First whitespace-delimited token is the overall species
        first_tok = body.split()[0] if body.split() else ""

        # Parse each GENE:VALUE token (VALUE = "missing" or "Species|Accession")
        for token in body.split():
            m = re.match(r"^(" + "|".join(GENES) + r"):(.+)$", token)
            if not m:
                continue
            gene = m.group(1)
            value = m.group(2)
            if value == "missing":
                # Gene absent for this sample
                continue

            # value is "Species|Accession"
            if "|" in value:
                species_gene, accession = value.split("|", 1)
            else:
                species_gene, accession = value, ""
            accession = accession.strip()

            genus = species_gene.split("_")[0] if "_" in species_gene else species_gene

            # Look up metadata
            meta = None
            tsv_path = os.path.join(mito_dir, f"{gene}_{genus}.tsv")
            table = load_metadata_tsv(tsv_path)
            if table is not None and accession and accession != "NA":
                meta = table.get(accession)

            record = {"Gene": gene, "Species": species_gene, "Accession": accession}
            for col in META_COLS:
                record[col] = (meta.get(col, "") if meta else "") or "NA"
            rows_out.append(record)

# Write consolidated table
fieldnames = ["Gene", "Species", "Accession"] + META_COLS
with open(out_table, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="$")
    writer.writeheader()
    for r in rows_out:
        writer.writerow(r)

print(f"Wrote {len(rows_out)} sample metadata rows to: {out_table}")
PYEOF

echo ""
echo "==================================================================="
echo "Metadata aggregation completed: $(date)"
echo "==================================================================="
echo ""
echo "Table saved to: ${OUT_TABLE}"