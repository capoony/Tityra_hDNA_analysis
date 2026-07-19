#!/usr/bin/env python3
"""
Create per-locus FASTA files combining all species from NEXUS alignments
with the Tityra_leucura consensus sequence, then build supermatrix + partitions.
"""

import os
import sys
import re
from collections import OrderedDict

# Paths
WD = "/media/inter/mkapun/projects/Tityra"
NEXUS_DIR = os.path.join(WD, "data/openwings_tityrinae/mafft-nexus-clean-54taxa-95per")
CONSENSUS_FILE = os.path.join(WD, "results/tityra_uce/mapping/Tityra_leucura_uce_consensus.fasta")
RESULTS = os.path.join(WD, "results/tityra_uce/phylogeny")
GENE_SEQS_DIR = os.path.join(RESULTS, "gene_sequences")
ALIGNMENTS_DIR = os.path.join(RESULTS, "alignments")
SUPERMATRIX_DIR = os.path.join(RESULTS, "supermatrix")
LOGS = os.path.join(WD, "logs")

# Reference genome data (from the original NEXUS files)
TITYRA_PATTERN = re.compile(r'^(Tityra_|Tit_)')


def parse_nexus(nexus_path):
    """Parse a NEXUS file and return a dict of {taxon_name: sequence_string}."""
    with open(nexus_path, 'r') as f:
        content = f.read()

    matrix_start = content.find('matrix')
    if matrix_start == -1:
        return {}

    matrix_text = content[matrix_start + len('matrix'):]

    end_idx = matrix_text.find(';')
    if end_idx == -1:
        return {}

    matrix_text = matrix_text[:end_idx].strip()

    sequences = {}
    current_name = None
    current_seq = ""

    for line in matrix_text.split('\n'):
        line = line.strip()
        if not line:
            if current_name and current_seq:
                sequences[current_name] = current_seq
                current_name = None
                current_seq = ""
            continue

        if current_name is not None:
            match = re.match(r'^(\S+)\s+', line)
            if match and match.group(1) and not all(c in 'ACGTacgt?- ' for c in match.group(1)):
                if current_name and current_seq:
                    sequences[current_name] = current_seq
                parts = line.split(None, 1)
                current_name = parts[0]
                current_seq = parts[1] if len(parts) > 1 else ""
            else:
                current_seq += line
        else:
            parts = line.split(None, 1)
            if len(parts) >= 1:
                current_name = parts[0]
                current_seq = parts[1] if len(parts) > 1 else ""

    if current_name and current_seq:
        sequences[current_name] = current_seq

    return sequences


def clean_sequence(seq):
    """Remove spaces, convert to uppercase."""
    return seq.replace(' ', '').replace('\t', '').replace('\n', '').replace('\r', '').upper()


def count_valid_bases(seq):
    """Count non-missing, non-gap characters."""
    return sum(1 for c in seq if c not in '?- \t\n\r')


def find_best_tityra(tityra_seqs):
    """Among Tityra sequences, find the one with most valid bases. Return (name, cleaned_seq)."""
    best_name = None
    best_count = -1
    best_seq = None

    for name, seq in tityra_seqs.items():
        clean = clean_sequence(seq)
        valid = count_valid_bases(clean)
        if valid > best_count:
            best_count = valid
            best_name = name
            best_seq = clean

    return best_name, best_seq


def main():
    # Create directories
    for d in [GENE_SEQS_DIR, ALIGNMENTS_DIR, SUPERMATRIX_DIR, LOGS]:
        os.makedirs(d, exist_ok=True)

    # Step 1: Load consensus sequences
    print("Loading Tityra_leucura consensus sequences...", file=sys.stderr)
    consensus = {}
    current = None
    with open(CONSENSUS_FILE) as f:
        for line in f:
            if line.startswith('>'):
                current = line[1:].strip()
            else:
                consensus[current] = consensus.get(current, '') + line.strip()

    print(f"  Loaded {len(consensus)} consensus sequences", file=sys.stderr)

    # Step 2: Create a mapping from consensus header -> Tityra_leucura
    # The consensus headers are like "uce-1005_Tityra_inquisitor_UFAC_422"
    # The Tityra naming behind the uce-XXXX_ may vary per locus
    # We need to extract the locus prefix to find the best match
    cons_mapping = {}  # locus -> Tityra_leucura sequence
    for header, seq in consensus.items():
        # Extract uce-XXXX prefix
        m = re.match(r'^(uce-\d+)', header)
        if m:
            locus = m.group(1)
            cons_mapping[locus] = seq
        else:
            print(f"  WARNING: cannot parse consensus header: {header}", file=sys.stderr)

    print(f"  Mapped {len(cons_mapping)} consensus sequences by locus", file=sys.stderr)

    # Step 3: Process each NEXUS file
    nexus_files = sorted([f for f in os.listdir(NEXUS_DIR)
                         if f.endswith('.nexus') and f.startswith('uce-')])

    print(f"\nProcessing {len(nexus_files)} NEXUS files...", file=sys.stderr)

    all_species = set()
    loci_used = []

    for filename in nexus_files:
        nexus_path = os.path.join(NEXUS_DIR, filename)
        locus = filename.replace('.nexus', '')

        # Parse NEXUS
        sequences = parse_nexus(nexus_path)
        if not sequences:
            continue

        # Get Tityra_leucura consensus for this locus
        leucura_seq = cons_mapping.get(locus)
        if leucura_seq is None:
            print(f"  SKIP {locus}: no consensus sequence available", file=sys.stderr)
            continue

        # Build FASTA for this locus
        # All non-Tityra samples from the reference + Tityra_leucura
        out_path = os.path.join(GENE_SEQS_DIR, f"{locus}.fasta")
        species_in_locus = []

        with open(out_path, 'w') as out:
            for name, seq in sequences.items():
                clean = clean_sequence(seq)
                out.write(f">{name}\n{clean}\n")
                species_in_locus.append(name)

            # Add Tityra_leucura consensus
            out.write(f">Tityra_leucura\n{leucura_seq}\n")
            species_in_locus.append("Tityra_leucura")

            all_species.update(species_in_locus)
            loci_used.append(locus)

    species_sorted = sorted(all_species)
    print(f"\nFound {len(species_sorted) - 1} reference species + Tityra_leucura", file=sys.stderr)
    print(f"Produced {len(loci_used)} per-locus FASTA files", file=sys.stderr)

    # Write species list
    species_list_path = os.path.join(RESULTS, "species_list.txt")
    with open(species_list_path, 'w') as f:
        for sp in species_sorted:
            f.write(f"{sp}\n")

    print(f"\nSpecies list: {species_list_path}", file=sys.stderr)
    print(f"Gene sequences: {GENE_SEQS_DIR}/", file=sys.stderr)

    # Step 4: Generate alignment and partition commands
    # Write a list of loci to process
    loci_list_path = os.path.join(RESULTS, "loci_list.txt")
    with open(loci_list_path, 'w') as f:
        for locus in loci_used:
            f.write(f"{locus}\n")

    # Write partition info (to be used after alignment)
    # Each locus becomes one partition
    partition_path = os.path.join(SUPERMATRIX_DIR, "partitions_template.txt")
    with open(partition_path, 'w') as f:
        for locus in loci_used:
            f.write(f"{locus}\n")

    print(f"\nDone. Now run MAFFT alignments and supermatrix construction (see UCEs.sh)", file=sys.stderr)
    print(f"  Loci: {len(loci_used)}", file=sys.stderr)
    print(f"  Species: {len(species_sorted)}", file=sys.stderr)


if __name__ == "__main__":
    main()