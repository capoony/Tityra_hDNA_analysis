#!/usr/bin/env python3
"""
Extract from each UCE nexus file the Tityra sample with the longest sequence
(without missing '?' or gap '-' characters) and concatenate into a single multifasta file.
"""

import os
import re
import sys

# Configuration
NEXUS_DIR = "/media/inter/mkapun/projects/Tityra/data/openwings_tityrinae/mafft-nexus-clean-54taxa-95per"
OUTPUT_FASTA = "/media/inter/mkapun/projects/Tityra/data/Tityra_concat_UCEs.fasta"

# Tityra samples matching pattern
TITYRA_PATTERN = re.compile(r'^(Tityra_|Tit_)')


def parse_nexus(nexus_path):
    """Parse a NEXUS file and return a dict of {taxon_name: sequence_string}."""
    with open(nexus_path, 'r') as f:
        content = f.read()

    # Find the matrix block
    matrix_start = content.find('matrix')
    if matrix_start == -1:
        print(f"  [WARN] No 'matrix' keyword found in {nexus_path}", file=sys.stderr)
        return {}

    # Content after 'matrix' keyword - skip the 'matrix' line itself
    matrix_text = content[matrix_start + len('matrix'):]

    # Find end of matrix (semicolon)
    end_idx = matrix_text.find(';')
    if end_idx == -1:
        print(f"  [WARN] No ';' semicolon found in {nexus_path}", file=sys.stderr)
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

        # If we're accumulating sequence for a current taxon
        if current_name is not None:
            # Check if this line starts a new taxon (starts with a letter and contains whitespace before ACGT)
            # A new taxon name is followed by sequence data - taxon names don't have digits/ACGT at start typically
            # Better heuristic: if line starts with a word that looks like a name (letters, underscores)
            # followed by whitespace and then sequence chars
            match = re.match(r'^(\S+)\s+', line)
            if match and match.group(1) and not all(c in 'ACGTacgt?- ' for c in match.group(1)):
                # Save previous
                if current_name and current_seq:
                    sequences[current_name] = current_seq

                # Start new taxon
                parts = line.split(None, 1)
                current_name = parts[0]
                current_seq = parts[1] if len(parts) > 1 else ""
            else:
                current_seq += line
        else:
            # First taxon
            parts = line.split(None, 1)
            if len(parts) >= 1:
                current_name = parts[0]
                current_seq = parts[1] if len(parts) > 1 else ""

    # Don't forget the last one
    if current_name and current_seq:
        sequences[current_name] = current_seq

    return sequences


def count_valid_bases(seq):
    """Count non-missing, non-gap characters in a sequence."""
    return sum(1 for c in seq if c not in '?- \t\n\r')


def process_all_nexus():
    """Main processing function."""
    nexus_files = sorted([f for f in os.listdir(NEXUS_DIR) if f.endswith('.nexus') and f.startswith('uce-')])

    if not nexus_files:
        print("ERROR: No nexus files found!", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(nexus_files)} nexus files to process.")

    total_processed = 0
    total_with_tityra = 0
    total_no_tityra = 0

    with open(OUTPUT_FASTA, 'w') as out:
        for filename in nexus_files:
            nexus_path = os.path.join(NEXUS_DIR, filename)
            locus = filename.replace('.nexus', '')

            # Parse nexus file
            sequences = parse_nexus(nexus_path)

            if not sequences:
                total_no_tityra += 1
                continue

            # Filter for Tityra sequences
            tityra_seqs = {name: seq for name, seq in sequences.items()
                          if TITYRA_PATTERN.match(name)}

            if not tityra_seqs:
                total_no_tityra += 1
                continue

            # Find the Tityra sequence with the most valid bases
            best_name = None
            best_count = -1
            best_seq = None

            for name, seq in tityra_seqs.items():
                # Remove whitespace from sequence
                clean_seq = seq.replace(' ', '').replace('\t', '').replace('\n', '').replace('\r', '')
                valid = count_valid_bases(clean_seq)
                if valid > best_count:
                    best_count = valid
                    best_name = name
                    best_seq = clean_seq

            if best_name is not None and best_seq is not None:
                # Remove gaps and missing characters
                clean_seq = best_seq.replace('-', '').replace('?', '')
                out.write(f">{locus}_{best_name}\n")
                # Write sequence in wrapped format (80 chars per line)
                for i in range(0, len(clean_seq), 80):
                    out.write(clean_seq[i:i+80] + '\n')
                total_with_tityra += 1
            else:
                total_no_tityra += 1

            total_processed += 1

            if total_processed % 100 == 0:
                print(f"  Processed {total_processed}/{len(nexus_files)}...")

    print(f"\nDone! Processed {total_processed} files.")
    print(f"  Loci with Tityra sequences: {total_with_tityra}")
    print(f"  Loci without Tityra sequences: {total_no_tityra}")
    print(f"Output written to: {OUTPUT_FASTA}")


if __name__ == "__main__":
    process_all_nexus()