#!/usr/bin/env python3
"""
Alternative analysis: For each aligned UCE locus, remove columns where
Tityra_leucura has an N (ambiguous base). Then build supermatrix + partition file.
"""

import os
import sys
from collections import OrderedDict
from Bio import AlignIO

WD = "/media/inter/mkapun/projects/Tityra"
RESULTS = os.path.join(WD, "results/tityra_uce/phylogeny")
ALIGNMENTS_DIR = os.path.join(RESULTS, "alignments")
OUTPUT_DIR = os.path.join(RESULTS, "trimmed")
SPECIES_LIST = os.path.join(RESULTS, "species_list.txt")
LOCI_LIST = os.path.join(RESULTS, "loci_list.txt")

OUT_ALIGNMENTS = os.path.join(OUTPUT_DIR, "alignments")
SUPERMATRIX_DIR = os.path.join(OUTPUT_DIR, "supermatrix")
TITYRA_LEUCURA = "Tityra_leucura"


def main():
    os.makedirs(OUT_ALIGNMENTS, exist_ok=True)
    os.makedirs(SUPERMATRIX_DIR, exist_ok=True)

    # Load species list
    with open(SPECIES_LIST) as f:
        all_species = [l.strip() for l in f if l.strip()]

    # Load loci list
    with open(LOCI_LIST) as f:
        loci = [l.strip() for l in f if l.strip()]

    print(f"Processing {len(loci)} aligned loci, trimming N columns from Tityra_leucura...", file=sys.stderr)

    trimmed_loci = []
    for locus in loci:
        aln_file = os.path.join(ALIGNMENTS_DIR, f"{locus}_aln.fasta")
        if not os.path.isfile(aln_file):
            print(f"  SKIP {locus}: alignment not found", file=sys.stderr)
            continue

        try:
            aln = AlignIO.read(aln_file, "fasta")
        except Exception as e:
            print(f"  SKIP {locus}: cannot read ({e})", file=sys.stderr)
            continue

        aln_len = aln.get_alignment_length()

        # Build dict of sequences
        seqs = {}
        for rec in aln:
            name = rec.id.replace('_R_', '')
            seqs[name] = str(rec.seq)

        if TITYRA_LEUCURA not in seqs:
            print(f"  SKIP {locus}: Tityra_leucura not in alignment", file=sys.stderr)
            continue

        leucura = seqs[TITYRA_LEUCURA]

        # Find columns where leucura does NOT have N
        keep_cols = [i for i, base in enumerate(leucura) if base.upper() != 'N']

        if len(keep_cols) == 0:
            print(f"  SKIP {locus}: all columns trimmed (100% N in leucura)", file=sys.stderr)
            continue

        # Trim all sequences to keep only those columns
        trimmed = OrderedDict()
        for name in seqs:
            trimmed[name] = ''.join(seqs[name][i] for i in keep_cols)

        # Write trimmed alignment
        out_file = os.path.join(OUT_ALIGNMENTS, f"{locus}_trim.fasta")
        with open(out_file, 'w') as f:
            for name in all_species:
                seq = trimmed.get(name)
                if seq is None:
                    seq = '-' * len(keep_cols)
                f.write(f">{name}\n{seq}\n")

        trimmed_loci.append(locus)

        if len(trimmed_loci) % 200 == 0:
            print(f"  Processed {len(trimmed_loci)} loci...", file=sys.stderr)

    print(f"\nTrimmed {len(trimmed_loci)} loci (removed N-columns from Tityra_leucura)", file=sys.stderr)

    if len(trimmed_loci) == 0:
        print("ERROR: no loci survived trimming!", file=sys.stderr)
        sys.exit(1)

    # Build supermatrix + partition file
    print("Building supermatrix...", file=sys.stderr)
    supermatrix = OrderedDict((sp, []) for sp in all_species)
    partitions = []
    pos = 1
    skipped = 0

    for locus in trimmed_loci:
        aln_file = os.path.join(OUT_ALIGNMENTS, f"{locus}_trim.fasta")
        try:
            aln = AlignIO.read(aln_file, "fasta")
        except:
            skipped += 1
            continue

        aln_len = aln.get_alignment_length()
        aln_dict = {rec.id: str(rec.seq) for rec in aln}

        for sp in all_species:
            seq = aln_dict.get(sp, "-" * aln_len)
            supermatrix[sp].append(seq)

        partitions.append(f"DNA, {locus}_trim = {pos}-{pos + aln_len - 1}")
        pos += aln_len

    # Remove all-gap species
    to_remove = [sp for sp, seqs in supermatrix.items()
                 if all(c == '-' for c in ''.join(seqs))]
    for sp in to_remove:
        del supermatrix[sp]
        print(f"  REMOVE {sp}: all gaps after trimming", file=sys.stderr)

    # Write supermatrix
    concat_path = os.path.join(SUPERMATRIX_DIR, "supermatrix_trimmed.fasta")
    with open(concat_path, 'w') as f:
        for sp, seqs in supermatrix.items():
            f.write(f">{sp}\n{''.join(seqs)}\n")

    # Write partition file
    part_path = os.path.join(SUPERMATRIX_DIR, "partitions_trimmed.txt")
    with open(part_path, 'w') as f:
        f.write("\n".join(partitions) + "\n")

    print(f"\n=== Alternative analysis (N-columns trimmed) ===", file=sys.stderr)
    print(f"  Loci used:       {len(trimmed_loci) - skipped}", file=sys.stderr)
    print(f"  Species:         {len(supermatrix)} (removed {len(to_remove)} all-gap)", file=sys.stderr)
    print(f"  Supermatrix len: {pos - 1} nt", file=sys.stderr)
    print(f"  Partitions:      {len(trimmed_loci) - skipped}", file=sys.stderr)
    print(f"  Output:          {concat_path}", file=sys.stderr)


if __name__ == "__main__":
    main()