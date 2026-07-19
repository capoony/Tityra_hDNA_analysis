#!/usr/bin/env python3
"""
Build supermatrix and partition file from aligned UCE loci.
"""

import os
import sys
from collections import OrderedDict

WD = "/media/inter/mkapun/projects/Tityra"
RESULTS = os.path.join(WD, "results/tityra_uce/phylogeny")
GENE_SEQS_DIR = os.path.join(RESULTS, "gene_sequences")
ALIGNMENTS_DIR = os.path.join(RESULTS, "alignments")
SUPERMATRIX_DIR = os.path.join(RESULTS, "supermatrix")
SPECIES_LIST = os.path.join(RESULTS, "species_list.txt")
LOCI_LIST = os.path.join(RESULTS, "loci_list.txt")


def main():
    os.makedirs(SUPERMATRIX_DIR, exist_ok=True)

    # Load species list
    with open(SPECIES_LIST) as f:
        all_species = [l.strip() for l in f if l.strip()]

    # Load loci list
    with open(LOCI_LIST) as f:
        loci = [l.strip() for l in f if l.strip()]

    print(f"Building supermatrix from {len(loci)} loci, {len(all_species)} species...", file=sys.stderr)

    # Build supermatrix
    supermatrix = OrderedDict((sp, []) for sp in all_species)
    partitions = []
    pos = 1
    skipped = 0

    for locus in loci:
        aln_file = os.path.join(ALIGNMENTS_DIR, f"{locus}_aln.fasta")
        if not os.path.isfile(aln_file):
            # Try gene_sequences dir as fallback (unaligned)
            aln_file = os.path.join(GENE_SEQS_DIR, f"{locus}.fasta")
            if not os.path.isfile(aln_file):
                print(f"  SKIP {locus}: file not found", file=sys.stderr)
                skipped += 1
                continue

        try:
            from Bio import AlignIO
            aln = AlignIO.read(aln_file, "fasta")
        except Exception as e:
            print(f"  SKIP {locus}: cannot read ({e})", file=sys.stderr)
            skipped += 1
            continue

        aln_len = aln.get_alignment_length()
        aln_dict = {rec.id.replace('_R_', ''): str(rec.seq) for rec in aln}

        for sp in all_species:
            # Try with and without _R_ prefix
            seq = aln_dict.get(sp)
            if seq is None:
                seq = aln_dict.get(f"_R_{sp}")
            if seq is None:
                seq = "-" * aln_len
            supermatrix[sp].append(seq)

        # One partition per locus (DNA model, will be optimized by IQ-TREE MFP)
        partitions.append(f"DNA, {locus} = {pos}-{pos + aln_len - 1}")
        pos += aln_len

    # Filter out species with all gaps
    to_remove = []
    for sp, seqs in supermatrix.items():
        concat = ''.join(seqs)
        if all(c == '-' for c in concat):
            to_remove.append(sp)

    for sp in to_remove:
        del supermatrix[sp]
        print(f"  REMOVE {sp}: all gaps", file=sys.stderr)

    # Write supermatrix FASTA
    concat_path = os.path.join(SUPERMATRIX_DIR, "supermatrix.fasta")
    with open(concat_path, 'w') as f:
        for sp, seqs in supermatrix.items():
            concat = ''.join(seqs)
            f.write(f">{sp}\n{concat}\n")

    # Write partition file
    part_path = os.path.join(SUPERMATRIX_DIR, "partitions.txt")
    with open(part_path, 'w') as f:
        f.write("\n".join(partitions) + "\n")

    used_loci = len(loci) - skipped

    print(f"\nDone!", file=sys.stderr)
    print(f"  Loci used:       {used_loci}", file=sys.stderr)
    print(f"  Loci skipped:    {skipped}", file=sys.stderr)
    print(f"  Species:         {len(supermatrix)} (removed {len(to_remove)} all-gap)", file=sys.stderr)
    print(f"  Supermatrix len: {pos - 1} nt", file=sys.stderr)
    print(f"  Partitions:      {used_loci}", file=sys.stderr)
    print(f"  Output:          {concat_path}", file=sys.stderr)


if __name__ == "__main__":
    main()