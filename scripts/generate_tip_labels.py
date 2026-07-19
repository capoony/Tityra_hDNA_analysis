#!/usr/bin/env python3
"""Generate a tip-label mapping TSV for a phylogeny tree.

For each tree tip whose label matches the `Code` column (column 7) of the
provided xlsx sample sheet, build a human-readable label:

    Genus Species Subspecies

When more than one individual in the *tree* maps to the same
Genus/Species/Subspecies (or Genus Species when subspecies is missing), the
collection ID (xlsx column 3) is appended to keep labels unique. If the
collection ID is missing, the trailing numeric/code suffix from the original
tip label is used instead.

Outputs a two-column, tab-separated file:  old_tip_label<TAB>new_label
Lines for tips not appearing in the xlsx are omitted (they keep their original
label, only with underscores converted to spaces by the R script).
"""

import argparse
import sys
from collections import Counter

import openpyxl


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree", required=True, help="Newick tree file (tips = xlsx codes)")
    p.add_argument("--xlsx", required=True, help="Sample sheet xlsx")
    p.add_argument("--out", required=True, help="Output TSV (old<TAB>new)")
    return p.parse_args()


def read_tree_tips(tree_file):
    """Extract tip labels from a Newick file (single line expected)."""
    with open(tree_file) as fh:
        text = fh.read().strip()
    # Tips are leaf names immediately preceding a ':' or end of a clade.
    tips = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if c == '(':
            i += 1
            continue
        if c in ':,);':
            i += 1
            continue
        # accumulate a token until delimiter
        j = i
        while j < n and text[j] not in ':,);':
            j += 1
        token = text[i:j].strip()
        if token:
            tips.append(token)
        i = j
    return tips


def main():
    args = parse_args()

    # --- read xlsx mapping: code -> (base_name, coll_id) ---
    wb = openpyxl.load_workbook(args.xlsx, read_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    wb.close()

    code_to_info = {}
    for r in rows[1:]:
        # columns: Material, Institution, Tissue/Coll.#, Genus, Species,
        #          Subspecies, Code, ND2, cyt b, ...
        if len(r) < 7:
            continue
        genus, species, subspecies, code = r[3], r[4], r[5], r[6]
        coll = r[2]
        if code is None:
            continue
        code = str(code).strip()
        genus = str(genus).strip() if genus is not None else ""
        species = str(species).strip() if species is not None else ""
        if subspecies is not None and str(subspecies).strip():
            base = f"{genus} {species} {str(subspecies).strip()}"
        else:
            base = f"{genus} {species}".strip()
        code_to_info[code] = (base, coll)

    # --- intersect with tree tips ---
    tree_tips = read_tree_tips(args.tree)

    # Build list of (code, base, coll) for tips present in the tree
    present = []
    for tip in tree_tips:
        if tip in code_to_info:
            base, coll = code_to_info[tip]
            present.append((tip, base, coll))

    # Count how many times each base name occurs among present tips
    base_counts = Counter(b for _, b, _ in present)

    # --- write mapping (ensure globally unique labels) ---
    written = 0
    used = set()
    with open(args.out, "w") as out:
        for tip, base, coll in present:
            new_label = base
            if base_counts[base] > 1:
                # append collection id (or fallback to trailing token of tip)
                if coll is not None and str(coll).strip():
                    suffix = str(coll).strip()
                    # Disambiguate when the same collection id is shared by
                    # different vouchers: prepend the institution token from
                    # the tip. Institutions may appear either as a separate
                    # underscore token (e.g. ..._LSU_18568) or joined with the
                    # id (e.g. ..._LSUMZ_18568).
                    tokens = tip.split("_")
                    if suffix in tokens:
                        pos = tokens.index(suffix)
                        if pos > 0:
                            inst = tokens[pos - 1]
                            if inst:
                                suffix = f"{inst} {suffix}"
                    else:
                        for tok in tokens:
                            if tok == f"{tok.rsplit('_', 1)[0]}_{suffix}":
                                inst = tok.rsplit("_", 1)[0]
                                if inst:
                                    suffix = f"{inst} {suffix}"
                                break
                else:
                    # use last underscore-delimited token of the original tip
                    suffix = tip.split("_")[-1]
                new_label = f"{base} {suffix}"
            # Guarantee uniqueness: append a numeric index if still duplicated
            if new_label in used:
                k = 2
                while f"{new_label} {k}" in used:
                    k += 1
                new_label = f"{new_label} {k}"
            used.add(new_label)
            out.write(f"{tip}\t{new_label}\n")
            written += 1

    print(f"Wrote {written} tip-label mappings to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())