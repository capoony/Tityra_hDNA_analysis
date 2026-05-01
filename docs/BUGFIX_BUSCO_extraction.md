# Bug Fix: BUSCO Sequence Extraction for Graph Alignment

## Date

April 7, 2026

## Problem

The Python script in Step 9 (graph alignment) was extracting entire genomic regions (including introns) from reference genomes using coordinates from BUSCO's `full_table.tsv`. This resulted in:

- Extremely large files (7.7K to 424K for single genes)
- Mixed uppercase/lowercase sequences (introns vs exons)
- Sequences that wouldn't align properly for graph construction
- File: `busco_regions/10020at8782.fasta` was 424K (should be ~2-5Kb for CDS)

## Root Cause

The script was attempting to extract sequences using genomic coordinates:

```python
for record in SeqIO.parse(ref_file, "fasta"):
    if record.id == scaffold:
        seq = str(record.seq[start-1:end])  # Extracts entire genomic region
```

This extracted the full gene locus including:

- Exons (coding sequences)
- Introns (non-coding, shown as lowercase in soft-masked genomes)
- Regulatory regions
- Total size: tens to hundreds of kilobases

## Solution

BUSCO already extracts clean coding sequences during its analysis! These are stored in:

```
busco_refs/[ref_name]_busco/run_aves_odb10/busco_sequences/single_copy_busco_sequences/
├── [gene_id].fna  # nucleotide coding sequence (CDS)
├── [gene_id].faa  # protein sequence
└── [gene_id].gff  # annotation
```

The fix: Use the pre-extracted `.fna` files directly instead of extracting from genomic coordinates.

### Updated Code

```python
for busco_result in os.listdir(busco_dir):
    if not busco_result.endswith('_busco'):
        continue
    
    ref_name = busco_result.replace('_busco', '')
    seq_dir = f"{busco_dir}/{busco_result}/run_aves_odb10/busco_sequences/single_copy_busco_sequences"
    
    # Find all .fna (nucleotide coding sequence) files
    for fname in os.listdir(seq_dir):
        if fname.endswith('.fna'):
            gene_id = fname.replace('.fna', '')
            busco_genes.add(gene_id)

# For each gene, collect pre-extracted sequences
for gene in sorted(busco_genes):
    for ref_name, seq_dir in busco_refs.items():
        fna_file = f"{seq_dir}/{gene}.fna"
        
        # Read the BUSCO-extracted coding sequence
        for record in SeqIO.parse(fna_file, "fasta"):
            gene_seqs.append(f">{ref_name}\n{str(record.seq)}\n")
            break
```

## Results After Fix

```bash
Found 8222 BUSCO genes across 3 references
Created 11393at8782.fasta with 3 sequences
```

File sizes now appropriate:

- 11393at8782.fasta: 7.9K (was potentially 100K+)
- Sequences are clean coding sequences
- No mixed case (no introns)
- Suitable for multiple sequence alignment and graph construction

## Files Modified

1. `/media/inter/mkapun/projects/Tityra/shell/main.sh` (lines 428-480)
2. `/media/inter/mkapun/projects/Tityra/shell/graph_alignment_busco.sh` (lines 140-195)

## Impact

This fix enables:

- Proper graph construction with vg toolkit
- Aligned orthologous coding sequences across references
- Correct phylogenetic inference
- Reduced file sizes and processing time

## Testing

```bash
cd /media/inter/mkapun/projects/Tityra
rm -rf results/graph_alignment/busco_regions
mkdir -p results/graph_alignment/busco_regions

# Run corrected extraction
bash shell/main.sh  # Run Step 9

# Verify output
ls -lh results/graph_alignment/busco_regions/ | head
head -10 results/graph_alignment/busco_regions/11393at8782.fasta
```

Expected output: Multi-FASTA files with 2-3 sequences per gene, each ~1-10Kb CDS length.
