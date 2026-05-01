# Graph-Based Phylogeny Integration

## Overview

This document explains how the graph-based alignment pipeline integrates with the phylogeny workflow to avoid reference bias in ancient DNA analysis.

## Directory Structure

### Graph Alignment Outputs (Step 9)

Created by: `graph_alignment_busco.sh`

```
results/
├── graph_alignment/
│   ├── busco_refs/                    # Reference genomes
│   │   ├── Oxyruncus_cristatus.fna.gz
│   │   ├── Tyrannus_savana.fna.gz
│   │   └── Pachyramphus_minor.fna.gz
│   ├── busco_results/                 # BUSCO runs on references
│   │   ├── Oxyruncus_cristatus/
│   │   ├── Tyrannus_savana/
│   │   └── Pachyramphus_minor/
│   ├── busco_regions/                 # Extracted BUSCO CDS sequences
│   │   ├── 10000at8782.fasta
│   │   ├── 10001at8782.fasta
│   │   └── ...
│   ├── aligned/                       # MAFFT-aligned sequences
│   │   ├── 10000at8782_aligned.fasta
│   │   ├── 10001at8782_aligned.fasta
│   │   └── ...
│   ├── graphs/                        # Variation graphs
│   │   ├── 10000at8782.xg
│   │   ├── 10000at8782.gcsa
│   │   └── ...
│   ├── mapped/                        # Read mappings to graphs
│   │   └── *.gam files (large, cleaned after processing)
│   └── consensus/                     # Consensus sequences from graph
│       ├── 10000at8782_Tityra.fasta
│       ├── 10001at8782_Tityra.fasta
│       └── ...
│
└── phylogeny_graph/
    └── busco_aa/                      # Translated protein sequences
        ├── 10000at8782.faa
        ├── 10001at8782.faa
        └── ...
```

### Phylogeny Integration (Step 10)

Created by: `main.sh` Step 10

```
results/
└── phylogeny/
    ├── data/
    │   ├── genomes.txt                # List of reference genomes
    │   ├── genomes.names              # List of genome names
    │   ├── Pachyramphus_minor.fna.gz  # Reference genome
    │   ├── Tityra_leucura.fna.gz      # De novo assembly
    │   └── Tityra_leucura_graph/      # Graph-based sequences
    │       └── busco_set.fasta        # Combined BUSCO proteins
    ├── BUSCO/                         # BUSCO runs on all genomes
    │   ├── Pachyramphus_minor/
    │   ├── Tityra_leucura/
    │   └── Tityra_leucura_graph/      # Graph-based BUSCO structure
    ├── concatenated/                  # Shared BUSCO genes
    ├── prealigned/                    # Pre-alignment sequences
    ├── mafft/                         # Aligned sequences
    └── phylogeny/                     # Final phylogenetic trees
        ├── alignment.fa               # Concatenated alignment
        ├── RAxML_bestTree.Tityra
        ├── RAxML_bootstrap.bootrep
        ├── RAxML_bipartitions.FINAL
        ├── Tityra_BUSCO.pdf
        └── Tityra_BUSCO.png
```

### Individual Gene Phylogenies

Created by: `BUSCO_phylogeny_graph_indgenes.sh`

```
results/
└── phylogeny_graph/
    ├── busco_aa/                      # Protein sequences (from Step 9)
    ├── gene_list.txt                  # List of genes to process
    └── phylogeny/                     # Per-gene phylogenies
        ├── 10000at8782/
        │   ├── 10000at8782_combined.fasta
        │   ├── 10000at8782_aligned.fasta
        │   ├── RAxML_bestTree.Tityra_graph
        │   ├── RAxML_bootstrap.bootrep
        │   ├── RAxML_bipartitions.FINAL
        │   ├── 10000at8782_graph_BUSCO.pdf
        │   └── 10000at8782_graph_BUSCO.png
        └── ...
```

## Workflow

### Step 9: Graph-Based Alignment

Script: `graph_alignment_busco.sh`

1. **Download References**: Get 2-3 related reference genomes
2. **Run BUSCO**: Identify single-copy orthologs in references
3. **Extract CDS**: Get pre-extracted BUSCO .fna files (coding sequences)
4. **Align with MAFFT**: Create MSA for each BUSCO gene
5. **Build Graphs**: Construct variation graphs with vg
6. **Map Reads**: Map Tityra ancient DNA reads to graphs
7. **Call Variants**: Extract variants with vg pack/call
8. **Generate Consensus**: Create consensus sequences with bcftools
9. **Translate**: Convert nucleotide to protein sequences

**Output**: Protein sequences in `phylogeny_graph/busco_aa/`

### Step 10: Phylogeny with Graph Integration

Script: `main.sh` Step 10

1. **Prepare Data**:
   - Copy reference genomes
   - Copy de novo assembly
   - **Combine graph-based BUSCO proteins** into single file

2. **Run BUSCO**: On all genomes (references, de novo, graph-based)
   - For graph-based: Skip BUSCO run, use existing sequences
   - Create BUSCO-like directory structure for compatibility

3. **Identify Shared Genes**: Find BUSCO genes present in all genomes

4. **Concatenate**: Combine sequences for each shared gene

5. **Align**: MAFFT alignment of concatenated sequences

6. **Build Tree**: RAxML maximum likelihood phylogeny with bootstrap

7. **Visualize**: Plot tree with ggtree

**Output**: Combined phylogeny comparing:

- Reference genomes
- De novo assembly (potential contamination)
- **Graph-based sequences (reduced reference bias)**

### Individual Gene Phylogenies

Script: `BUSCO_phylogeny_graph_indgenes.sh`

1. **Check Availability**: Verify graph-based BUSCO sequences exist

2. **For Each Gene**:
   - Combine graph-based Tityra with reference sequences
   - Align with MAFFT
   - Build ML tree with RAxML
   - Bootstrap (25 replicates)
   - Visualize with ggtree

**Output**: Per-gene phylogenies in `phylogeny_graph/phylogeny/`

## Key Features

### Reference Bias Mitigation

- **Problem**: Mapping to single reference creates bias toward that reference
- **Solution**: Variation graph incorporates multiple references
- **Benefit**: Captures true variation, not reference artifacts

### Ancient DNA Considerations

- **Short reads**: 50bp median length (graphs handle this well)
- **Low coverage**: vg call handles sparse data
- **Damage patterns**: mapDamage analysis done separately

### Comparison Framework

Three phylogenetic approaches:

1. **Graph-based**: Reduced bias (new method)
2. **De novo assembly**: Independent but may have contamination
3. **Direct mapping**: Baseline with known bias

## Execution Order

1. **First run Step 10 regular phylogeny** (BUSCO on references)
2. **Then run Step 9** (graph alignment reuses BUSCO results)
3. **Step 10 automatically integrates** graph-based sequences
4. **Finally build individual gene trees** (optional, for detailed comparison)

## Error Handling

### Missing Graph Sequences

If Step 9 hasn't been run:

- Step 10 detects missing `phylogeny_graph/busco_aa/`
- Prints warning but continues without graph-based sequences
- Only includes de novo assembly and references

### Empty Consensus

If bcftools consensus fails:

- Check VCF is compressed with bgzip
- Verify tabix index exists
- Ensure reference FASTA is indexed with samtools faidx

### Missing BUSCO Genes

If gene missing in regular phylogeny:

- Graph phylogeny script skips that gene
- Continues with available genes
- Prints warning message

## Validation

Check successful integration:

```bash
# Verify graph-based BUSCO sequences exist
ls -lh ${WD}/results/phylogeny_graph/busco_aa/*.faa | wc -l

# Check they're included in phylogeny
grep "Tityra_leucura_graph" ${WD}/results/phylogeny/data/genomes.names

# Verify protein sequences are non-empty
for f in ${WD}/results/phylogeny_graph/busco_aa/*.faa; do
    wc -l $f
done

# Check final tree includes graph-based sample
grep "Tityra_leucura_graph" ${WD}/results/phylogeny/phylogeny/RAxML_bipartitions.FINAL
```

## Expected Results

### Phylogenetic Placement

- **Graph-based Tityra**: Should cluster with true Tityra position
- **De novo assembly**: May diverge if contaminated
- **Direct mapping**: Likely biased toward reference species

### Branch Lengths

- Graph-based should show realistic divergence
- Direct mapping may show artificially short branch to reference

## Troubleshooting

### Directory Not Found

**Error**: `phylogeny_graph does not exist`
**Solution**: Run Step 9 before Step 10

### Empty Sequences

**Error**: Consensus sequences are all N's
**Solution**: Check VCF processing (bgzip, tabix, samtools faidx)

### Translation Errors

**Error**: Protein sequences have stop codons
**Solution**: BUSCO CDS should be in-frame; check extraction

### Missing Dependencies

**Error**: `mafft: command not found`
**Solution**: Activate conda environment (mafft-7.487)

## Performance

### Time Requirements

- Step 9 (graph alignment): ~4-8 hours (depends on number of genes)
- Step 10 (phylogeny): ~6-12 hours (includes BUSCO runs)
- Individual gene trees: ~1-2 hours per gene

### Storage Requirements

- Graph files (*.xg,*.gcsa): ~100MB per gene
- Mapped reads (*.gam): ~500MB-1GB per gene (cleaned after processing)
- Consensus sequences: ~5KB per gene
- Total: ~50-100GB for complete pipeline

## References

- vg toolkit: <https://github.com/vgteam/vg>
- BUSCO: <https://busco.ezlab.org/>
- RAxML: <https://cme.h-its.org/exelixis/web/software/raxml/>

## Last Updated

April 8, 2026
