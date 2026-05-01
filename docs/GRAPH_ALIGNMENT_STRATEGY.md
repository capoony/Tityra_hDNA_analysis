# Graph-Based Alignment Strategy for Tityra Phylogenomics

## Problem Statement

Low-coverage ancient DNA (median 50bp read length) from Tityra leucura (unknown species) needs to be aligned to related reference genome (Pachyramphus minor) for BUSCO gene extraction and phylogenetic analysis. **Reference bias** is a major concern that could distort phylogenetic placement.

## Solution: Pangenome Variation Graphs

### Concept

Instead of mapping reads to a single linear reference genome (which introduces bias toward that reference), we:

1. **Build variation graphs** from multiple closely-related reference genomes
2. **Map reads to the graph structure**, allowing reads to follow paths that best match the unknown sample
3. **Extract consensus sequences** that represent the actual Tityra leucura sequence, not biased toward any single reference

### Implementation in Your Pipeline

The graph-based approach has been integrated into your pipeline in **Step 9**, implemented in two ways:

#### Option 1: Inline Implementation (currently in main.sh lines 327-546)

- Embedded directly in the main pipeline
- Uses vg toolkit for graph construction and alignment
- Processes BUSCO gene regions specifically (computational efficiency)

#### Option 2: Modular Script (graph_alignment_busco.sh)

- Cleaner, reusable implementation
- Can be called with: `bash ${WD}/shell/graph_alignment_busco.sh ${WD} 50`
- Better for testing and debugging

### Workflow Steps

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Prepare Multiple Reference Genomes                      │
│    - Pachyramphus minor (primary, closest relative)        │
│    - Oxyruncus cristatus (family Cotingidae)              │
│    - Tyrannus savana (family Tyrannidae)                  │
│    - Piprites chloris (family Pipritidae)                 │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. Run BUSCO on Each Reference                             │
│    - Identify complete single-copy ortholog genes          │
│    - Extract coordinates for each BUSCO gene               │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. Extract BUSCO Gene Sequences from References            │
│    - For each BUSCO gene present in ≥2 references          │
│    - Create multi-FASTA with all reference sequences       │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. Build Variation Graphs (vg construct)                   │
│    - One graph per BUSCO gene                              │
│    - Graph encodes variation between references            │
│    - Index for efficient read mapping                      │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 5. Map Tityra Reads to Graphs (vg map)                     │
│    - Paired-end reads + merged reads                       │
│    - Reads follow optimal path through graph               │
│    - No single-reference bias                              │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 6. Call Variants and Extract Consensus (vg call)           │
│    - Identify Tityra-specific variants                     │
│    - Generate consensus sequence per BUSCO gene            │
│    - Convert to protein sequences                          │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 7. Integrate with Phylogeny Pipeline (Step 10)             │
│    - Add "Tityra_leucura_graph" as new sample              │
│    - Compare with de novo assembly-based sequences         │
│    - Run RAxML phylogenetic analysis                       │
└─────────────────────────────────────────────────────────────┘
```

### Key Advantages

1. **Reduced Reference Bias**: Reads can align to paths representing variation across multiple references, not forced to match a single genome

2. **Better for Divergent Sequences**: Graph structure accommodates polymorphisms and structural variation naturally

3. **Quality Control**: Compare graph-based vs. de novo assembly vs. direct mapping results to assess bias

4. **Validation**: In final phylogeny, "Tityra_leucura_graph" should cluster with "Tityra_leucura" (de novo). If they differ significantly, this indicates reference bias in one approach.

### Expected Outputs

```
results/
├── graph_alignment/
│   ├── references/                     # Multiple reference genomes
│   ├── busco_refs/                     # BUSCO results for each reference
│   ├── busco_regions/                  # Extracted BUSCO gene sequences
│   │   └── [gene_id].fasta            # Multi-FASTA per gene
│   ├── graphs/                         # Variation graphs
│   │   ├── [gene_id].xg               # Indexed graph
│   │   └── [gene_id].gcsa             # GCSA index for mapping
│   ├── mapped/                         # Read alignments
│   │   ├── [gene_id].gam              # Graph alignment
│   │   ├── [gene_id].vcf              # Called variants
│   │   └── [gene_id].bam              # BAM for QC
│   └── consensus/                      # Final sequences
│       └── [gene_id]_Tityra.fasta     # Consensus per gene
│
├── phylogeny_graph/
│   └── busco_aa/                       # Protein sequences for phylogeny
│       └── [gene_id].faa              # One file per BUSCO gene
│
└── phylogeny/
    ├── BUSCO/
    │   ├── Tityra_leucura/            # De novo assembly BUSCOs
    │   └── Tityra_leucura_graph/      # Graph-based BUSCOs
    └── phylogeny/
        ├── alignment.fa                # Concatenated alignment
        ├── RAxML_bestTree.Tityra      # Best ML tree
        └── Tityra_BUSCO.png           # Final phylogeny with both samples
```

### Tools Required

- **vg** (v1.59.0+): Variation graph toolkit - auto-downloaded if not present
- **BUSCO** (v5.4.3): Single-copy ortholog identification
- **samtools/bcftools**: Consensus calling
- **BioPython**: Sequence extraction (in Python scripts)

### Validation Strategy

Compare three approaches:

1. **Graph-based** (Step 9): `Tityra_leucura_graph`
   - Least reference bias
   - Should be most accurate for phylogenetics

2. **De novo assembly** (Step 8): `Tityra_leucura`
   - No reference bias
   - May be fragmented due to low coverage

3. **Direct mapping** (other scripts): `Tityra_leucura_mapped`
   - Most reference bias
   - Expected to cluster artifactually close to reference

### Interpreting Phylogeny Results

**Ideal outcome**:

```
├── Outgroup (Acanthisitta)
└── Ingroup
    ├── Reference (Pachyramphus_minor)
    └── Tityra clade
        ├── Tityra_leucura_graph    ← Graph-based (unbiased)
        └── Tityra_leucura          ← De novo assembly
```

**If seeing reference bias**:

```
├── Outgroup
└── Ingroup
    ├── Pachyramphus_minor
    │   └── Tityra_leucura_mapped   ← Too close = bias!
    └── True Tityra clade
        ├── Tityra_leucura_graph
        └── Tityra_leucura
```

### Performance Considerations

- **Runtime**: ~1-2 hours per BUSCO gene for graph construction + mapping (parallelizable)
- **Memory**: ~10-50 GB RAM per graph (depending on gene length and number of references)
- **Disk**: ~100-500 GB for all graphs and alignments

### Recommendations

1. **Start with modular script**: Run `graph_alignment_busco.sh` independently first
2. **Test on subset**: Process 10-20 BUSCO genes initially to validate approach
3. **Compare results**: Check concordance between graph-based and de novo before full run
4. **Monitor coverage**: Ensure adequate read depth for graph-based consensus (recommend ≥5x per gene)

### Troubleshooting

**Issue**: Few BUSCO genes extracted from graph

- **Solution**: Lower minimum reference count (use genes in ≥1 ref instead of ≥2)

**Issue**: Graph construction fails

- **Solution**: Check reference sequence formatting, ensure BUSCO coordinates are valid

**Issue**: Low mapping rate to graphs  

- **Solution**: Relax vg map parameters, check read quality, verify graph construction

### Alternative Approaches (if graph method not feasible)

1. **Iterative consensus mapping**:

   ```bash
   # Round 1: Map to Pachyramphus → call consensus
   # Round 2: Map to round 1 consensus → call consensus
   # Round 3: Map to round 2 consensus → final
   ```

2. **Competitive mapping**:

   ```bash
   # Map to 3 references simultaneously
   # Keep only reads with similar mapping quality to all
   ```

3. **Assembly-only**: Skip mapping entirely, use only de novo assembly BUSCO genes

## Summary

The graph-based approach provides a **robust, innovative solution** to reference bias in phylogenomics of divergent species. By encoding variation from multiple references, it allows your Tityra sequences to be accurately reconstructed without being artificially pulled toward any single reference genome.

This is particularly important for your 50bp aDNA reads where:

- Assembly may be fragmented (low coverage)
- Direct mapping would introduce strong bias
- Graph approach balances both concerns

The integration with your existing pipeline (Step 9→10) allows direct comparison of all three methods, providing strong validation of phylogenetic placement.
