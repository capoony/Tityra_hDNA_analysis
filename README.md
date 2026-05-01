# *Tityra leucura* Mitochondrial Phylogenomics Pipeline

This pipeline provides an end-to-end workflow for analyzing historical DNA from *Tityra leucura* museum specimens, focusing on mitochondrial genome assembly and phylogenetic analysis. The pipeline processes raw sequencing data through quality control, contamination assessment, mitochondrial genome assembly, and multi-gene phylogenetic reconstruction.

## Overview

The pipeline extracts and analyzes mitochondrial DNA from degraded historical specimens to:

- Assess DNA quality and contamination levels
- Reconstruct mitochondrial genes
- Build phylogenetic trees using multiple mitochondrial markers
- Compare *Tityra leucura* with related Tyrannidae species

## Quick Start

### Running the Complete Pipeline

```bash
WD=/media/inter/mkapun/projects/Tityra
bash ${WD}/shell/main.sh
```

### Pipeline Scripts

- **Main pipeline**: [`shell/main.sh`](shell/main.sh) - Master orchestration script
- **Mitochondrial assembly**: [`shell/mitofinder_assembly.sh`](shell/mitofinder_assembly.sh)
- **Gene combination**: [`shell/combine_mito_genes_assembly.sh`](shell/combine_mito_genes_assembly.sh)
- **Phylogenetic analysis**: [`shell/phylogenetic_analysis_assembly.sh`](shell/phylogenetic_analysis_assembly.sh)
- **Gene concatenation**: [`shell/concatenate_genes_assembly.sh`](shell/concatenate_genes_assembly.sh)
- **Concatenated tree**: [`shell/phylogenetic_analysis_concatenated_assembly.sh`](shell/phylogenetic_analysis_concatenated_assembly.sh)

Archived scripts from previous analyses are in [`shell/old/`](shell/old/).

## Pipeline Steps

### 1. Raw Data Preparation

Raw sequencing files (Illumina paired-end, 2×150 bp) are copied from the central repository to the working directory.

**Input**: Raw FASTQ files  
**Output**: `data/Illumina/`

---

### 2. Read Trimming and Quality Control

**Tool**: fastp v0.23.4

Quality trimming, adapter removal, and read merging are performed with two protocols:

- **Standard trimming**: Basic quality filtering (Q20), adapter removal, deduplication
- **Enhanced trimming**: Additional 5' end trimming to remove DNA damage artifacts

**Command**:

```bash
fastp -i R1.fastq.gz -I R2.fastq.gz \
    --merge --length_required 25 --dedup \
    --trim_poly_g --detect_adapter_for_pe
```

**Outputs**:

- Trimmed reads: `data/trimmed/Tityra_leucura_*_trimmed.fastq.gz`
- Merged reads: `data/trimmed/Tityra_leucura_merged.fastq.gz`
- QC report: [`data/trimmed/Tityra_leucura.html`](data/trimmed/Tityra_leucura.html)

**Interpretation**: Review HTML report for adapter content, quality scores, and insert size distribution. Typical results: ~124M raw reads, 93.8% passing filters, 36.6% successfully merged.

---

### 3. Taxonomic Classification

**Tool**: Kraken2 v2.1.2 (PlusPFP database)

Reads are classified against a comprehensive taxonomic database to identify contamination sources.

**Command**:

```bash
kraken2 --db pluspfp_20240904 --threads 150 \
    --report report.txt --paired R1 R2
```

**Outputs**:

- Classification reports: `results/kraken2/report_*.txt`
- Summary CSV: [`results/kraken2/kraken_summary.csv`](results/kraken2/kraken_summary.csv)

**Interpretation**: Typical contamination profile:

- 26% of reads classified
- Human: 2.1% (lab contamination)
- Bacteria: 14.8% (environmental)
- Fungi: 2.0% (*Penicillium*, *Aspergillus*)
- Target organism: 73.99% unclassified (not in database)

---

### 4. ECMSD Metagenomic Analysis

**Tool**: ECMSD pipeline

Detailed metagenomic profiling targeting mitochondrial DNA to assess contamination and read length distributions.

**Command**:

```bash
bash ECMSD.sh --fwd R1 --rev R2 --merged merged \
    --threads 200 --Binsize 1000
```

**Outputs**:

- Mapping statistics: `results/ECMSD/mapping/*.txt`
- Read length plots: [`results/ECMSD/mapping/Mito_summary_genus_ReadLengths.png`](results/ECMSD/mapping/Mito_summary_genus_ReadLengths.png)

**Interpretation**: Identifies *Pachyramphus* (Tyrannidae) as ninth-most abundant taxon, confirming endogenous DNA presence. Short read lengths typical for historical DNA.

---

### 5. Reference Genome Mapping

**Tool**: minimap2 v2.24, SAMtools v1.15

Reads are mapped to the closest available reference (*Tityra cayana*, NCBI: GCA_013397135.1) to assess genome coverage.

**Command**:

```bash
minimap2 -ax sr --secondary=no -t 200 reference.fna.gz R1 R2 | \
    samtools view -bS -F 4 - | \
    samtools sort -o output.bam
```

**Outputs**:

- BAM files: `results/minimap2/Tityra_leucura*.bam`
- Coverage statistics: [`results/minimap2/Tityra_leucura.coverage.txt`](results/minimap2/Tityra_leucura.coverage.txt)
- Coverage plot: [`results/minimap2/Tityra_leucura.coverage_plot.png`](results/minimap2/Tityra_leucura.coverage_plot.png)

**Interpretation**:

- 18.5M reads mapped
- Genome coverage: 17.29%
- Mean depth: 0.91×
- Base quality: Q39.5

Low coverage typical for degraded historical specimens.

---

### 6. DNA Damage Assessment

**Tool**: mapDamage2 v2.2.1

Ancient DNA damage patterns (cytosine deamination, depurination) are quantified to confirm historical origin.

**Command**:

```bash
mapDamage -i input.bam -r reference.fna.gz \
    --merge-reference-sequences
```

**Outputs**:

- Damage plots: [`results/mapDamage/Tityra_leucura/Fragmisincorporation_plot.png`](results/mapDamage/Tityra_leucura/Fragmisincorporation_plot.png)
- Statistics: `results/mapDamage/Tityra_leucura/misincorporation.txt`

**Interpretation**:

- C→T transitions at 5' ends: DNA deamination
- G→A transitions at 3' ends: complementary strand deamination
- δS = 0.014 (damage signal)
- Mean fragment length: 156 bp
- Confirms historical DNA origin

---

### 7. Contaminant Removal

**Process**: Reference-based filtering

Reads are mapped against known contaminant genomes. Only unmapped reads are retained for downstream analysis.

**Contaminant references**:

- *Homo sapiens* (GCA_000001405.15)
- *Penicillium coprophilum* (GCA_001890745.1)
- *Vanrija pseudolonga* (GCA_024498055.1)
- *Malassezia restricta* (GCA_000286055.1)
- *Aspergillus cristatus* (GCA_003344705.1)

**Command**:

```bash
minimap2 -ax sr --secondary=no -t 200 contaminants.fna.gz reads.fq.gz | \
    awk '($5 < 20 || and($2,4)) || $1 ~ /^@/' | \
    samtools view -bS - > filtered.bam
```

**Outputs**:

- Cleaned reads: `data/trimmed2/Tityra_leucura_contaminants_merged_unmapped.fastq.gz`
- Mapping statistics: `results/contaminants/mappings/*.bam`

**Interpretation**: ~181K reads identified as contaminants and removed.

---

### 8. Mitochondrial Genome Assembly

**Tool**: MitoFinder v1.4.1 (Singularity container)

Mitochondrial genome is assembled from cleaned reads using a reference-guided approach.

**Command**:

```bash
mitofinder -j Tityra_leucura -a assembly.fasta \
    -r reference.gb -o aves --adjust-direction
```

**Outputs**:

- Annotated mitochondrial genome: [`results/mitofinder_assembly/Tityra_leucura_mitogenome.fa`](results/mitofinder_assembly/)
- GenBank format: `results/mitofinder_assembly/Tityra_leucura_mitogenome.gb`
- Assembly statistics: `results/mitofinder_assembly/Tityra_leucura_mitofinder_output/`

**Interpretation**:

- Assembly time: ~6 minutes (from genome assembly)
- Contigs: 3 (15+2+0 genes identified)
- Total length: ~31.5 kb
- 12 protein-coding genes extracted

---

### 9. Gene Sequence Retrieval

**Tool**: EfetchThePython

Additional mitochondrial gene sequences from related Tyrannidae species are downloaded from GenBank to provide phylogenetic context.

**Taxa sampled** (8 genera):

- *Onychorhynchus coronatus* (outgroup, royal flycatcher)
- *Tityra* (3 species)
- *Schiffornis*
- *Laniocera*
- *Iodopleura*
- *Laniisoma*
- *Xenopsaris*
- *Pachyramphus*

**Genes downloaded** (18 mitochondrial):
COI, CO2, CO3, ND1, ND2, ND3, ND4, ND4L, ND5, ND6, CYTB, ATP6, ATP8, COX1, COX2, COX3, rrnL (16S), rrnS (12S)

**Command**:

```bash
python EfetchThePython.py --email your@email.com \
    --Term '"Tityra"[Organism], COI' \
    --Output COI_Tityra --FASTA
```

**Outputs**:

- Gene sequences: `data/Mito/*_*.fasta`
- Total: 1,247 sequences retrieved

---

### 10. Gene Sequence Integration

**Script**: `combine_mito_genes_assembly.sh`

Downloaded GenBank sequences are combined with *Tityra leucura* mitochondrial genes extracted from MitoFinder assembly.

**Process**:

1. Parse MitoFinder output for gene sequences
2. Combine with GenBank downloads
3. Handle gene name synonyms (COI/COX1, CO2/COX2, CO3/COX3)
4. Remove duplicate sequences

**Outputs**:

- Combined gene FASTA files: `results/combined_mito_genes_assembly/*.fasta`
- 12 protein-coding genes prepared for alignment

---

### 11. Individual Gene Phylogenies

**Tools**: MAFFT v7.487, Gblocks v0.91b, IQ-TREE v2.3.3

**Process**:

1. **Multiple sequence alignment** (MAFFT):

```bash
mafft --auto --adjustdirection input.fasta > aligned.fasta
```

- Algorithm: Auto-selection (FFT-NS-2 for <200 seqs, L-INS-i for <50 seqs)
- Direction adjustment: Handles reverse complement sequences
- Conserved blocks extracted with Gblocks

1. **Phylogenetic inference** (IQ-TREE):

```bash
iqtree2 -s aligned_gblocks.fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO
```

- Model selection: ModelFinder evaluates 286 DNA models
- Bootstrap: 1,000 ultrafast bootstrap replicates
- Branch support: SH-aLRT test with 1,000 replicates

1. **Tree visualization** (R/ggtree):

```R
tree <- read.tree("gene.treefile")
tree <- root(tree, outgroup="Onychorhynchus_coronatus")
ggtree(tree, layout="roundrect") + geom_tiplab() + geom_nodelab()
```

**Outputs**:

- Alignments: `results/phylogenetic_analysis_assembly/alignments/*_aligned_gblocks.fasta`
- Trees: `results/phylogenetic_analysis_assembly/trees/*.treefile`
- Plots: `results/phylogenetic_analysis_assembly/plots/*.pdf`

**Genes analyzed** (12):
ATP6, ATP8, CO2, CO3, COI, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5

**Interpretation**:

- Individual gene trees show varying phylogenetic signal
- Bootstrap support values indicate confidence in each node
- *Tityra* species form well-supported clade
- Gene trees useful for detecting recombination or horizontal transfer

---

### 12. Concatenated Multi-Gene Phylogeny

**Script**: `concatenate_genes_assembly.sh`

**Strategy**: Species-level concatenation with flexible voucher selection

To maximize taxonomic sampling while maintaining data quality, different voucher specimens can contribute different genes for each species.

**Concatenation criteria**:

1. **Gene selection**: 4 mitochondrial protein-coding genes
   - COI (cytochrome c oxidase I)
   - CYTB (cytochrome b)
   - ND2 (NADH dehydrogenase 2)
   - CO2 (cytochrome c oxidase II)

2. **Species inclusion filters**:
   - Genus priority: ≥1 species per genus (ensures broad taxonomic coverage)
   - Data completeness: ≥50% gene presence (≥2/4 genes per species)
   - Sequence quality: Fewest gaps/ambiguities selected when multiple sequences available

3. **Voucher provenance tracking**:
   - Each gene tracks which GenBank accession was used
   - FASTA headers include full provenance: `>Species genes=X/4 COI:voucher|accession CYTB:voucher|accession`
   - Missing genes filled with gap characters

**Process**:

```python
# For each species:
for species in all_species:
    for gene in [COI, CYTB, ND2, CO2]:
        sequences = get_all_sequences(species, gene)
        best_seq = select_best_quality(sequences)  # fewest gaps
        accession = extract_accession(best_seq.description)
        concat_parts.append((gene, best_seq, accession))
```

**Outputs**:

- Concatenated alignment: [`results/concatenated_genes_assembly/concatenated_all_genes.fasta`](results/concatenated_genes_assembly/)
- Species summary: Lists which genes present for each taxon

**Example header**:

```
>Tityra_cayana genes=4/4 COI:ZMUC141678|MN356369.1 CYTB:ZMUC141678|MN356420.1 ND2:ZMUC141678|MN356495.1 CO2:ZMUC141678|MN356445.1
```

---

### 13. Partitioned Phylogenetic Analysis

**Script**: `phylogenetic_analysis_concatenated_assembly.sh`

**Model**: Partitioned analysis with gene-specific substitution models

**Partition file**:

```
DNA, COI = 1-651
DNA, CYTB = 652-1673
DNA, ND2 = 1674-2714
DNA, CO2 = 2715-3398
```

**IQ-TREE command**:

```bash
iqtree2 -s concatenated_all_genes.fasta \
    -p partitions.txt \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO
```

**Parameters**:

- Partitioned model: Each gene gets independent substitution model
- ModelFinder: Tests 286 models per partition (1,144 total)
- Ultrafast bootstrap: 1,000 replicates
- Branch support: SH-aLRT test

**Tree rooting and visualization**:

```R
# Root with Onychorhynchus coronatus (Tyrannidae, royal flycatcher)
tree <- root(tree, outgroup="Onychorhynchus_coronatus", resolve.root=TRUE)
tree <- ladderize(tree)

# Visualization with rounded edges
ggtree(tree, layout="roundrect", ladderize=TRUE, right=TRUE) +
    geom_nodepoint(color="darkblue", size=2, shape=21, fill="lightblue") +
    geom_tippoint(aes(subset=is_tityra), color="red", size=3, shape=17) +
    geom_tippoint(aes(subset=is_outgroup), color="darkgreen", size=3, shape=15)
```

**Outputs**:

- Partitioned tree: [`results/phylogenetic_analysis_concatenated_assembly/trees/concatenated_all_genes.treefile`](results/phylogenetic_analysis_concatenated_assembly/trees/)
- Model report: `concatenated_all_genes.iqtree`
- Tree plots:
  - Rectangular layout: `plots/concatenated_tree_rectangular.pdf`
  - Circular layout: `plots/concatenated_tree_circular.pdf`
  - Bootstrap support: `plots/concatenated_tree_support.pdf`

**Final dataset**:

- Taxa: 24 species across 8 genera
- Total alignment length: 3,398 bp
- Missing data: 28.84% (due to incomplete gene coverage)
- Parsimony-informative sites: Variable by partition

**Interpretation**:

- Concatenated analysis provides robust phylogenetic hypothesis
- Partitioned models account for different evolutionary rates across genes
- Bootstrap support values ≥70 indicate well-supported clades
- Voucher provenance allows reproducibility and data validation

---

## Output Directory Structure

```
Tityra/
├── data/
│   ├── Illumina/              # Raw sequencing reads
│   ├── trimmed/               # Quality-trimmed reads
│   ├── trimmed2/              # Contaminant-filtered reads
│   ├── ref/                   # Reference genomes
│   ├── contaminants/          # Contaminant reference genomes
│   └── Mito/                  # Downloaded mitochondrial genes
├── results/
│   ├── kraken2/               # Taxonomic classification
│   ├── ECMSD/                 # Metagenomic profiling
│   ├── minimap2/              # Reference mapping & coverage
│   ├── mapDamage/             # DNA damage patterns
│   ├── contaminants/          # Contaminant mapping
│   ├── mitofinder_assembly/   # Mitochondrial genome assembly
│   ├── combined_mito_genes_assembly/      # Combined gene sequences
│   ├── phylogenetic_analysis_assembly/    # Individual gene trees
│   ├── concatenated_genes_assembly/       # Concatenated alignment
│   └── phylogenetic_analysis_concatenated_assembly/  # Final phylogeny
├── shell/
│   ├── main.sh                            # Master pipeline script
│   ├── mitofinder_assembly.sh             # Mitochondrial assembly
│   ├── combine_mito_genes_assembly.sh     # Gene combination
│   ├── phylogenetic_analysis_assembly.sh  # Individual trees
│   ├── concatenate_genes_assembly.sh      # Gene concatenation
│   ├── phylogenetic_analysis_concatenated_assembly.sh  # Final tree
│   └── old/                               # Archived scripts
├── docs/
│   └── Methods_and_Results.md             # Detailed documentation
├── logs/                                  # Pipeline log files
└── scripts/                               # Utility scripts
```

---

## Key Results Summary

### Data Quality

- **Raw reads**: 123.9M reads (2×150 bp)
- **Post-QC**: 116.2M reads (93.78%)
- **Merged reads**: 45.3M (36.56%)
- **Mean quality**: Q36.8

### Contamination Assessment

- **Classified reads**: 26.01%
  - Human: 2.12%
  - Bacteria: 14.79%
  - Fungi: 2.0%
- **Target DNA**: 73.99% (unclassified, presumably *Tityra*)

### Reference Mapping

- **Mapped reads**: 18.5M
- **Genome coverage**: 17.29%
- **Mean depth**: 0.91×
- **Base quality**: Q39.5

### DNA Damage

- **C→T at 5' end**: 3.2% (deamination signal)
- **Damage index (δS)**: 0.014
- **Mean fragment**: 156 bp
- **Conclusion**: Authentic historical DNA

### Mitochondrial Assembly

- **Assembly method**: MitoFinder (reference-guided)
- **Contigs**: 3 contigs
- **Genes identified**: 12 protein-coding genes
- **Total length**: ~31.5 kb

### Phylogenetic Results

- **Individual gene trees**: 12 mitochondrial genes
- **Concatenated dataset**: 4 genes, 24 taxa, 3,398 bp
- **Outgroup**: *Onychorhynchus coronatus*
- **Support**: Ultrafast bootstrap + SH-aLRT
- **Topology**: *Tityra* forms well-supported clade within Tyrannidae

---

## Software Requirements

### Core Tools

- **fastp** v0.23.4 - Quality control
- **Kraken2** v2.1.2 - Taxonomic classification
- **minimap2** v2.24 - Read mapping
- **SAMtools** v1.15 - BAM manipulation
- **mapDamage2** v2.2.1 - DNA damage analysis
- **MitoFinder** v1.4.1 - Mitochondrial assembly
- **MAFFT** v7.487 - Multiple alignment
- **Gblocks** v0.91b - Alignment refinement
- **IQ-TREE** v2.3.3 - Phylogenetic inference

### Scripting

- **Python** 3.8+ (with BioPython)
- **R** 4.0+ (ggtree, ape, ggplot2, dplyr)
- **Bash** 4.0+

### Databases

- Kraken2 PlusPFP database (Sept 2024)
- NCBI GenBank (nucleotide)

---


## Citation

If you use this pipeline, please cite:

**Tools**:

- Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884-i890.
- Wood et al. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology* 20:257.
- Li (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094-3100.
- Jónsson et al. (2013). mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics* 29(13):1682-1684.
- Allio et al. (2020). MitoFinder: Efficient automated large‐scale extraction of mitogenomic data in target enrichment phylogenomics. *Molecular Ecology Resources* 20(4):892-905.
- Katoh & Standley (2013). MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution* 30(4):772-780.
- Nguyen et al. (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Molecular Biology and Evolution* 32(1):268-274.
- Yu et al. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees. *Methods in Ecology and Evolution* 8(1):28-36.

---

## Troubleshooting

**Low genome coverage?**

- Expected for historical DNA (0.5-2× typical)
- Focus on high-copy mitochondrial DNA
- Consider target enrichment for nuclear loci

**High contamination?**

- Review ECMSD and Kraken2 reports
- Adjust contaminant filtering thresholds
- Check for cross-contamination during extraction

**Poor tree resolution?**

- Increase gene sampling (add more loci)
- Use longer alignments or more taxa
- Consider partitioned models or codon-based models

**Missing sequences in concatenation?**

- Adjust species inclusion threshold (default: ≥50% genes)
- Check GenBank for additional sequences
- Verify gene name synonyms (COI=COX1, etc.)

---
