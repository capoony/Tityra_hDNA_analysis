# *Tityra leucura* Mitochondrial Phylogenomics Pipeline

This pipeline provides an end-to-end workflow for analyzing historical DNA from *Tityra leucura* museum specimens, focusing on mitochondrial genome assembly and phylogenetic analysis. The pipeline processes raw sequencing data through quality control, contamination assessment, mitochondrial genome assembly, and multi-gene phylogenetic reconstruction.

## Overview

This repository contains an end-to-end historical DNA workflow centered on the master script [`shell/main.sh`](shell/main.sh).

## Quick Start

```bash
WD=/media/inter/mkapun/projects/Tityra
bash ${WD}/shell/main.sh
```

## Pipeline Steps

The master workflow in [`shell/main.sh`](shell/main.sh) runs the following numbered sections in order.

### 1. Copy raw data into the working tree

The raw FASTQ files are staged locally so the rest of the analysis can run reproducibly from the project workspace.

### 2. Trim reads with fastp

The paired-end reads are adapter-trimmed, quality-filtered, deduplicated, and merged into a cleaned read set.

### 3. Run ECMSD taxonomic profiling

The trimmed reads are screened for broad metagenomic composition before more targeted mapping and assembly steps.

### 4. Run Kraken2 taxonomic classification

Kraken2 is used to identify likely contaminants and the dominant biological signal in the historical DNA library.

### 5. Map reads to the reference genome and generate summary figures

The cleaned reads are aligned to the reference genome to estimate coverage, depth, and read mapping quality.

### 6. Map trimmed reads to contaminant references and retain unmapped reads

Known contaminant genomes are used to filter out exogenous signal, leaving a decontaminated read set for downstream work.

### 7. Run the UCE workflow

The UCE companion workflow is called from the master pipeline to run the Tityra UCE mapping and phylogeny steps.

### 8. Run AutDeNovo assembly and annotation

The cleaned read set is assembled de novo and annotated to produce a contig-level assembly base for later mitochondrial work.

### 9. Assemble the mitochondrial genome

The decontaminated reads are used to reconstruct the mitochondrial genome for gene-level analysis.

### 10. Estimate sequencing depth on the assembled mitochondrial genome

Reads are mapped back to the assembled mitochondrial contig to confirm the recovered mitogenome is supported by the expected depth.

### 11. Retrieve additional mitochondrial marker sequences from related Tityrinae taxa

Additional mitochondrial loci are downloaded from related taxa to extend the comparative phylogenetic context.

### 12. Combine mitochondrial gene sequences

The downloaded marker sequences are merged with the Tityra-derived mitochondrial loci into one gene-by-gene dataset.

### 13. Run species-level phylogenetic analysis on the combined mitochondrial genes

The combined mitochondrial gene set is analyzed to recover the initial species-level phylogenetic signal.

### 14. Concatenate genes by species for the assembly-based matrix

The assembled mitochondrial markers are concatenated into a species-level supermatrix for joint analysis.

### 15. Aggregate mitochondrial gene metadata for all analysis samples

The accession and taxon metadata associated with the concatenated mitochondrial genes are consolidated into a single summary table.

### 16. Run phylogenetic analysis on the concatenated gene matrix

The concatenated mitochondrial matrix is analyzed to infer the main species-level phylogenetic relationships.

### 17. Run phylogenetic analysis on ND2 gene sequences only using the 2 July 2026 alignment

The ND2 alignment is analyzed separately to provide a comparison with the concatenated mitochondrial result.

### 18. Download complete mitogenomes from GenBank/RefSeq

Published complete mitochondrial genomes are downloaded to provide a broader comparative backbone for whole-mitogenome work.

### 19. Reorient and circularize the Tityra mitochondrial genome

The assembled Tityra mitochondrial genome is reoriented into a circularized representation suitable for downstream plotting.

### 20. Recreate the circular mitogenome map overview

The circular mitochondrial genome is visualized through the OGDRAW-style overview workflow to summarize genome architecture.

### 21. Align the complete mitochondrial genomes

The complete mitogenome references are aligned to prepare a global comparison across taxa.

### 22. Visualize the complete mitogenome alignments

The aligned mitogenome data are converted into plot-ready summaries for inspection of structural consistency and divergence.

### 23. Run phylogenetic analysis on the complete mitogenomes

A final whole-mitogenome tree is inferred to compare the Tityra assembly with complete mitogenome references.

### 24. Run mitochondrial genome synteny analysis

The gene order of the mitogenomes is compared to check for major structural rearrangements or concordance with the expected gene order.

### 25. Plot the GBIF occurrence map

A GBIF-based occurrence map is generated to place the molecular results in a broader geographic context.

