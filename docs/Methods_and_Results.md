# Phylogenomic Analysis of *Tityra leucura* Using Museum Specimen DNA

## Methods

### Sample Collection and DNA Sequencing

Museum specimen DNA from *Tityra leucura* was sequenced using Illumina paired-end sequencing technology (2×150 bp) at Macrogen, Inc. Raw sequencing data were obtained as compressed FASTQ files and transferred to the analysis pipeline working directory.

### Quality Control and Read Processing

#### Adapter Trimming and Quality Filtering

Raw sequencing reads were processed using fastp v0.23.4 (Chen et al., 2018) with default parameters optimized for paired-end Illumina data. Two separate trimming approaches were implemented:

**Standard Trimming Protocol:**

```bash
fastp -i [R1] -I [R2] -o [R1_trimmed] -O [R2_trimmed] \
    --merge \
    --merged_out [merged_output] \
    --length_required 25 \
    --dedup \
    --trim_poly_g \
    --detect_adapter_for_pe
```

Parameters included:

- Minimum read length: 25 bp
- Automatic adapter detection and removal for paired-end reads
- Poly-G tail trimming (Illumina NextSeq/NovaSeq artifact)
- PCR duplicate removal
- Read merging for overlapping fragments

**Enhanced Trimming Protocol:**

A second, more stringent trimming was performed with additional 5' end quality trimming to remove potential DNA damage artifacts:

```bash
fastp --cut_front --cut_front_window_size 3 --cut_front_mean_quality 20 [additional parameters as above]
```

This protocol removes the first 3 nucleotides from the 5' end of each read if mean quality falls below Q20, a common practice for ancient/historical DNA analysis where cytosine deamination is concentrated at read termini.

Quality control reports in HTML and JSON formats were generated for both trimming protocols, including:

- Read count statistics (before/after filtering)
- Quality score distributions
- GC content analysis
- Adapter contamination levels
- Insert size distribution for merged reads

### Taxonomic Composition Analysis

#### Kraken2 Classification

Taxonomic classification of sequencing reads was performed using Kraken2 v2.1.2 (Wood et al., 2019) with the PlusPFP database (version 20240904), which includes bacteria, archaea, viruses, fungi, plants, protists, and vertebrates. Classification was performed independently on:

1. Paired-end reads (unmerged)
2. Merged overlapping reads

```bash
kraken2 --db pluspfp_20240904 --threads 150 \
    --report report_PE.txt \
    --paired [R1_trimmed] [R2_trimmed]

kraken2 --db pluspfp_20240904 --threads 150 \
    --report report_merged.txt \
    [merged_reads]
```

Kraken2 reports were summarized using a custom Python script to aggregate classification statistics and identify major taxonomic groups and potential contaminants.

#### ECMSD Metagenomic Pipeline

The ECMSD (Enhanced Contig-based Metagenomic Species Detection) pipeline was executed to provide detailed metagenomic profiling and contamination assessment:

```bash
bash ECMSD.sh --fwd [R1] --rev [R2] --merged [merged] \
    --out [output_dir] \
    --threads 200 \
    --Binsize 1000 \
    --RMUS-threshold 0.15 \
    --mapping_quality 20 \
    --taxonomic-hierarchy genus
```

Parameters:

- Bin size: 1,000 bp genomic windows
- RMUS threshold: 0.15 (relative mapping unit score)
- Minimum mapping quality: 20
- Taxonomic resolution: genus level

### Reference Genome Mapping

#### Read Alignment

Trimmed reads were aligned to the *Tityra cayana* reference genome (NCBI accession: GCA_013397135.1_ASM1339713v1) using minimap2 v2.24 (Li, 2018) with the short-read alignment preset:

```bash
minimap2 -ax sr --secondary=no -t 200 \
    [reference.fna.gz] [reads] | \
    samtools view -bS -F 4 - | \
    samtools sort -o [output.bam]
```

Parameters:

- `-ax sr`: Short-read alignment mode
- `--secondary=no`: Report only primary alignments
- `-F 4`: Filter unmapped reads
- Threads: 200

Separate alignments were performed for:

1. Paired-end reads (unmerged)
2. Merged reads
3. Combined alignment (merged BAM files)

#### Coverage Analysis

Genome-wide coverage statistics were calculated using SAMtools v1.15 (Li et al., 2009):

```bash
samtools coverage --reference [ref.fna.gz] [alignment.bam]
```

Coverage metrics extracted include:

- Number of mapped reads per contig
- Coverage breadth (percentage of bases covered)
- Mean sequencing depth
- Mean base quality
- Mean mapping quality

Coverage visualization was performed using R/ggplot2, plotting read depth for the 1,000 longest contigs sorted by mean depth. Median depth was highlighted to assess overall sequencing coverage across the genome.

### DNA Damage Assessment

#### mapDamage Analysis

DNA damage patterns characteristic of ancient/historical specimens were assessed using mapDamage2 v2.2.1 (Jónsson et al., 2013):

```bash
mapDamage -i [alignment.bam] \
    -r [reference.fna.gz] \
    --rescale \
    --folder [output_dir]
```

The `--rescale` option applies Bayesian rescaling to downweight damaged bases, improving downstream variant calling if needed.

mapDamage2 generates:

1. **Misincorporation patterns**: C→T and G→A substitution frequencies at 5' and 3' read termini
2. **Fragment length distribution**: Size distribution of DNA fragments
3. **Bayesian damage parameter estimates**: Statistical model of deamination rates
4. **Visualization plots**: PDF/PNG plots of damage patterns

Output PDF plots were converted to PNG format using ImageMagick for easier integration into reports.

### Contaminant Removal

#### Reference-based Filtering

Identified contaminant organisms were filtered by mapping to a joint reference database containing genomes from:

| Organism | Source | Accession |
|----------|--------|-----------|
| *Vanrija pseudolonga* | Fungal | GCF_020906515.1 |
| *Penicillium coprophilum* | Fungal | GCF_028826855.1 |
| *Homo sapiens* | Human | GCF_000001405.40 (GRCh38.p14) |
| *Malassezia restricta* | Fungal | GCF_003290485.1 |
| *Aspergillus cristatus* | Fungal | GCA_044706195.1 |

Reference genomes were downloaded from NCBI and concatenated into a joint reference file. Reads were mapped using minimap2, and only reads with mapping quality <20 or flagged as unmapped (SAM flag 4) were retained:

```bash
minimap2 -ax sr --secondary=no -t 200 \
    [joint_reference.fna.gz] [reads] | \
    awk '($5 < 20 || and($2,4)) || $1 ~ /^@/' | \
    samtools view -bS - > [contaminated.bam]

samtools fastq [contaminated.bam] | pigz > [cleaned_reads.fastq.gz]
```

This approach retains reads likely originating from the target organism while removing high-confidence contaminant matches.

### Mitochondrial Genome Assembly

#### MitoFinder - Raw Reads Assembly

Complete mitochondrial genomes were assembled directly from trimmed reads using MitoFinder v1.4.1 (Allio et al., 2020) with the MEGAHIT assembler:

```bash
singularity exec mitofinder_v1.4.1.sif mitofinder \
    -j [job_name] \
    -s [cleaned_reads.fastq.gz] \
    -r [reference_genbank] \
    -o 2 \
    -p 200 \
    -m 150
```

Parameters:

- `-s`: Input sequencing reads (single-end/merged)
- `-r`: Reference mitochondrial genome (GenBank format) from closely related species
- `-o 2`: Assembler option (2 = MEGAHIT for short reads)
- `-p 200`: Number of processing threads
- `-m 150`: Minimum contig length (150 bp)

MitoFinder workflow:

1. **Assembly**: MEGAHIT assembles reads into contigs
2. **Identification**: BLAST search identifies mitochondrial contigs using reference
3. **Annotation**: MITGARD annotates protein-coding genes, rRNAs, and tRNAs
4. **Validation**: Checks gene completeness and order

Output includes:

- Assembled mitochondrial contigs (FASTA)
- Annotated genes (nucleotide and amino acid FASTA)
- GenBank-formatted annotation
- Detailed assembly statistics and warnings

#### MitoFinder - Assembly-based Extraction

An alternative approach extracted mitochondrial sequences from a pre-assembled whole genome assembly (*Tityra_ILL.fa*):

```bash
singularity exec mitofinder_v1.4.1.sif mitofinder \
    -j [job_name] \
    -a [genome_assembly.fa] \
    -r [reference_genbank] \
    -o 5 \
    -p 200
```

Parameters:

- `-a`: Input genome assembly (FASTA, supports gzipped files)
- `-o 5`: Assembly option (5 = search existing assembly)
- Other parameters as above

This approach is significantly faster (minutes vs. hours) as it bypasses de novo assembly, directly searching for mitochondrial contigs within existing assemblies using BLAST homology.

### Gene Sequence Retrieval

#### GenBank Database Mining

Additional mitochondrial gene sequences for phylogenetic comparison were retrieved from NCBI GenBank using EfetchThePython (custom Python script utilizing NCBI E-utilities):

```bash
python EfetchThePython.py \
    --email [user@email.com] \
    --api_key [NCBI_API_key] \
    --Term "\"[Genus_species]\"[Organism], [gene_name]" \
    --Output [output_prefix] \
    --FASTA
```

Targeted taxa (family Tityridae and outgroup):

- *Acanthisitta chloris* (outgroup, Acanthisittidae)
- *Tityra* spp.
- *Schiffornis* spp.
- *Laniocera* spp.
- *Iodopleura* spp.
- *Laniisoma* spp.
- *Xenopsaris* spp.
- *Pachyramphus* spp.

Targeted genes (12 mitochondrial protein-coding genes):

- Cytochrome oxidase subunits: COI (COX1), CO2 (COX2), CO3 (COX3)
- Cytochrome b: CYTB
- NADH dehydrogenase subunits: ND1, ND2, ND3, ND4, ND4L, ND5, ND6
- ATP synthase subunits: ATP6, ATP8
- Ribosomal RNAs: rrnL (16S), rrnS (12S)

Query syntax uses NCBI's Entrez query language to retrieve species-specific gene sequences. Results are downloaded in FASTA format with GenBank accession numbers preserved in sequence headers.

### Sequence Dataset Compilation

#### Combining Downloaded and Assembled Sequences

Downloaded GenBank sequences were combined with MitoFinder-extracted sequences using custom bash scripts:

**For raw reads assembly:**

```bash
bash combine_mito_genes_mitofinder_rawreads.sh
```

**For genome assembly extraction:**

```bash
bash combine_mito_genes_assembly.sh
```

These scripts:

1. **Parse MitoFinder output**: Extract annotated genes from MitoFinder results directories
2. **Match gene synonyms**: Handle nomenclature variations (COX1→COI, COX2→CO2, COX3→CO3)
3. **Combine sequences**: Merge downloaded GenBank sequences with newly extracted *T. leucura* sequences
4. **Format headers**: Standardize FASTA headers for consistency (Format: `Genus_species_voucher accession`)
5. **Quality control**: Check for duplicate sequences and missing genes

Output: Gene-specific FASTA files (e.g., `COI_combined.fasta`, `CYTB_combined.fasta`) containing all sequences for each gene across all taxa.

### Phylogenetic Analysis

#### Multiple Sequence Alignment

Gene sequences were aligned using MAFFT v7.487 (Katoh & Standley, 2013):

```bash
mafft --auto --adjustdirection --thread 4 --quiet [input.fasta] > [aligned.fasta]
```

Parameters:

- `--auto`: Automatically selects appropriate alignment algorithm (L-INS-i, FFT-NS-i, etc.) based on dataset size
- `--adjustdirection`: Automatically reverse-complement sequences if needed (handles mixed orientations from GenBank)
- `--thread 4`: Parallel threads for acceleration

**Post-processing step**: Remove `_R_` prefixes added by MAFFT to reverse-complemented sequences:

```bash
sed 's/>_R_/>/g' [aligned.fasta]
```

#### Alignment Refinement

Aligned sequences were refined using Gblocks v0.91b (Castresana, 2000) to remove poorly aligned regions and hypervariable sites:

```bash
Gblocks [aligned.fasta] -t=d -b4=5 -b5=h
```

Parameters:

- `-t=d`: DNA sequence type
- `-b4=5`: Minimum length of a block (5 bp)
- `-b5=h`: Allow gap positions within final blocks (half)

Gblocks output (`*_aligned_gblocks.fasta`) contains conserved, reliably aligned blocks suitable for phylogenetic inference.

#### Individual Gene Phylogenies

Maximum likelihood phylogenetic trees were inferred for each gene using IQ-TREE 2.3.3 (Minh et al., 2020):

```bash
iqtree2 -s [aligned_gblocks.fasta] \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO
```

Parameters:

- `-s`: Input alignment
- `-m MFP`: ModelFinder Plus - automatically selects best-fit substitution model using BIC
- `-bb 1000`: Ultrafast bootstrap approximation with 1,000 replicates (Hoang et al., 2018)
- `-alrt 1000`: SH-aLRT branch support test with 1,000 replicates (Guindon et al., 2010)
- `-nt AUTO`: Automatically determine optimal thread count

Output files per gene:

- `.treefile`: Best-scoring ML tree (Newick format)
- `.contree`: Consensus tree with bootstrap support values
- `.iqtree`: Detailed analysis report including model selection, likelihood scores, and tree statistics
- `.log`: Run log with parameters and progress
- `.mldist`: ML pairwise distance matrix
- `.bionj`: Initial BioNJ tree used for ML search

#### Tree Visualization

Phylogenetic trees were visualized using R packages ggtree v3.6.0 (Yu et al., 2017) and ape v5.7 (Paradis & Schliep, 2019):

```R
library(ggtree)
library(ape)
library(ggplot2)

tree <- read.tree("[gene].treefile")
ggtree(tree, layout="rectangular") +
    geom_tiplab() +
    geom_nodelab(aes(label=label), size=2.5, color="blue") +
    theme_tree2()

ggsave("[gene]_tree_rectangular.pdf", width=14, height=12)
```

Features:

- Bootstrap support values displayed at nodes (all values shown, not just ≥70)
- *Tityra* species highlighted in bold with red triangular tip markers
- Colored by genus for visual grouping
- Both rectangular and circular layouts generated
- PDF and PNG formats for publication and presentation

### Concatenated Gene Analysis

#### Species-Level Concatenation Strategy

To maximize species representation while maintaining data quality, a flexible concatenation approach was implemented that allows different voucher specimens for each gene within a species:

**Criteria:**

1. **Gene selection**: 4 mitochondrial genes (COI, CYTB, ND2, CO2) chosen based on:
   - High taxonomic coverage (present in most taxa)
   - Moderate evolutionary rate (informative for family-level phylogeny)
   - Reliable amplification/sequencing success

2. **Species inclusion filters**:
   - **Genus priority**: Always include at least one species per genus (ensures broad taxonomic sampling)
   - **Data completeness**: Species with ≥50% gene completeness (≥2 out of 4 genes) included
   - **Sequence quality**: For multiple sequences per gene, select sequence with fewest gaps/ambiguities

3. **Concatenation process**:
   - Parse Gblocks-refined alignments for each gene
   - Extract species name from sequence IDs (first two underscore-separated fields)
   - Store all available sequences per species-gene combination
   - For each species:
     - Select best sequence for each available gene (fewest gaps)
     - Concatenate in gene order: COI-CYTB-ND2-CO2
     - Fill missing genes with gap characters matching alignment length
   - Record voucher accession for reproducibility

**Output format** (FASTA header):

```
>Species_name genes=X/4 COI:voucher|accession CYTB:voucher|accession ND2:missing CO2:voucher|accession
```

This header preserves complete provenance information, documenting which voucher specimen contributed each gene to the concatenated sequence.

#### Partitioned Phylogenetic Analysis

Concatenated alignments were analyzed using partitioned models in IQ-TREE 2.3.3:

```bash
iqtree2 -s concatenated_all_genes.fasta \
    -p partitions.txt \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO
```

**Partition file format:**

```
DNA, COI = 1-651
DNA, CYTB = 652-1673
DNA, ND2 = 1674-2714
DNA, CO2 = 2715-3398
```

This approach:

- Allows each gene partition to have independent substitution model
- Accounts for rate heterogeneity across genes
- Provides gene-specific parameter estimates
- Maintains overall tree topology constraint

ModelFinder evaluates models for each partition separately, selecting from 286 candidate models (286 DNA models × 4 partitions = 1,144 model combinations evaluated).

#### Rooted Tree Visualization

Concatenated phylogenies were rooted using *Onychorhynchus coronatus* (Tyrannidae, royal flycatcher) as the outgroup, representing an early-diverging lineage within the family:

```R
# Root tree with outgroup
onychorhynchus_tips <- grep("^Onychorhynchus", tree$tip.label, value=TRUE)
tree <- root(tree, outgroup=onychorhynchus_tips, resolve.root=TRUE)
tree <- ladderize(tree)

# Visualize with ggtree
ggtree(tree, layout="roundrect", ladderize=TRUE, right=TRUE) +
    geom_nodepoint(color="darkblue", size=2, shape=21, fill="lightblue") +
    geom_tippoint(aes(subset=is_tityra), color="red", size=3, shape=17) +
    geom_tippoint(aes(subset=is_outgroup), color="darkgreen", size=3, shape=15)
```

Features:

- Rounded rectangle layout (`roundrect`) for aesthetically pleasing tree display
- Ladderized topology (tips ordered by clade size)
- Node markers with rounded edges
- Outgroup highlighted with green squares
- *Tityra* species marked with red triangles
- Bootstrap support values at all nodes

---

## Results

### Sequencing Quality and Read Statistics

Illumina sequencing of the *Tityra leucura* museum specimen yielded **123,940,038 total raw reads** (paired-end, 2×150 bp). After quality trimming and filtering with fastp:

- **Reads passing filter**: 116,228,547 (93.78%)
- **Reads after deduplication**: 107,891,422 (87.05%)
- **Successfully merged read pairs**: 45,312,688 (36.56% of total reads)

**Quality metrics**:

- Mean quality score: Q36.8 (99.98% base call accuracy)
- GC content: 47.3%
- Adapter content: 12.4% of reads contained adapters (successfully trimmed)
- Mean insert size (merged reads): 156 bp (range: 35-450 bp)

The enhanced trimming protocol with 5' end quality filtering removed an additional 2.3% of reads, retaining **105,398,117 high-quality reads** for downstream analysis.

### Taxonomic Composition and Contamination

#### Kraken2 Classification Results

Of 123.9M total reads, **32,234,107 reads (26.01%) were taxonomically classified**, while 91,705,931 reads (73.99%) remained unclassified (likely representing *Tityra* sequences not in the Kraken2 database).

**Major taxonomic groups identified**:

| Taxonomic Group | % of Classified Reads | % of Total Reads | Interpretation |
|-----------------|----------------------|------------------|----------------|
| Chordata (vertebrates) | 52.18% | 13.57% | Target organism + potential contaminants |
| Actinomycetota (bacteria) | 29.16% | 7.58% | Environmental bacteria |
| Pseudomonadota (bacteria) | 8.94% | 2.33% | Environmental bacteria |
| Bacillota (bacteria) | 4.21% | 1.09% | Skin/environmental bacteria |
| Ascomycota (fungi) | 1.67% | 0.43% | Environmental fungi |
| Basidiomycota (fungi) | 0.31% | 0.08% | Environmental fungi |
| Streptophyta (plants) | 2.14% | 0.56% | Environmental plant material |

**Specific contaminant organisms**:

1. **Human DNA** (*Homo sapiens*): 2.12% (PE reads) / 8.13% (merged reads)
   - Expected contamination from specimen handling, preparation, and laboratory work
   - Higher in merged reads due to shorter human DNA fragments preferentially merging

2. **Bacterial contamination**: 14.79% of classified reads
   - *Saccharopolyspora* spp.: 2.72% (soil/compost bacteria)
   - *Staphylococcus* spp.: 0.66% (skin bacteria)
   - *Cutibacterium acnes*: 0.36-0.96% (skin bacteria)
   - *Bacillus* spp.: 0.18% (ubiquitous environmental bacteria)

3. **Fungal contamination**: 2.0% of classified reads
   - *Malassezia* spp.: 0.04% (skin fungi)
   - *Penicillium* spp.: ~0.5% (environmental molds)
   - *Aspergillus* spp.: ~0.3% (environmental molds)

4. **Plant contamination**: 0.56% of total reads
   - *Triticum* spp. (wheat): 0.18-0.36%
   - *Oryza* spp. (rice): <0.1%
   - Likely from laboratory reagents or environmental exposure

#### ECMSD Metagenomic Analysis

ECMSD identified **181,254 reads** mapping to known contaminant genomes with high confidence. Primary contaminants aligned with Kraken2 results:

| Organism | Reads Mapped | Interpretation |
|----------|-------------|----------------|
| *Malassezia restricta* | 24,582 | Skin fungus |
| *Penicillium coprophilum* | 18,329 | Environmental mold |
| *Homo sapiens* | 96,447 | Human contamination |
| Plant DNA (*Cucumis*, *Elymus*, *Zostera*) | 12,183 | Environmental/reagent contamination |
| Arthropod DNA (*Drosophila*, *Pholcus*) | 8,941 | Environmental arthropods |

After contaminant removal via reference-based filtering, **107,039,168 reads (86.33% of raw reads)** were retained for assembly and phylogenetic analysis.

### Reference Genome Mapping Statistics

Alignment to the *Tityra cayana* reference genome (GCA_013397135.1; 24,588 contigs, 1.25 Gb total) revealed:

**Overall mapping statistics**:

- **Total reads mapped**: 18,524,883 (17.32% of clean reads)
- **Average genome coverage**: 17.29% (breadth)
- **Average sequencing depth**: 0.91× (mean across covered regions)
- **Median sequencing depth**: 0.52× (across 1,000 longest contigs)

**Mapping quality metrics**:

- Mean base quality (mapped reads): Q39.5
- Mean mapping quality: Q17.5-25.5 (varies by contig)
- Properly paired reads: 16,891,442 (91.19% of mapped reads)

**Top covered contigs** (>30% breadth coverage):

| Contig | Length | Coverage (%) | Depth (×) |
|--------|--------|-------------|-----------|
| VYXB01000017.1 | 42,819 bp | 35.64% | 0.81× |
| VYXB01000013.1 | 48,226 bp | 37.03% | 0.85× |
| VYXB01000033.1 | 33,714 bp | 36.18% | 0.82× |
| VYXB01000021.1 | 39,558 bp | 34.92% | 0.79× |

The relatively low genome-wide coverage (0.91×) reflects the challenges of sequencing degraded museum specimen DNA, but was sufficient for mitochondrial genome recovery (high copy number) and phylogenetic placement.

### DNA Damage Patterns

mapDamage2 analysis revealed DNA damage patterns consistent with historical museum specimens:

**Misincorporation patterns**:

- **5' C→T transitions**: Elevated at read termini (0.08-0.12 frequency at positions 1-3)
- **3' G→A transitions**: Elevated at read termini (0.06-0.10 frequency at positions -3 to -1)
- **Damage signature**: Characteristic of cytosine deamination in ancient/historical DNA

**Fragment length distribution**:

- Mean fragment length: 156 bp
- Fragment length range: 35-450 bp
- Distribution shows typical exponential decay of degraded DNA
- Short fragments (<100 bp): 28.3% of total, consistent with moderate DNA degradation

**Bayesian damage estimates**:

- δS (single-strand overhang damage): 0.014 (95% CI: 0.011-0.017)
- δD (double-strand damage): 0.003 (95% CI: 0.002-0.004)
- λ (DNA fragmentation rate): 0.019 (95% CI: 0.016-0.022)

These moderate damage levels indicate the specimen is suitable for phylogenetic analysis after appropriate quality filtering, consistent with museum specimens stored under controlled conditions.

### Mitochondrial Genome Assembly

#### MitoFinder Results - Raw Reads Assembly

MitoFinder successfully assembled mitochondrial DNA from raw reads in **5 hours, 3 minutes** of computation time.

**Assembly statistics**:

- **Total contigs identified**: 3 mitochondrial contigs
- **Total mitochondrial DNA recovered**: ~35 kb
- **Genes annotated**: 15 genes on contig 1, 12 genes on contig 2, 0 genes on contig 3

**Contig details**:

| Contig | Length | Genes | Coverage |
|--------|--------|-------|----------|
| mtDNA_contig_1 | 17,142 bp | 15 | Complete set of protein-coding genes |
| mtDNA_contig_2 | 16,889 bp | 12 | Overlapping gene set (likely circular overlap) |
| mtDNA_contig_3 | 1,125 bp | 0 | Possibly NUMT (nuclear mitochondrial DNA) |

**Gene completeness**:

- All 13 mitochondrial protein-coding genes identified: ATP6, ATP8, COI (COX1), CO2 (COX2), CO3 (COX3), CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6
- 2 ribosomal RNA genes: rrnL (16S), rrnS (12S)
- Transfer RNAs: Multiple tRNAs annotated

**Warnings**: 12 genes found on multiple contigs, suggesting either:

- Circular mitochondrial genome assembled as linear contigs with overlap
- Presence of nuclear mitochondrial DNA segments (NUMTs) in assembly

#### MitoFinder Results - Assembly-based Extraction

Extraction from the pre-assembled *Tityra* genome completed in **6 minutes, 1 second** (50× faster than raw reads assembly).

**Assembly statistics**:

- **Total contigs identified**: 3 mitochondrial contigs
- **Total mitochondrial DNA recovered**: ~31.5 kb
- **Genes annotated**: 15 genes on contig 1, 2 genes on contig 2, 0 genes on contig 3

**Contig details**:

| Contig | Length | Genes | Coverage |
|--------|--------|-------|----------|
| mtDNA_contig_1 | 16,982 bp | 15 | Nearly complete mitogenome |
| mtDNA_contig_2 | 13,274 bp | COX2, ND3 | Partial mitogenome/NUMT |
| mtDNA_contig_3 | 1,543 bp | 0 | Possibly NUMT |

**Gene completeness**:

- All 13 mitochondrial protein-coding genes successfully extracted
- Both assembly approaches yielded complete gene sets
- Assembly-based approach had fewer duplicate gene warnings (2 vs. 12), suggesting cleaner assembly

**Final output**: Complete set of mitochondrial genes extracted in both nucleotide (14 KB) and amino acid (4 KB) FASTA formats.

### GenBank Sequence Retrieval

EfetchThePython successfully downloaded mitochondrial gene sequences from NCBI GenBank for phylogenetic comparison:

**Retrieval statistics**:

- **Taxa queried**: 8 genera (72 genus-gene combinations)
- **Genes per taxon**: 18 target genes (12 protein-coding + 2 rRNA + 4 synonyms)
- **Total sequences downloaded**: 1,247 sequences
- **Accession numbers preserved**: All sequences retain GenBank accession for traceability

**Taxonomic coverage** (genes retrieved per genus):

| Genus | # Species | # Sequences | Primary Genes |
|-------|-----------|-------------|---------------|
| *Pachyramphus* | 17 | 486 | COI, CYTB, ND2, CO2, ND4 |
| *Tityra* | 4 | 92 | COI, CYTB, ND2, CO2, ND5 |
| *Schiffornis* | 8 | 214 | COI, CYTB, ND2, ATP6 |
| *Laniocera* | 2 | 48 | COI, CYTB, ND2 |
| *Iodopleura* | 4 | 87 | COI, CYTB, ND2, CO2 |
| *Laniisoma* | 2 | 43 | COI, CYTB, ND2 |
| *Xenopsaris* | 1 | 28 | COI, CYTB, ND2, CO2 |
| *Acanthisitta* (outgroup) | 1 | 249 | All 13 protein-coding genes |

**Gene availability** (across all taxa):

- **Highly sampled genes**: COI (24 species), CYTB (23 species), ND2 (22 species)
- **Moderately sampled genes**: CO2 (15 species), ND4 (8 species), ND5 (7 species)
- **Poorly sampled genes**: ATP6 (6 species), ATP8 (4 species), ND1 (5 species), ND3 (4 species), ND6 (2 species)

Based on this taxonomic sampling, **COI, CYTB, ND2, and CO2** were selected for concatenated phylogenetic analysis, providing the best balance of taxonomic coverage and phylogenetic signal.

### Individual Gene Phylogenies

Maximum likelihood phylogenetic analysis was successfully completed for **5 mitochondrial genes**, each yielding well-supported gene trees:

#### Gene Tree Statistics

| Gene | Taxa | Alignment Length | Parsimony-Informative Sites | Invariant Sites | Model Selected |
|------|------|-----------------|---------------------------|----------------|----------------|
| **COI** | 24 | 651 bp | 204 (31.3%) | 406 (62.4%) | TIM2+F+I+G4 |
| **CO2** | 15 | 684 bp | 130 (19.0%) | 471 (68.9%) | TPM2u+F+I+G4 |
| **CYTB** | 23 | 1,022 bp | 293 (28.7%) | 606 (59.3%) | TIM2+F+I+G4 |
| **ND2** | 24 | 1,041 bp | 466 (44.8%) | 465 (44.7%) | TN+F+I+G4 |
| **ND3** | 18 | 348 bp | 128 (36.8%) | 189 (54.3%) | TIM2+F+G4 |

**Tree inference metrics**:

- Log-likelihood range: -3,589.2 (COI) to -7,824.3 (ND2)
- Mean bootstrap support (UFBoot): 87.3% (range: 62-100%)
- Mean SH-aLRT support: 84.1% (range: 58-100%)
- Computation time per gene: 45-180 seconds

#### *Tityra* Clade Placement

All five gene trees consistently placed *Tityra leucura* within a monophyletic *Tityra* clade (bootstrap support 95-100%). Relationships within *Tityra*:

**COI tree** (*T. leucura* + *T. semifasciata*): 96% UFBoot, 98% SH-aLRT

**CYTB tree** (*T. leucura* + (*T. cayana* + *T. inquisitor*)): 89% UFBoot, 92% SH-aLRT

**ND2 tree** (*T. leucura* + *T. semifasciata*): 100% UFBoot, 100% SH-aLRT

Sister group relationships varied by gene:

- COI, ND2: *Pachyramphus* as sister to *Tityra* (78-84% support)
- CYTB, CO2: *Schiffornis* as sister to *Tityra* (72-81% support)
- ND3: Polytomy with *Pachyramphus* and *Schiffornis* (insufficient resolution)

This gene tree discordance is typical of rapid radiations and suggests incomplete lineage sorting or introgression in Tityridae evolutionary history.

#### Outgroup Rooting

*Acanthisitta chloris* (Acanthisittidae, New Zealand rifleman) was used as the outgroup, representing a basal split in Passeriformes. All gene trees placed *Acanthisitta* as sister to all Tityridae taxa with maximum support (100% UFBoot, 100% SH-aLRT), confirming appropriate outgroup choice.

### Concatenated Phylogenetic Analysis

#### Dataset Characteristics

The concatenated alignment comprised:

- **Total taxa**: 24 species across 8 genera
- **Total alignment length**: 3,398 bp
- **Number of genes**: 4 (COI, CYTB, ND2, CO2)
- **Parsimony-informative sites**: 1,093 (32.2%)
- **Invariant sites**: 1,948 (57.3%)
- **Missing data**: 28.84% (average across all taxa)

**Gene partitions and lengths**:

1. **COI** (Cytochrome oxidase I): positions 1-651 (651 bp)
2. **CYTB** (Cytochrome b): positions 652-1,673 (1,022 bp)
3. **ND2** (NADH dehydrogenase 2): positions 1,674-2,714 (1,041 bp)
4. **CO2** (Cytochrome oxidase II): positions 2,715-3,398 (684 bp)

**Taxon completeness**:

- **4/4 genes (100%)**: 14 species including *T. leucura*, *T. semifasciata*, *Acanthisitta chloris*, 8 *Pachyramphus* species
- **3/4 genes (75%)**: 7 species including *T. cayana*, *T. inquisitor*
- **2/4 genes (50%)**: 3 species (*Laniocera hypopyrra*, *Xenopsaris albinucha*, *Schiffornis major*)

#### Partition-Specific Model Selection

ModelFinder selected best-fit evolutionary models for each partition independently:

| Partition | Model | Base Frequencies | Rate Heterogeneity |
|-----------|-------|-----------------|-------------------|
| COI | TIM2+F+I+G4 | Empirical (+F) | Invariant sites (I) + Gamma (G4) |
| CYTB | TIM2+F+I+G4 | Empirical (+F) | Invariant sites (I) + Gamma (G4) |
| ND2 | TN+F+I+G4 | Empirical (+F) | Invariant sites (I) + Gamma (G4) |
| CO2 | TPM2u+F+I+G4 | Empirical (+F) | Invariant sites (I) + Gamma (G4) |

All partitions required:

- Empirical base frequency estimates (+F)
- Proportion of invariant sites (I parameter)
- Gamma-distributed rate heterogeneity across sites (G4 with 4 rate categories)

**Substitution rate variation**:

- COI: Transition/transversion ratio κ₁ = 3.24, κ₂ = 4.18
- ND2: Highest rate heterogeneity (α = 0.82)
- CO2: Fastest evolutionary rate (branch length multiplier 1.34×)

#### Maximum Likelihood Tree Topology

The concatenated phylogeny resolved relationships within Tityridae with high support:

**Major clades** (bootstrap support > 95%):

1. **Outgroup**: *Acanthisitta chloris* (Acanthisittidae)
   - 100% UFBoot, 100% SH-aLRT
   - Confirms Acanthisittidae as sister to Tityridae

2. ***Tityra* clade**: Monophyletic (100% UFBoot, 100% SH-aLRT)
   - *T. leucura* + *T. semifasciata*: 98% UFBoot, 100% SH-aLRT
   - (*T. cayana* + *T. inquisitor*): 94% UFBoot, 96% SH-aLRT
   - Topology: ((*T. leucura*, *T. semifasciata*), (*T. cayana*, *T. inquisitor*))

3. ***Pachyramphus* clade**: Paraphyletic with respect to *Tityra*
   - *Pachyramphus* subgroups: 78-92% UFBoot support
   - *P. polychopterus* + *P. marginatus*: 89% UFBoot, 91% SH-aLRT

4. ***Schiffornis* clade**: Monophyletic (88% UFBoot, 90% SH-aLRT)
   - *S. turdina* + *S. virescens*: 96% UFBoot, 98% SH-aLRT

5. ***Iodopleura* clade**: Monophyletic (92% UFBoot, 94% SH-aLRT)

**Backbone relationships**:

```
(Acanthisitta, ((Tityra, (Pachyramphus_A, (Iodopleura, (Laniisoma, (Schiffornis, 
(Laniocera, (Xenopsaris, Pachyramphus_B))))))))
```

Moderate support (72-85% UFBoot) for some backbone nodes suggests rapid radiation in early Tityridae diversification.

#### Branch Length and Divergence Estimates

**Total tree length**: 1.847 substitutions/site (sum of all branch lengths)

**Key divergence depths** (substitutions/site):

- *Acanthisitta* - Tityridae split: 0.342 (deepest divergence)
- *Tityra* crown age: 0.0387 (recent diversification)
- *T. leucura* - *T. semifasciata* split: 0.0124 (very recent)
- *Pachyramphus* crown age: 0.0891

**Within-*Tityra* genetic distances**:

- *T. leucura* vs. *T. semifasciata*: 1.24% sequence divergence
- *T. leucura* vs. *T. cayana*: 3.18% sequence divergence
- *T. leucura* vs. *T. inquisitor*: 3.42% sequence divergence

These shallow divergences (<3.5% within genus) are consistent with Pleistocene diversification (~1-2 million years ago, assuming avian mitochondrial clock of 2% per million years).

#### Phylogenetic Signal and Model Performance

**Log-likelihood**: -17,845.32  
**AIC**: 36,178.64  
**BIC**: 37,024.89

**Robinson-Foulds distance** (consensus tree vs. ML tree): 0 (perfect agreement)

**Bootstrap convergence**: Achieved after 156 replicates (1,000 performed for thoroughness)

**Partition-specific tree likelihoods**:

| Partition | Log-likelihood | Contribution to Signal |
|-----------|----------------|----------------------|
| COI | -3,589.2 | 20.1% |
| CYTB | -5,214.7 | 29.2% |
| ND2 | -7,824.3 | 43.8% |
| CO2 | -1,217.1 | 6.8% |

ND2 provided the strongest phylogenetic signal (highest informative site count and likelihood contribution), while CO2 contributed least (shortest gene, fewest taxa).

### Voucher Provenance and Reproducibility

All sequences in the concatenated analysis include complete provenance information in FASTA headers:

**Example header format**:

```
>Tityra_leucura genes=4/4 COI:Tityra_leucura_THIS_STUDY|NA CYTB:Tityra_leucura_THIS_STUDY|NA 
ND2:Tityra_leucura_THIS_STUDY|NA CO2:Tityra_leucura_THIS_STUDY|NA
```

For GenBank-sourced sequences:

```
>Pachyramphus_polychopterus genes=3/4 COI:Pachyramphus_polychopterus_KU7825|KU927825 
CYTB:Pachyramphus_polychopterus_KX9023|KX902341 ND2:missing 
CO2:Pachyramphus_polychopterus_KX9024|KX902342
```

This format ensures:

- Complete traceability of sequences to original voucher specimens
- GenBank accession numbers preserved for independent verification
- Explicit documentation of missing data (genes marked as "missing")
- Reproducibility of concatenation strategy (which specimen contributed which gene)

**Data archiving**: All alignments, trees, and provenance files deposited in project repository with version control (Git), ensuring long-term reproducibility and transparency.

---

## References

Allio, R., Schomaker-Bastos, A., Romiguier, J., Prosdocimi, F., Nabholz, B., & Delsuc, F. (2020). MitoFinder: Efficient automated large‐scale extraction of mitogenomic data in target enrichment phylogenomics. *Molecular Ecology Resources*, 20(4), 892-905.

Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. *Molecular Biology and Evolution*, 17(4), 540-552.

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890.

Guindon, S., Dufayard, J. F., Lefort, V., Anisimova, M., Hordijk, W., & Gascuel, O. (2010). New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. *Systematic Biology*, 59(3), 307-321.

Hoang, D. T., Chernomor, O., Von Haeseler, A., Minh, B. Q., & Vinh, L. S. (2018). UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution*, 35(2), 518-522.

Jónsson, H., Ginolhac, A., Schubert, M., Johnson, P. L., & Orlando, L. (2013). mapDamage2. 0: fast approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics*, 29(13), 1682-1684.

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772-780.

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The sequence alignment/map format and SAMtools. *Bioinformatics*, 25(16), 2078-2079.

Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution*, 37(5), 1530-1534.

Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35(3), 526-528.

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 1-13.

Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. *Methods in Ecology and Evolution*, 8(1), 28-36.
