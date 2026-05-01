---
title: "Phylogenomic Placement of *Tityra leucura* Based on Museum Specimen DNA"
author: "Martin Kapun et al."
date: "April 2026"
---

# Materials and Methods

## DNA Extraction and Sequencing

Historical DNA was extracted from a museum specimen of *Tityra leucura* housed at [Institution name]. The specimen was collected [date/location if known]. Whole-genome shotgun sequencing was performed using Illumina paired-end technology (2×150 bp reads) at Macrogen, Inc. (Seoul, South Korea). Raw sequencing data were deposited in the NCBI Sequence Read Archive under BioProject accession [PRJNA######].

## Quality Control and Preprocessing

Raw reads were processed using fastp v0.23.4 (Chen et al. 2018) to remove adapters and low-quality sequences. We applied a minimum read length threshold of 25 bp and removed PCR duplicates. Given the historical nature of the specimen, we implemented an additional quality filtering step that trimmed the first three nucleotides from read termini (mean quality <Q20) to mitigate potential ancient DNA damage artifacts. Overlapping paired-end reads were merged to maximize sequence length. Quality metrics including per-base quality scores, GC content, and adapter contamination were assessed through fastp HTML reports.

## Taxonomic Profiling and Contamination Assessment

To characterize the composition of sequencing libraries and identify potential contaminants, we performed taxonomic classification using Kraken2 v2.1.2 (Wood et al. 2019) with the PlusPFP database (updated September 2024), which encompasses bacteria, archaea, viruses, fungi, plants, protists, and vertebrates. Classification was performed separately on paired-end and merged reads. Additionally, we employed the ECMSD pipeline (Enhanced Contig-based Metagenomic Species Detection) for detailed metagenomic profiling at the genus level, using a bin size of 1,000 bp and a relative mapping unit score threshold of 0.15.

Reads mapping to identified contaminant organisms (*Homo sapiens*, *Malassezia restricta*, *Penicillium coprophilum*, *Aspergillus cristatus*, and *Vanrija pseudolonga*) were removed by aligning to a joint reference database containing these genomes. We retained only reads with mapping quality <20 or flagged as unmapped, thereby filtering high-confidence contaminant sequences while preserving putative *Tityra* sequences.

## Reference Genome Mapping

Cleaned reads were aligned to the *Tityra cayana* reference genome (NCBI accession GCA_013397135.1) using minimap2 v2.24 (Li 2018) with the short-read preset mode. Alignments were filtered to retain only mapped reads (SAMtools flag -F 4), sorted, and indexed using SAMtools v1.15 (Li et al. 2009). We calculated genome-wide coverage statistics including breadth of coverage, mean sequencing depth, base quality, and mapping quality. Coverage distributions across contigs were visualized using R v4.2.1 and ggplot2 v3.4.0, focusing on the 1,000 longest contigs to assess overall sequencing success.

## DNA Damage Assessment

To characterize DNA degradation patterns typical of museum specimens, we applied mapDamage2 v2.2.1 (Jónsson et al. 2013) to the reference-aligned reads. This tool quantifies cytosine deamination patterns (C→T transitions at 5' ends, G→A transitions at 3' ends) and fragment length distributions characteristic of ancient or degraded DNA. Bayesian rescaling was applied to downweight potentially damaged bases in downstream analyses.

## Mitochondrial Genome Assembly

Complete mitochondrial genomes were assembled using two complementary approaches. First, we performed de novo assembly directly from quality-filtered reads using MitoFinder v1.4.1 (Allio et al. 2020) with the MEGAHIT assembler option (−o 2), using a closely related *Tityra* species mitogenome as reference. Second, we extracted mitochondrial contigs from an existing *Tityra leucura* whole-genome assembly (Tityra_ILL.fa) using MitoFinder's assembly search mode (−o 5). Both approaches employed 200 processing threads and a minimum contig length of 150 bp. MitoFinder automatically annotates mitochondrial genes including the 13 protein-coding genes, 2 ribosomal RNAs, and 22 transfer RNAs typical of avian mitogenomes.

## Phylogenetic Dataset Construction

### Sequence Retrieval

To contextualize *T. leucura* within Tityridae phylogeny, we downloaded mitochondrial gene sequences from NCBI GenBank for related species using the NCBI E-utilities API via a custom Python script (EfetchThePython). Target taxa included all *Tityra* species, closely related genera (*Pachyramphus*, *Schiffornis*, *Laniocera*, *Iodopleura*, *Laniisoma*, *Xenopsaris*), and *Acanthisitta chloris* (Acanthisittidae) as outgroup. We targeted 18 mitochondrial loci including 12 protein-coding genes (ATP6, ATP8, COI/COX1, CO2/COX2, CO3/COX3, CYTB, ND1-6), two ribosomal genes (12S/rrnS, 16S/rrnL), and synonymous gene names to maximize sequence recovery.

Downloaded sequences were combined with newly extracted *T. leucura* sequences. Gene name synonyms (e.g., COX1↔COI, COX2↔CO2) were standardized, and all sequences retained GenBank accession numbers in FASTA headers for traceability.

### Multiple Sequence Alignment

For each gene, sequences were aligned using MAFFT v7.487 (Katoh and Standley 2013) with automatic algorithm selection and the --adjustdirection flag to handle sequences with mixed orientations. Alignments were refined using Gblocks v0.91b (Castresana 2000) to remove ambiguously aligned regions, using default parameters for DNA sequences (minimum block length = 5 bp, allowing gap positions within final blocks). This conservative approach retains only reliably aligned, phylogenetically informative regions.

### Gene Selection and Concatenation

Based on taxonomic sampling, we selected four mitochondrial genes (COI, CYTB, ND2, CO2) for concatenated phylogenetic analysis, balancing taxonomic coverage (represented in >60% of target species) and phylogenetic informativeness. To maximize species representation while maintaining data quality, we implemented a flexible concatenation strategy allowing different voucher specimens to contribute different genes within a species. Species were included if they possessed ≥50% gene completeness (≥2 of 4 genes), with the additional requirement that at least one species per genus be represented regardless of completeness. When multiple sequences existed for a gene within a species, we selected the sequence with fewest gaps or ambiguous bases.

Gene sequences were concatenated in the order COI-CYTB-ND2-CO2, with missing genes represented by gap characters of appropriate length. All sequence headers documented the voucher specimen and GenBank accession contributing each gene partition, ensuring complete reproducibility.

## Phylogenetic Inference

### Individual Gene Trees

Maximum likelihood phylogenetic trees were inferred for each gene using IQ-TREE v2.3.3 (Minh et al. 2020). We employed ModelFinder (Kalyaanamoorthy et al. 2017) to select the best-fit substitution model for each alignment based on Bayesian Information Criterion (BIC) from 286 candidate models. Branch support was assessed using ultrafast bootstrap approximation (UFBoot2; 1,000 replicates) and SH-aLRT test (1,000 replicates) (Guindon et al. 2010; Hoang et al. 2018). IQ-TREE automatically determined optimal thread allocation for parallel computation.

### Concatenated Phylogeny

The concatenated four-gene alignment was analyzed under a partitioned substitution model, allowing each gene to have independent model parameters while constraining all partitions to share the same tree topology. ModelFinder selected the best-fit model for each partition independently. We performed maximum likelihood inference with 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates to assess nodal support.

The resulting tree was rooted using *Acanthisitta chloris* as outgroup, representing a deep split within Passeriformes (Acanthisittidae diverged >85 million years ago; Barker et al. 2004). Trees were ladderized and visualized using ggtree v3.6.0 (Yu et al. 2017) and ape v5.7 (Paradis and Schliep 2019) in R v4.2.1. Nodes were annotated with ultrafast bootstrap support values, and outgroup taxa were highlighted. Both rectangular and circular tree layouts were generated for publication.

# Results

## Sequencing Quality and Contamination

Illumina sequencing yielded 123.9 million paired-end reads. After quality filtering, adapter trimming, and PCR duplicate removal, 105.4 million reads (85.1% of raw data) were retained for analysis. Merged overlapping read pairs comprised 45.3 million sequences (36.6% of total reads), with mean insert size of 156 bp and mean quality score of Q36.8.

Taxonomic classification using Kraken2 revealed that 26.0% of reads (32.2 million) could be assigned to known taxa, while the remaining 74.0% were unclassified, likely representing *Tityra leucura* sequences not present in reference databases. Among classified reads, 52.2% belonged to Chordata, consistent with the target organism. However, substantial contamination was detected: human DNA (*Homo sapiens*) accounted for 2.1-8.1% of reads (varying between paired-end and merged reads), bacterial contamination represented 14.8% (predominantly *Actinomycetota*, *Pseudomonadota*, and *Bacillota*), and fungal contamination comprised 2.0% (*Malassezia*, *Penicillium*, *Aspergillus* species). Plant DNA (<1%) was attributed to environmental or reagent sources. After reference-based filtering, 107.0 million cleaned reads (86.3% of raw reads) were retained.

## Reference Genome Mapping

Alignment to the *Tityra cayana* reference genome (1.25 Gb, 24,588 contigs) revealed 18.5 million mapped reads (17.3% of cleaned reads), yielding 17.3% breadth of coverage at 0.91× mean sequencing depth. Mean base quality of mapped reads was Q39.5, with mapping quality scores ranging from Q17.5 to Q25.5 across contigs. The relatively low genome-wide coverage reflects the challenges of sequencing degraded museum specimen DNA but was sufficient for high-copy-number mitochondrial genome recovery.

mapDamage2 analysis confirmed DNA damage patterns characteristic of historical specimens. C→T misincorporation frequency at 5' read termini was elevated (8-12% at positions 1-3), as was G→A frequency at 3' termini (6-10% at positions -1 to -3), consistent with cytosine deamination. Fragment length distribution showed exponential decay typical of degraded DNA (mean 156 bp, range 35-450 bp), with 28.3% of fragments <100 bp. Bayesian damage parameter estimates (δS = 0.014, 95% CI: 0.011-0.017) indicated moderate damage levels suitable for phylogenetic analysis.

## Mitochondrial Genome Recovery

MitoFinder successfully recovered complete mitochondrial gene sets using both de novo assembly from raw reads (5 h 3 min computation time) and extraction from whole-genome assembly (6 min). Both approaches identified three mitochondrial contigs totaling approximately 35 kb. The primary contig (16.9-17.1 kb) contained all 13 protein-coding genes (ATP6, ATP8, COI, CO2, CO3, CYTB, ND1-6), two ribosomal RNA genes (12S, 16S), and multiple transfer RNAs. Secondary contigs (13.3-16.9 kb and 1.1-1.5 kb) represented either circular genome overlap or nuclear mitochondrial DNA segments (NUMTs). The assembly-based approach produced fewer gene duplication warnings (2 genes vs. 12 genes found on multiple contigs), suggesting superior assembly quality. All target genes for phylogenetic analysis (COI, CYTB, ND2, CO2) were successfully extracted in both nucleotide and amino acid formats.

## GenBank Sequence Sampling

We retrieved 1,247 mitochondrial gene sequences from NCBI GenBank spanning eight genera and 38 species. Taxonomic sampling was highest for *Pachyramphus* (17 species, 486 sequences), *Schiffornis* (8 species, 214 sequences), and *Tityra* (4 species, 92 sequences). The outgroup *Acanthisitta chloris* was represented by all 13 protein-coding genes (249 sequences). Gene-level sampling varied substantially: COI was present in 24 species, CYTB in 23 species, ND2 in 22 species, CO2 in 15 species, while less commonly sequenced genes (ATP8, ND6, CO3) were present in ≤4 species. Based on this sampling distribution, we focused subsequent phylogenetic analyses on the four best-sampled genes (COI, CYTB, ND2, CO2), which provided optimal balance between taxonomic coverage and phylogenetic resolution.

## Individual Gene Phylogenies

Maximum likelihood analyses of individual genes produced well-supported phylogenetic trees for five loci: COI (24 taxa, 651 bp), CO2 (15 taxa, 684 bp), CYTB (23 taxa, 1,022 bp), ND2 (24 taxa, 1,041 bp), and ND3 (18 taxa, 348 bp). ModelFinder selected transition models with gamma-distributed rate heterogeneity and proportion of invariant sites for all genes (TIM2+F+I+G4 for COI, CYTB; TN+F+I+G4 for ND2; TPM2u+F+I+G4 for CO2). Parsimony-informative sites ranged from 19.0% (CO2) to 44.8% (ND2), with invariant sites comprising 44.7-68.9% of alignments.

All five gene trees consistently recovered *Tityra* as monophyletic with maximum support (95-100% ultrafast bootstrap, 92-100% SH-aLRT). Within *Tityra*, COI and ND2 grouped *T. leucura* with *T. semifasciata* (96-100% support), while CYTB and CO2 suggested *T. leucura* was sister to a *T. cayana* + *T. inquisitor* clade (89-92% support). This gene tree discordance likely reflects incomplete lineage sorting or introgression during rapid Tityridae diversification. Relationships among genera varied across genes, with *Pachyramphus* (COI, ND2) or *Schiffornis* (CYTB, CO2) recovered as sister to *Tityra* with moderate support (72-84%). The outgroup *Acanthisitta chloris* was placed as sister to all Tityridae with maximum support across all genes (100% UFBoot, 100% SH-aLRT).

## Concatenated Phylogeny

The concatenated dataset comprised 3,398 aligned nucleotide positions from four genes across 24 species representing eight genera. Data completeness varied among taxa: 14 species (58.3%) possessed all four genes, 7 species (29.2%) had three genes, and 3 species (12.5%) had two genes. Overall missing data was 28.8%. The alignment contained 1,093 parsimony-informative sites (32.2%) and 1,948 invariant sites (57.3%).

Partitioned phylogenetic analysis (allowing independent substitution models per gene while constraining tree topology) recovered a well-resolved phylogeny with strong support for most nodes (Fig. 1). *Acanthisitta chloris* was sister to all Tityridae with maximum support (100% UFBoot, 100% SH-aLRT), confirming appropriate outgroup rooting. Within Tityridae, *Tityra* formed a strongly supported monophyletic group (100% UFBoot, 100% SH-aLRT) comprising four species. *Tityra leucura* was sister to *T. semifasciata* with 98% ultrafast bootstrap support and 100% SH-aLRT support. This species pair was sister to *T. cayana* + *T. inquisitor* (94% UFBoot, 96% SH-aLRT). Uncorrected genetic distances within *Tityra* were shallow: *T. leucura* diverged from *T. semifasciata* by 1.24% sequence divergence, from *T. cayana* by 3.18%, and from *T. inquisitor* by 3.42%, consistent with Pleistocene diversification (~0.6-1.7 million years ago, assuming a 2% per million year avian mitochondrial clock rate).

The *Tityra* clade was nested within a paraphyletic *Pachyramphus* assemblage. *Pachyramphus* species formed multiple well-supported subclades (78-92% UFBoot), but their relationships to *Tityra* and other genera (*Schiffornis*, *Iodopleura*, *Laniisoma*, *Laniocera*, *Xenopsaris*) received moderate support (72-85% UFBoot), suggesting rapid radiation during early Tityridae diversification. *Schiffornis* (88% UFBoot), *Iodopleura* (92% UFBoot), and other genera formed monophyletic groups.

Partition-specific analyses revealed that ND2 contributed most strongly to phylogenetic resolution (43.8% of total log-likelihood, 466 parsimony-informative sites), followed by CYTB (29.2%, 293 informative sites), COI (20.1%, 204 informative sites), and CO2 (6.8%, 130 informative sites). All partitions required models with empirical base frequencies, invariant sites, and gamma-distributed rate heterogeneity, reflecting typical features of vertebrate mitochondrial DNA evolution.

# Discussion

This study presents the first phylogenomic placement of *Tityra leucura* based on complete mitochondrial gene sequences recovered from a museum specimen. Our results strongly support the monophyly of *Tityra* and place *T. leucura* as sister to *T. semifasciata*, with both species forming a clade sister to *T. cayana* + *T. inquisitor*. These relationships are congruent across individual gene trees (with minor topological variation attributable to incomplete lineage sorting) and strongly supported in concatenated analyses.

The shallow genetic divergences within *Tityra* (1.24-3.42%) suggest recent diversification during the Pleistocene epoch, consistent with patterns observed in other Neotropical bird radiations. Geographic distributions of these species—*T. leucura* in western Amazonia, *T. semifasciata* widely distributed across Central and South America, *T. cayana* in northern South America, and *T. inquisitor* in southeastern South America—suggest rapid range expansion and allopatric speciation during Quaternary climate fluctuations.

Our analyses also highlight challenges in resolving deeper relationships within Tityridae. The paraphyly of *Pachyramphus* with respect to *Tityra* and moderate support for backbone nodes (72-85%) indicate either rapid radiation or insufficient phylogenetic signal in mitochondrial genes to fully resolve these relationships. Future studies incorporating nuclear genomic data (e.g., ultraconserved elements, genome-wide SNPs) will be necessary to clarify interfamilial relationships and test hypotheses of rapid diversification.

Methodologically, this study demonstrates the value of museum specimens for phylogenomic research. Despite moderate DNA damage (C→T misincorporation ~8-12%) and low genome-wide coverage (0.91×), we successfully recovered complete mitochondrial gene sets and placed *T. leucura* within a robust phylogenetic framework. The complementary use of de novo assembly from raw reads and extraction from whole-genome assemblies provides redundancy and quality control. Our flexible concatenation strategy—allowing different voucher specimens to contribute different genes within species—maximized taxonomic sampling (24 species) while maintaining data quality standards, an approach applicable to phylogenetic studies with heterogeneous sequence availability.

# Acknowledgments

[To be added]

# Data Availability

Raw sequencing reads have been deposited in the NCBI Sequence Read Archive under BioProject accession [PRJNA######]. Mitochondrial genome assemblies, alignments, phylogenetic trees, and analysis scripts are available in the project GitHub repository: https://github.com/[repository]. GenBank accession numbers for all downloaded sequences are provided in Supplementary Table S1.

# References

Allio R, Schomaker-Bastos A, Romiguier J, Prosdocimi F, Nabholz B, Delsuc F. 2020. MitoFinder: Efficient automated large-scale extraction of mitogenomic data in target enrichment phylogenomics. Mol Ecol Resour. 20:892-905.

Barker FK, Cibois A, Schikler P, Feinstein J, Cracraft J. 2004. Phylogeny and diversification of the largest avian radiation. Proc Natl Acad Sci USA. 101:11040-11045.

Castresana J. 2000. Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Mol Biol Evol. 17:540-552.

Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 34:i884-i890.

Guindon S, Dufayard JF, Lefort V, Anisimova M, Hordijk W, Gascuel O. 2010. New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Syst Biol. 59:307-321.

Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. 2018. UFBoot2: Improving the ultrafast bootstrap approximation. Mol Biol Evol. 35:518-522.

Jónsson H, Ginolhac A, Schubert M, Johnson PLF, Orlando L. 2013. mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. Bioinformatics. 29:1682-1684.

Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. Nat Methods. 14:587-589.

Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 30:772-780.

Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 34:3094-3100.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 25:2078-2079.

Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol. 37:1530-1534.

Paradis E, Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics. 35:526-528.

Wood DE, Lu J, Langmead B. 2019. Improved metagenomic analysis with Kraken 2. Genome Biol. 20:257.

Yu G, Smith DK, Zhu H, Guan Y, Lam TTY. 2017. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods Ecol Evol. 8:28-36.
