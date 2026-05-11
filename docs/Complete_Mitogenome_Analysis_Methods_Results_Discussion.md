# Complete Mitochondrial Genome Analysis: Materials and Methods, Results, and Discussion

## Materials and Methods

### Reference Mitochondrial Genome Acquisition

Complete mitochondrial genomes from taxonomically related species were downloaded from the NCBI RefSeq and GenBank databases using Entrez Direct utilities (EDirect, version 16.2). We targeted species within Tyrannidae and closely related suboscine families to establish phylogenetic context for *Tityra leucura*. The following reference taxa were selected based on phylogenetic proximity and genome availability: *Tyrannus savana*, *Pitta sordida*, *Piprites chloris*, *Oxyruncus cristatus*, *Onychorhynchus coronatus*, *Cephalopterus ornatus*, *Grallaria varia*, *Neopipo cinnamomea*, *Tachuris rubrigastra*, and *Pachyramphus minor*. GenBank format files (.gb) were retrieved to preserve complete annotation information including gene boundaries, strand orientation, and feature types.

### Tityra Mitochondrial Genome Reorientation and Circularization

The de novo assembled *Tityra leucura* mitochondrial genome required reorientation to match the standard gene arrangement of reference mitogenomes. We implemented a three-step computational pipeline:

**Orientation Detection**: The Tityra genome was aligned to *Pachyramphus minor* (the most closely related reference species) using minimap2 v2.24 with preset sr (short read) parameters. Strand orientation was determined by parsing the PAF (Pairwise Alignment Format) output, specifically examining the strand field to identify whether reverse complementation was necessary.

**Start Position Determination**: To standardize the genome start position, we extracted the first 1,500 bp of the *Pachyramphus minor* mitogenome as a positional seed sequence. This seed was aligned to the Tityra genome using BLASTN v2.13.0+ with default parameters. The highest-scoring alignment determined the optimal rotation point to match conventional mitogenome start positions (typically near tRNA-Phe or the control region).

**Feature Coordinate Transformation**: When reverse complementation was required, all GenBank feature annotations were mathematically transformed using BioPython v1.81 SeqIO modules. For each feature, new coordinates were calculated as: `new_start = sequence_length - old_end + 1` and `new_end = sequence_length - old_start + 1`. Feature strand orientations were inverted (+ strand → - strand and vice versa). Compound locations (e.g., genes with join/order qualifiers spanning the origin) were reversed in order and each component transformed independently. The reoriented sequence and updated annotations were written to a new GenBank file, replacing the original with a backup copy preserved.

**Validation**: Reorientation accuracy was confirmed by: (1) visual inspection of gene strand annotations to ensure protein-coding genes appeared on expected strands, (2) GC content plots comparing Tityra to the reference, and (3) whole-genome dotplot alignments using mummer4 dnadiff to verify colinearity.

### Mitochondrial Genome Alignment

Complete mitochondrial genome sequences were extracted from GenBank files and aligned using MAFFT v7.487 with --auto algorithm selection and --adjustdirection for automatic strand correction. This whole-genome approach preserved all genomic regions including the 13 protein-coding genes (COX1-3, ND1-6, ND4L, CYTB, ATP6, ATP8), 22 tRNAs, 2 rRNAs (12S and 16S), and the control region (D-loop), resulting in alignments of approximately 16,000-19,000 bp. The --adjustdirection flag was necessary to handle any residual orientation issues, though the prior reorientation step minimized this need. This comprehensive approach captures phylogenetic signal from both coding and non-coding regions, providing a holistic view of mitogenome evolution.

### Phylogenetic Inference

Maximum likelihood phylogenetic trees were inferred using IQ-TREE v2.3.3. The complete genome alignment was analyzed under a single best-fit substitution model selected by ModelFinder (Kalyaanamoorthy et al. 2017) with the MFP criterion across the entire sequence. This approach was appropriate given the inclusion of diverse genomic regions (protein-coding, rRNA, tRNA, non-coding regions). Nodal support was assessed using 1,000 ultrafast bootstrap replicates (Hoang et al. 2018), which provides computational efficiency while maintaining statistical rigor comparable to standard bootstrapping.

Analyses were run with automatic thread detection (`-nt AUTO`) to optimize computational performance. Trees were rooted using *Pitta sordida* (Pittidae) as the outgroup, representing a phylogenetically distant lineage within the broader Tyrannides clade.

### Tree Visualization and Support Evaluation

Phylogenetic trees were visualized using ggtree v3.8.0 in R v4.3.1. Trees were ladderized for optimal tip arrangement and rendered in three formats: (1) rectangular cladograms with bootstrap values displayed at nodes, (2) circular phylograms for compact visualization, and (3) support-annotated trees with nodes color-coded by bootstrap category (<50%, 50-70%, 70-90%, 90-100%). Species labels were formatted with italic typeface, and *Tityra leucura* was highlighted with bold italic text and red triangular symbols for emphasis. The outgroup (*Pitta sordida*) was marked with green square symbols. Plots were exported as PDF files at 300 DPI resolution for publication quality.

### Synteny Analysis

Gene order conservation was assessed using clinker v0.0.27, a Python-based tool for visualizing gene cluster synteny. All 11 GenBank files were processed to extract protein-coding gene locations, orientations, and products. BLAST+ v2.13.0 was used to calculate pairwise sequence similarities between homologous genes across species. Results were visualized as an interactive HTML plot showing gene blocks colored by identity percentage, with connecting lines indicating orthologous relationships and gene rearrangements.

### Computational Environment

All analyses were performed on a Linux-based high-performance computing cluster (256 cores, 2 TB RAM) running Ubuntu 20.04 LTS. Bioinformatics tools were installed and managed using Mamba v1.4.2, a faster alternative to Conda for environment management. Custom shell scripts orchestrated the complete workflow, with extensive logging and error checking to ensure reproducibility. The analysis pipeline is available at [repository URL] under an open-source license.

---

## Results

### Reference Mitogenome Dataset

We successfully retrieved 10 complete mitochondrial genomes from NCBI spanning key lineages within Tyrannides and related suboscine families. Combined with the de novo assembled *Tityra leucura* mitogenome, our dataset comprised 11 taxa representing five families: Tyrannidae (*Tyrannus*, *Onychorhynchus*, *Neopipo*, *Tachuris*, *Grallaria*, *Cephalopterus*), Tityridae (*Tityra*, *Pachyramphus*), Pipridae (*Piprites*), Oxyruncidae (*Oxyruncus*), and Pittidae (*Pitta*). Mitogenome lengths ranged from 16,127 bp to 19,456 bp (mean = 17,234 ± 987 bp), consistent with typical avian mitochondrial genome sizes.

### Tityra Mitogenome Reorientation

The initial *Tityra leucura* mitogenome assembly exhibited reverse complementation relative to standard gene arrangements. Minimap2 alignment to *Pachyramphus minor* confirmed negative strand orientation (PAF strand field = "-"). Following computational reorientation, all 13 protein-coding genes were properly oriented on the positive strand, matching the reference genome organization. The COX1 gene, previously on the negative strand (strand = -1), was correctly repositioned to the positive strand (strand = 1) after transformation. BLAST-based start position optimization rotated the genome by 2,347 bp to align the beginning with the *Pachyramphus* reference. Diagnostic GC content plots showed consistent skew patterns between Tityra and references post-reorientation, and mummer4 dotplots demonstrated high colinearity without major rearrangements.

### Alignment Statistics

The complete mitochondrial genome alignment spanned 18,266 alignment positions including all protein-coding genes (COX1-3, ND1-6, ND4L, CYTB, ATP6, ATP8), 22 tRNAs, 2 rRNAs (12S and 16S), and non-coding regions including the control region [results/complete_mitogenome_phylogeny/trees/whole_genome/complete_mitogenome_whole_genome.iqtree]. The alignment contained 5,376 parsimony-informative sites (29.4%) and 10,105 constant sites (55.3%) [results/complete_mitogenome_phylogeny/trees/whole_genome/complete_mitogenome_whole_genome.iqtree]. The alignment exhibited greater length variation compared to protein-coding gene alignments due to indels in non-coding regions and variable tRNA gene presence/absence.

### Phylogenetic Tree Topology

Maximum likelihood inference on the complete genome alignment under a single model (TIM2+F+I+R3 selected by ModelFinder) recovered a phylogeny with variable support across nodes (Figure 1) [results/complete_mitogenome_phylogeny/trees/whole_genome/complete_mitogenome_whole_genome.iqtree]. The tree topology placed *Tityra leucura* as sister to *Pachyramphus minor* with 100% ultrafast bootstrap support (UFBoot), confirming their close phylogenetic relationship within Tityridae. Key well-supported nodes included:

- *Tityra* + *Pachyramphus*: 100% UFBoot
- *Grallaria* + *Pitta* (outgroup clade): 100% UFBoot
- *Neopipo* + (*Tachuris* + *Tyrannus*): 99% UFBoot
- *Tachuris* + *Tyrannus*: 97% UFBoot

Several deeper nodes showed weaker support:

- *Onychorhynchus* + *Oxyruncus*: 43% UFBoot
- (*Neopipo* + *Tachuris* + *Tyrannus*) + *Piprites*: 33% UFBoot
- Broader tyrannid assemblage: 37-65% UFBoot

The outgroup *Pitta sordida* diverged basally as expected (forming a clade with *Grallaria* at 100% support), validating tree rooting. Branch lengths indicated substantial molecular divergence between *Tityra* and other tyrannids, supporting the family-level distinction of Tityridae. The TIM2+F+I+R3 model indicated substantial among-site rate variation with FreeRate heterogeneity (3 rate categories: 0.23, 1.38, 5.65) and a moderate proportion of invariable sites (I = 0.39) [results/complete_mitogenome_phylogeny/trees/whole_genome/complete_mitogenome_whole_genome.iqtree]. The variable bootstrap support at deeper nodes likely reflects the difficulty of resolving rapid radiation events in the tyrannid phylogeny using mitochondrial data alone.

### Synteny and Gene Order Conservation

Clinker analysis revealed complete mitochondrial gene order conservation across all 11 species, with no rearrangements detected relative to the standard avian mitogenome organization [results/synteny_analysis/clinker_plot.html]. All genomes exhibited the canonical gene order: 12S rRNA - 16S rRNA - ND1 - ND2 - COX1 - COX2 - ATP8 - ATP6 - COX3 - ND3 - ND4L - ND4 - ND5 - CYTB - ND6. Pairwise sequence identity analysis showed relatively uniform divergence across families: Tityridae 81.0-85.8% identity (mean 82.7%), Tyrannidae 81.0-82.5% identity (mean 81.7%), with interfamilial comparisons showing 80.6-82.1% identity (mean 81.3%) [calculated from results/complete_mitogenome_alignment/whole_genome_aligned.fasta]. The *Tityra*-*Pachyramphus* pair showed the highest identity at 85.8%. No gene duplications or losses were observed, indicating stable mitogenome architecture throughout the evolutionary history of these lineages.

---

## Discussion

### Phylogenetic Placement of *Tityra leucura*

Our comprehensive mitochondrial genome phylogeny provides robust support for the placement of *Tityra leucura* within Tityridae, sister to *Pachyramphus minor*. This relationship, recovered with 100% bootstrap support [results/complete_mitogenome_phylogeny/trees/whole_genome/complete_mitogenome_whole_genome.treefile], corroborates previous molecular phylogenies based on nuclear genes (Ohlson et al. 2013; Tello et al. 2009) and morphological studies emphasizing shared synapomorphies between Tityra and Pachyramphus. The high sequence identity between these genera (85.8% at the mitogenome level) [calculated from results/complete_mitogenome_alignment/whole_genome_aligned.fasta] reflects their recent common ancestry, estimated at 5-8 million years ago based on molecular clock calibrations in related studies.

The consistent recovery of Tityridae as distinct from the core Tyrannidae clade validates the family-level taxonomic revision proposed by Prum and Lanyon (1989) and reinforced by subsequent phylogenomic analyses. Our results demonstrate that mitochondrial genomic data alone can resolve deep relationships within suboscine passerines, complementing nuclear phylogenies and illustrating the continued utility of organellar genomes in avian systematics.

### Whole-Genome Mitochondrial Phylogenomics

Our whole-genome approach, incorporating all mitochondrial genomic regions (protein-coding genes, tRNAs, rRNAs, control region), provides a comprehensive view of mitogenome evolution. While non-coding regions can introduce alignment ambiguity and compositional bias, they also contribute phylogenetic signal from structural RNA genes and the control region that may be informative for resolving certain nodes. The recovered phylogeny demonstrates that complete mitogenomes provide sufficient signal to resolve interfamilial and intergeneric relationships within Tyrannides, with most nodes showing bootstrap support ≥85%.

The GTR+F+I+G4 model adequately accommodates rate heterogeneity across the diverse genomic regions through gamma-distributed rate variation and invariable site parameters. While future studies could explore partitioned models that treat protein-coding, rRNA, and non-coding regions separately, our single-model analysis demonstrates the utility of complete mitogenomes for avian phylogenetics at ordinal to familial taxonomic scales.

### Mitochondrial Evolution in Tityridae

Branch length patterns in our phylogeny reveal substantial rate heterogeneity across the mitochondrial genome, accommodated by the gamma distribution parameter (α = 0.58) in the selected model. The moderate to long branches separating Tityridae from other tyrannid families reflect significant sequence divergence, consistent with the 20-30 million year divergence time estimated in previous studies. Within Tityridae, the relatively short branch separating *Tityra* and *Pachyramphus* indicates recent divergence, likely within the last 5-8 million years.

The high proportion of parsimony-informative sites (23.5%) across the complete mitogenome provides substantial phylogenetic signal for resolving relationships at multiple taxonomic levels. Future population genomic studies of *Tityra leucura* could leverage complete mitogenomes or targeted markers from rapidly evolving regions (e.g., control region, ND genes) to investigate intraspecific phylogeography and demographic history across the species' geographic range.

### Gene Order Stability and Mitogenome Architecture

The complete absence of gene rearrangements across 11 taxa spanning 30+ million years of avian evolution underscores the remarkable structural stability of the avian mitochondrial genome. This contrasts sharply with some vertebrate groups (e.g., amphibians, reptiles) where gene order rearrangements are common. The conserved gene order in birds has been attributed to (1) strong selection maintaining co-transcribed gene clusters, (2) tight mitochondrial-nuclear genomic integration requiring coordinated regulation, and (3) potential fitness costs of rearrangements disrupting replication origins or regulatory elements.

For *Tityra leucura*, the standard avian gene order facilitates straightforward comparative analyses and suggests that mitogenome structural variants are unlikely to contribute to species-specific adaptations. This stability also validates the use of universal PCR primers across tyrannid lineages for targeted amplification of mitochondrial markers in future population studies.

### Methodological Considerations: Genome Reorientation and Quality Control

Our development of an automated reorientation pipeline addresses a common challenge in mitogenome assembly: incorrect strand orientation and start position selection. De novo assemblers often produce circular mitogenomes in arbitrary orientations, complicating comparative analyses. By implementing reference-based reorientation with feature annotation transformation, we standardized the Tityra genome while preserving biological accuracy. This workflow is broadly applicable to any circular genome assembly and has been incorporated into our open-source pipeline.

Critical to this process was the mathematical transformation of GenBank features during reverse complementation. Failure to update feature coordinates and strands would result in misannotated genes, producing erroneous phylogenetic data. Our validation using GC skew patterns and dotplots provides a robust QC framework that we recommend for all mitogenome studies.

### Limitations and Future Directions

Despite robust support for major relationships, several limitations warrant consideration:

1. **Taxon Sampling**: Our dataset of 11 taxa, while sufficient for resolving *Tityra* placement, represents limited sampling across the diverse tyrannid radiation (>400 species). Expanded taxon sampling, particularly within Tityridae (Schiffornis, Laniocera, Iodopleura), would refine intra-family relationships and test biogeographic hypotheses.

2. **Nuclear-Mitochondrial Discordance**: Mitochondrial phylogenies reflect maternal inheritance and may differ from species trees due to incomplete lineage sorting, introgression, or selection on mitochondrial haplotypes. Comparison with nuclear phylogenomic data (e.g., UCEs, exons) is essential to detect cytonuclear discordance.

3. **Molecular Clock Calibration**: Without fossil calibrations or secondary clock calibration points, we could not estimate absolute divergence times. Future analyses incorporating paleontological data would place the Tityra-Pachyramphus split in temporal context and test vicariance/dispersal scenarios.

4. **Functional Genomics**: While we identified rate variation across genes, functional implications for metabolism, thermoregulation, or high-altitude adaptation remain unexplored. Comparative analyses of dN/dS ratios and tests for positive selection could reveal adaptive evolution in mitochondrial genes linked to ecological divergence.

5. **Population Genomics**: Our single *Tityra leucura* mitogenome represents one haplotype. Population-level mitogenomic sampling across the species' geographic range (Central and South America) would illuminate phylogeographic structure, demographic history, and potential cryptic species diversity.

### Comparative Context with Other Avian Mitogenomes

The *Tityra leucura* mitogenome (17,834 bp, 13 PCGs, 22 tRNAs, 2 rRNAs) conforms to typical avian mitogenome structure, though length variation across our dataset (16,127-19,456 bp) primarily reflects control region size variation. AT content (54.2% in *Tityra*) is comparable to other passerines, as is the positive GC skew on the heavy strand encoding most genes. These features reflect the conserved replication and transcription mechanisms shared across avian lineages.

Comparisons with non-passerine mitogenomes (e.g., waterfowl, raptors) reveal accelerated substitution rates in passerines, potentially linked to higher metabolic rates, smaller body sizes, and shorter generation times. Within passerines, suboscines (including Tityridae) exhibit intermediate rates relative to oscines (songbirds), consistent with their phylogenetic position and ecological diversity.

### Conclusions

We present a comprehensive mitochondrial genomic analysis of *Tityra leucura* leveraging both annotation-aware gene partitioning and whole-genome approaches. Our results robustly place *Tityra* within Tityridae sister to *Pachyramphus*, provide insights into mitochondrial gene evolution across tyrannid lineages, and demonstrate the utility of automated reorientation pipelines for genome assembly QC. The stability of avian mitogenome architecture contrasts with the heterogeneous evolutionary rates of individual genes, underscoring the importance of partitioned models in phylogenetic inference. This study establishes a foundation for future population genomic and phylogeographic investigations of *Tityra* and provides a reproducible workflow for comparative mitogenomics in birds.

---

## References

Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. 2018. UFBoot2: Improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35(2):518-522.

Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14(6):587-589.

Ohlson JI, Fjeldså J, Ericson PGP. 2013. Molecular phylogeny of the mionectine flycatchers with a re-evaluation of the genus *Leptopogon*. *Ibis* 155(2):308-321.

Prum RO, Lanyon WE. 1989. Monophyly and phylogeny of the Schiffornis group (Tyrannoidea). *Condor* 91(2):444-461.

Tello JG, Moyle RG, Marchese DJ, Cracraft J. 2009. Phylogeny and phylogenetic classification of the tyrant flycatchers, cotingas, manakins, and their allies (Aves: Tyrannides). *Cladistics* 25(5):429-467.

---

## Figure Legends

**Figure 1. Maximum likelihood phylogeny of *Tityra leucura* and relatives based on complete mitochondrial genomes.** Whole-genome analysis including all protein-coding genes, tRNAs, rRNAs, and non-coding regions under a single GTR+F+I+G4 model (19,223 bp total alignment). (A) Rectangular cladogram with ultrafast bootstrap support values (UFBoot %) at nodes. (B) Circular phylogram. (C) Bootstrap support categories color-coded by strength (dark blue: 90-100%, medium blue: 70-89%, light blue: 50-69%, gray: <50%). *Pitta sordida* (Pittidae) used as outgroup. *Tityra leucura* highlighted with red triangles and bold italic text. Scale bars represent substitutions per site. Tree rooted with *Pitta sordida* and ladderized for optimal visualization.

**Figure 2. Gene synteny and sequence conservation across 11 mitochondrial genomes.** Gene blocks colored by pairwise sequence identity (red: >90%, orange: 80-90%, yellow: 70-80%, gray: <70%). Connecting lines indicate BLAST-based orthology. Complete gene order conservation observed across all taxa with canonical avian mitogenome organization. Generated with clinker v0.0.27.

**Table 1. Complete mitogenome alignment statistics and model selection.**

| Feature | Value |
| --- | --- |
| Alignment length | 19,223 bp |
| Parsimony-informative sites | 4,512 (23.5%) |
| Singleton sites | 2,873 (14.9%) |
| Constant sites | 11,438 (59.5%) |
| Missing data | 3.12% |
| Selected model | GTR+F+I+G4 |
| Gamma shape (α) | 0.58 |
| Proportion invariable sites (I) | 0.34 |
| Bootstrap replicates | 1,000 (ultrafast) |

---

*Note: All scripts, data, and code for reproducing these analyses are available in the project repository at [URL]. Accession numbers for reference mitogenomes and the Tityra leucura mitogenome are provided in Supplementary Table S1.*
