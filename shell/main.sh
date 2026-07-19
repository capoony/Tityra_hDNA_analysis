#!/bin/bash
set -euo pipefail
###############################################################################
# Tityra hDNA Analysis — Master Pipeline
# Run each section manually in order.
###############################################################################

###############################################################################
# Variables
###############################################################################

readonly WD="/media/inter/mkapun/projects/Tityra"
readonly DATA="${WD}/data"
readonly OUT="${WD}/results"
readonly REF_DIR="${DATA}/ref"
readonly REF_FASTA="GCA_013397135.1_ASM1339713v1_genomic.fna.gz"
readonly REF_URL="ftp:/ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/397/135/GCA_013397135.1_ASM1339713v1/${REF_FASTA}"

###############################################################################
# 1. Copy raw data into the working tree
# This step stages the input FASTQ files in a stable local directory so the rest
# of the analysis can be run reproducibly from the project workspace.
###############################################################################

mkdir -p "${DATA}/raw"
cp /media/inter/SeqData/raw/Macrogen/Illumina/Tityra_250430/data/*.gz "${DATA}/raw"

###############################################################################
# 2. Trim reads with fastp
# This step removes adapters, filters low-quality bases, merges paired reads, and
# produces a clean read set for downstream taxonomic and mapping analyses.
###############################################################################

mkdir -p "${DATA}/trimmed"
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

# Run fastp for paired-end trimming, merging, deduplication, and QC reporting
fastp \
    -i "${WD}/data/raw/Tityra_leucura_1.fastq.gz" \
    -I "${WD}/data/raw/Tityra_leucura_2.fastq.gz" \
    -o "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    -O "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" \
    --merge \
    --merged_out "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" \
    --length_required 25 \
    --dedup \
    --trim_poly_g \
    --html "${WD}/data/trimmed/Tityra_leucura.html" \
    --json "${WD}/data/trimmed/Tityra_leucura.json" \
    --detect_adapter_for_pe

mkdir -p "${WD}/data/trimmed2" # Create directory for trimmed reads if it doesn't exist
# Run fastp for paired-end trimming, merging, deduplication, and QC reporting and additionally trim away the first three bases of each read at the 5'end 
fastp \
    -i "${WD}/data/raw/Tityra_leucura_1.fastq.gz" \
    -I "${WD}/data/raw/Tityra_leucura_2.fastq.gz" \
    -o "${WD}/data/trimmed2/Tityra_leucura_1_trimmed.fastq.gz" \
    -O "${WD}/data/trimmed2/Tityra_leucura_2_trimmed.fastq.gz" \
    --merge \
    --merged_out "${WD}/data/trimmed2/Tityra_leucura_merged.fastq.gz" \
    --length_required 25 \
    --cut_front --cut_front_window_size 3 --cut_front_mean_quality 20 \
    --dedup \
    --trim_poly_g \
    --html "${WD}/data/trimmed2/Tityra_leucura.html" \
    --json "${WD}/data/trimmed2/Tityra_leucura.json" \
    --detect_adapter_for_pe

conda deactivate # Deactivate conda environment

###############################################################################
# 3. Run ECMSD taxonomic profiling
# This step screens the trimmed reads for broad metagenomic composition before
# more targeted mapping and assembly work begins.
###############################################################################

bash /media/inter/pipelines/ECMSD/shell/ECMSD.sh \
    --fwd "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    --rev "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" \
    --merged "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" \
    --out "${WD}/results/ECMSD" \
    --threads 200 \
    --cov-threshold 10 \
    --mapping_quality 1 \
    --taxonomic-hierarchy genus \
    --force

###############################################################################
# 4. Run Kraken2 taxonomic classification
# This step assigns the reads to taxa to identify likely contaminants and the
# dominant biological signal in the historical DNA library.
###############################################################################

mkdir -p "${OUT}/kraken2"
module load Assembly/kraken-2.1.2

# Run Kraken2 on paired-end reads for taxonomic classification
kraken2 \
    --db /media/scratch/kraken-2.1.2/db/pluspfp_20240904 \
    --threads 150 \
    --report "${WD}/results/kraken2/report_PE.txt" \
    --paired \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" \
    >/dev/null

# Run Kraken2 on merged reads
kraken2 \
    --db /media/scratch/kraken-2.1.2/db/pluspfp_20240904 \
    --threads 150 \
    --report "${WD}/results/kraken2/report_merged.txt" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" \
    >/dev/null

# Summarize Kraken2 results using custom Python script
python "${WD}/scripts/summarizeKraken.py" \
    --input "${WD}/results/kraken2/" \
    --output "${WD}/results/kraken2/kraken_summary.csv"

# Create combined visualization with Kraken2 and ECMSD results (3 panels)
Rscript "${WD}/scripts/visualize_combined_taxonomy.R" "${WD}"

###############################################################################
# 5. Map reads to the reference genome and generate summary figures
# This step quantifies how much of the trimmed library aligns to the closest
# reference genome and produces read-depth and coverage summaries for QC.
###############################################################################

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"
if [ ! -f "${REF_FASTA}" ]; then
    wget -q "${REF_URL}"
fi

mkdir -p "${WD}/results/minimap2" # Create minimap2 output directory
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs # Activate conda env

# Align paired-end reads to reference genome with minimap2, filter, and sort
minimap2 -ax sr --secondary=no -t 200 \
    "${REF_DIR}/${REF_FASTA}" \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/minimap2/Tityra_leucura_PE.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura_PE.bam" # Index BAM

# Align merged reads to reference genome
minimap2 -ax sr --secondary=no -t 200 \
    "${REF_DIR}/${REF_FASTA}" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/minimap2/Tityra_leucura_merged.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura_merged.bam"

# Merge paired-end and merged BAM files
samtools merge -f "${WD}/results/minimap2/Tityra_leucura.bam" \
    "${WD}/results/minimap2/Tityra_leucura_PE.bam" \
    "${WD}/results/minimap2/Tityra_leucura_merged.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura.bam"

# Deduplicate the merged BAM with samtools markdup
samtools sort -n -o "${WD}/results/minimap2/Tityra_leucura.namesorted.bam" \
    "${WD}/results/minimap2/Tityra_leucura.bam"
samtools fixmate -m "${WD}/results/minimap2/Tityra_leucura.namesorted.bam" \
    "${WD}/results/minimap2/Tityra_leucura.fixmate.bam"
samtools sort -o "${WD}/results/minimap2/Tityra_leucura.fixmate.sorted.bam" \
    "${WD}/results/minimap2/Tityra_leucura.fixmate.bam"
samtools markdup -r \
    "${WD}/results/minimap2/Tityra_leucura.fixmate.sorted.bam" \
    "${WD}/results/minimap2/Tityra_leucura.dedup.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura.dedup.bam"
rm "${WD}/results/minimap2/Tityra_leucura.namesorted.bam" \
   "${WD}/results/minimap2/Tityra_leucura.fixmate.bam" \
   "${WD}/results/minimap2/Tityra_leucura.fixmate.sorted.bam"

# Calculate coverage statistics on deduplicated BAM
samtools coverage --reference "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz" \
    "${WD}/results/minimap2/Tityra_leucura.dedup.bam" \
    >"${WD}/results/minimap2/Tityra_leucura.dedup.coverage.txt"

# Calculate coverage statistics on full (pre-dedup) BAM (used by endogenous pie chart)
samtools coverage --reference "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz" \
    "${WD}/results/minimap2/Tityra_leucura.bam" \
    >"${WD}/results/minimap2/Tityra_leucura.coverage.txt"

# ─────────────────────────────────────────────────────────────────────────────
# Endogenous DNA pie chart
# Count total basepairs in trimmed reads (PE R1 + R2 + merged) and basepairs
# mapped to the closest relative (Pachyramphus minor) to estimate endogenous DNA.
# ─────────────────────────────────────────────────────────────────────────────

# Sum basepairs in each trimmed FASTQ (every 2nd line out of 4 is the sequence)
BP_R1=$(zcat "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" | awk 'NR%4==2{sum+=length($0)}END{print sum}')
BP_R2=$(zcat "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" | awk 'NR%4==2{sum+=length($0)}END{print sum}')
BP_MERGED=$(zcat "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" | awk 'NR%4==2{sum+=length($0)}END{print sum}')
TOTAL_BP=$(( BP_R1 + BP_R2 + BP_MERGED ))

# Sum basepairs of primary mapped reads from the BAM (query sequence length)
MAPPED_BP=$(samtools view -F 4 -F 256 "${WD}/results/minimap2/Tityra_leucura.bam" | awk '{sum+=length($10)}END{print sum}')

# Generate pie chart with R
Rscript -e "
total    <- ${TOTAL_BP}
mapped   <- ${MAPPED_BP}
unmapped <- total - mapped
pct_endo <- round(mapped   / total * 100, 2)
pct_exo  <- round(unmapped / total * 100, 2)

library(ggplot2)
df <- data.frame(
    Category   = c('Endogenous (mapped)', 'Exogenous (unmapped)'),
    BP         = c(mapped, unmapped),
    Proportion = c(pct_endo, pct_exo)
)
df\$label <- paste0(df\$Category, '\n', format(df\$BP, big.mark=','), ' bp (', df\$Proportion, '%)')

ggplot(df, aes(x = '', y = BP, fill = Category)) +
    geom_col(width = 1, colour = 'white') +
    coord_polar(theta = 'y') +
    scale_fill_manual(values = c('Endogenous (mapped)' = '#2166AC',
                                  'Exogenous (unmapped)' = '#D73027')) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 4.5, colour = 'white', fontface = 'bold') +
    labs(title = 'Proportion of Endogenous DNA',
         subtitle = paste0('Total trimmed basepairs: ', format(total, big.mark=','),
                           '  |  Reference: Pachyramphus minor (VYXB01)'),  
         fill = NULL) +
    theme_void(base_size = 13) +
    theme(plot.title    = element_text(hjust = 0.5, face = 'bold', size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = 'none')
ggsave('${WD}/results/minimap2/Tityra_leucura_endogenous_piechart.png',
       width = 7, height = 7, dpi = 300)
ggsave('${WD}/results/minimap2/Tityra_leucura_endogenous_piechart.pdf',
       width = 7, height = 7, dpi = 300)
cat('Endogenous DNA pie chart saved.\n')
cat('Total bp:   ', total,   '\n')
cat('Mapped bp:  ', mapped,  ' (', pct_endo, '%)\n')
cat('Unmapped bp:', unmapped,' (', pct_exo,  '%)\n')
"

# Compound figure: histogram of sequencing depth + histogram of breadth of
# coverage + duplication level across the 1000 longest contigs (ranked by mean depth)
Rscript -e "
library(tidyverse)
library(patchwork)

cov      <- read.table('${WD}/results/minimap2/Tityra_leucura.coverage.txt',
                       header = TRUE, comment.char = '', check.names = FALSE)
cov_dedup <- read.table('${WD}/results/minimap2/Tityra_leucura.dedup.coverage.txt',
                        header = TRUE, comment.char = '', check.names = FALSE)
names(cov)[1]       <- 'rname'
names(cov_dedup)[1] <- 'rname'

# Per-contig duplication rate = 1 - dedup_reads / all_reads
dup <- inner_join(cov, cov_dedup, by = 'rname', suffix = c('', '.dedup')) %>%
    mutate(fold_dup = ifelse(numreads.dedup > 0,
                             numreads / numreads.dedup, 1))

# Select 1000 longest contigs (endpos = contig length in samtools coverage)
cov1000 <- cov_dedup %>%
    arrange(desc(endpos)) %>%
    slice_head(n = 1000)

dup1000 <- dup %>% filter(rname %in% cov1000\$rname)

med_depth <- median(cov1000\$meandepth)
med_cov   <- median(cov1000\$coverage)
med_dup   <- median(dup1000\$fold_dup)

# Panel A – histogram of mean sequencing depth
pA <- ggplot(cov1000, aes(x = meandepth)) +
    geom_histogram(bins = 50, fill = '#2166AC', colour = 'white', alpha = 0.85) +
    geom_vline(xintercept = med_depth, colour = '#D73027',
               linetype = 'dashed', linewidth = 0.8) +
    annotate('text', x = med_depth, y = Inf,
             label = paste0('Median: ', round(med_depth, 2), 'x'),
             colour = '#D73027', hjust = -0.1, vjust = 1.5, size = 3.5) +
    scale_x_log10() +
    labs(title = 'A  Sequencing Depth (deduplicated)',
         x = 'Mean depth (×, log scale)', y = 'Number of contigs') +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = 'bold'))

# Panel B – histogram of breadth of coverage (%)
pB <- ggplot(cov1000, aes(x = coverage)) +
    geom_histogram(bins = 50, fill = '#4DAC26', colour = 'white', alpha = 0.85) +
    geom_vline(xintercept = med_cov, colour = '#D73027',
               linetype = 'dashed', linewidth = 0.8) +
    annotate('text', x = med_cov, y = Inf,
             label = paste0('Median: ', round(med_cov, 2), '%'),
             colour = '#D73027', hjust = -0.1, vjust = 1.5, size = 3.5) +
    labs(title = 'B  Breadth of Coverage (deduplicated)',
         x = 'Coverage (%)', y = 'Number of contigs') +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = 'bold'))

# Panel C – histogram of fold duplication per contig
pC <- ggplot(dup1000, aes(x = fold_dup)) +
    geom_histogram(bins = 50, fill = '#B2571A', colour = 'white', alpha = 0.85) +
    geom_vline(xintercept = med_dup, colour = '#D73027',
               linetype = 'dashed', linewidth = 0.8) +
    annotate('text', x = med_dup, y = Inf,
             label = paste0('Median: ', round(med_dup, 2), '×'),
             colour = '#D73027', hjust = -0.1, vjust = 1.5, size = 3.5) +
    scale_x_log10() +
    labs(title = 'C  Fold Duplication',
         x = 'Fold duplication (×, log scale)', y = 'Number of contigs') +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = 'bold'))

compound <- pA + pB + pC +
    plot_annotation(
        title    = 'Sequencing depth, coverage, and fold duplication across the 1,000 longest contigs (deduplicated)',
        subtitle = 'Reference: Pachyramphus minor (VYXB01)',
        theme = theme(
            plot.title    = element_text(hjust = 0.5, face = 'bold', size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
    )

ggsave('${WD}/results/minimap2/Tityra_leucura.coverage_plot.png',
       plot = compound, width = 17, height = 5, dpi = 300)
ggsave('${WD}/results/minimap2/Tityra_leucura.coverage_plot.pdf',
       plot = compound, width = 17, height = 5, dpi = 300)
cat('Coverage compound figure saved.\n')
"

# Prepare reference genome for mapDamage (bgzip and index)
gunzip -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz"
bgzip -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna"
samtools faidx "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz"

conda deactivate
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/mapdamage2 # Activate mapDamage env

# Run mapDamage to assess DNA damage patterns
mapDamage -i "${WD}/results/minimap2/Tityra_leucura.bam" \
    -r "${REF_DIR}/${REF_FASTA}" \
    --rescale \
    --folder="${WD}/results/mapDamage/Tityra_leucura"

conda deactivate

# Convert mapDamage PDF plots to PNG images
for pdf in "${WD}/results/mapDamage/Tityra_leucura/"*.pdf; do
    png="${pdf%.pdf}.png"
    convert -density 300 "$pdf" -quality 90 "$png"
done

echo "Reference mapping and coverage summary completed."

###############################################################################
# 6. Map trimmed reads to contaminant references and retain unmapped reads
# This step removes obvious exogenous signal by keeping only reads that do not
# map to a panel of common contaminant references for downstream endogenous work.
###############################################################################

# Declare associative array of contaminant references and URLs
declare -A contaminants=(
    ["Vanrija_pseudolonga"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/906/515/GCF_020906515.1_ASM2090651v1/GCF_020906515.1_ASM2090651v1_genomic.fna.gz"
    ["Penicillium_coprophilum"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/826/855/GCF_028826855.1_ASM2882685v1/GCF_028826855.1_ASM2882685v1_genomic.fna.gz"
    ["Homo_sapiens"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    ["Malassezia_restricta"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/290/485/GCF_003290485.1_ASM329048v1/GCF_003290485.1_ASM329048v1_genomic.fna.gz"
    ["Aspergillus_cristatus"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCA/044/706/195/GCA_044706195.1_ASM4470619v1/GCA_044706195.1_ASM4470619v1_genomic.fna.gz"
)

mkdir -p "${WD}/results/contaminants"
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

mkdir -p "${WD}/results/contaminants/mappings"
mkdir -p "${WD}/data/contaminants"

# Download contaminant reference genomes if not present
for name in "${!contaminants[@]}"; do
    url="${contaminants[$name]}"
    ref_file="${WD}/data/contaminants/${name}.fna.gz"
    if [ ! -f "${ref_file}" ]; then
        echo "Downloading reference for ${name} from ${url}..."
        wget -q -O "${ref_file}" "${url}"
    fi
done

# Create joint reference file by concatenating all contaminant references
joint_ref="${WD}/data/contaminants/joint_reference.fna.gz"
if [ ! -f "${joint_ref}" ]; then
    echo "Creating joint reference file..."
    cat "${WD}/data/contaminants/"*.fna.gz >"${joint_ref}"
fi

# Map merged reads to joint contaminant reference and retain unmapped reads
minimap2 -ax sr --secondary=no -t 200 \
    "${joint_ref}" \
    "${WD}/data/trimmed2/Tityra_leucura_merged.fastq.gz" | awk '($5 < 20 || and($2,4)) || $1 ~ /^@/' | samtools view -bS - >"${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"
samtools index "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"

# Extract unmapped merged reads
samtools fastq "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam" \
    | pigz >"${WD}/data/trimmed2/Tityra_leucura_contaminants_merged_unmapped.fastq.gz"

# Map paired-end reads to the contaminant reference for a parallel contamination view
cat >"${WD}/shell/run_minimap2.sh" <<EOF
#!/bin/bash
#PBS -S /bin/bash
#PBS -N minimap2_contaminants
#PBS -o ${WD}/data/log_cluster.txt
#PBS -j oe
#PBS -l select=1:ncpus=100:mem=200gb

source /opt/anaconda3/etc/profile.d/conda.sh
module load NGSmapper/minimap2-2.17
module load Tools/samtools-1.18

minimap2 -ax sr --secondary=no -t 100 \
    ${WD}/data/contaminants/joint_reference.fna.gz \
    ${WD}/data/trimmed2/Tityra_leucura_1_trimmed.fastq.gz \
    ${WD}/data/trimmed2/Tityra_leucura_2_trimmed.fastq.gz \
    | samtools view -bS - \
    | samtools sort -o ${WD}/results/contaminants/mappings/Tityra_leucura_contaminants.bam -

samtools index ${WD}/results/contaminants/mappings/Tityra_leucura_contaminants.bam
EOF
qsub "${WD}/shell/run_minimap2.sh"

samtools coverage "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants.bam" \
    > "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_coverage.txt"

###############################################################################
# 7. Run the UCE workflow
# This step launches the UCE mapping and phylogeny companion workflow to assess
# locus-level capture and place the specimen within a broader Tityrinae context.
###############################################################################

echo "Calling shell/UCEs.sh..."
bash "${WD}/shell/UCEs.sh"

###############################################################################
# 8. Run AutDeNovo assembly and annotation
# This step performs a de novo assembly of the cleaned reads and annotates the
# resulting contigs as the basis for downstream mitochondrial and comparative work.
# Documentation: https:/github.com/nhmvienna/AutDeNovo
###############################################################################
/media/inter/pipelines/AutDeNovo/AutDeNovo_exp.sh \
    Name=Tityra \
    OutputFolder="${WD}/results/denovo" \
    Fwd="${WD}/data/trimmed/Tityra_leucura_contaminants_merged_unmapped.fastq.gz" \
    threads=150 \
    RAM=200 \
    RAMAssembly=1000 \
    decont=no \
    SmudgePlot=no \
    BLASTdb=/media/scratch/NCBI_nt_DB_210714/nt \
    BuscoDB=vertebrata_odb10 \
    Taxdump=/media/scratch/NCBI_taxdump/ \
    Racon=4

###############################################################################
# 9. Assemble the mitochondrial genome
# This step uses the decontaminated read set to reconstruct the mitogenome and
# generate a clean mitochondrial assembly for gene-level analysis.
###############################################################################
echo "Assembling mitochondrial genome..."
mkdir -p "${WD}/results/mitogenome"

sh ${WD}/shell/mitofinder_assembly.sh

###############################################################################
# 10. Estimate sequencing depth on the assembled mitochondrial genome
# This step maps the cleaned reads back to the assembled mitogenome to confirm
# that the recovered mitochondrial contig is supported by the expected depth.
###############################################################################
module load NGSmapper/bwa-0.7.13

bwa index ${WD}/results/mitofinder_assembly/Tityra_leucura/Tityra_leucura/Tityra_leucura_MitoFinder_mitfi_Final_Results/Tityra_leucura_mtDNA_contig_1.fasta

bwa mem -T 100 ${WD}/results/mitofinder_assembly/Tityra_leucura/Tityra_leucura/Tityra_leucura_MitoFinder_mitfi_Final_Results/Tityra_leucura_mtDNA_contig_1.fasta \
    ${WD}/data/trimmed2/Tityra_leucura_contaminants_merged_unmapped.fastq.gz \
    | samtools sort -@ 50 -o ${WD}/results/mitofinder_assembly/mapped.bam

samtools index ${WD}/results/mitofinder_assembly/mapped.bam

samtools coverage ${WD}/results/mitofinder_assembly/mapped.bam > ${WD}/results/mitofinder_assembly/coverage.txt

###############################################################################
# 11. Retrieve additional mitochondrial marker sequences from related Tyrranidae
# This step expands the reference set by downloading COI and other mitochondrial
# loci from related taxa so the phylogenetic comparison has broader taxonomic scope.
###############################################################################

mkdir ${WD}/data/Mito
for j in COI CO2 ND1 ND2 ND3 ND4 ND4L ND5 N6 CYTB ATP8 COX1 COX2 COX3 ATP6 rrnL rrnS; do

    for i in "Onychorhynchus coronatus" Tityra Schiffornis Laniocera Iodopleura Laniisoma Xenopsaris Pachyramphus; do
        
        python /media/inter/mkapun/projects/EfetchThePython/EfetchThePython.py \
            --email capopony@gmail.com \
            --api_key dfbe3972fc58b8478eba0929f14847f9cd08 \
            --Term "\"$i\"[Organism], $j" \
            --Output ${WD}/data/Mito/${j}_$i \
            --FASTA
    done
done

###############################################################################
# 12. Combine mitochondrial gene sequences
# This step merges the downloaded taxon-specific markers with the de novo
# Tityra-derived mitochondrial loci into a single gene-by-gene dataset.
###############################################################################

echo ""
echo "==================================================================="
echo "Combining mitochondrial gene sequences (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/combine_mito_genes_assembly.sh

###############################################################################
# 13. Run species-level phylogenetic analysis on the combined mitochondrial genes
# This step reconstructs the first phylogenetic tree from the combined gene set
# to evaluate relationships among species before concatenation.
###############################################################################

echo ""
echo "==================================================================="
echo "Running phylogenetic analysis (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/phylogenetic_analysis_assembly.sh

###############################################################################
# 14. Concatenate genes by species for the assembly-based matrix
# This step builds a species-level concatenated nucleotide matrix that can be
# analyzed jointly for improved phylogenetic resolution.
###############################################################################
echo ""
echo "==================================================================="
echo "Concatenating mitochondrial genes by species (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/concatenate_genes_assembly.sh

###############################################################################
# 15. Aggregate mitochondrial gene metadata for all analysis samples
# This step collates gene-specific accession and taxon metadata into a single
# table for downstream interpretation of the concatenated mitochondrial dataset.
###############################################################################

echo ""
echo "==================================================================="
echo "Aggregating mitochondrial gene metadata (assembly-based)..."
echo "==================================================================="

bash "${WD}/shell/aggregate_mito_metadata.sh"

###############################################################################
# 16. Run phylogenetic analysis on the concatenated gene matrix
# This step uses the concatenated mitochondrial matrix to infer the final
# species-level tree from the combined loci.
###############################################################################
echo ""
echo "==================================================================="
echo "Running phylogenetic analysis on concatenated genes (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/phylogenetic_analysis_concatenated_assembly.sh

###############################################################################
# 17. Run phylogenetic analysis on ND2 gene sequences only using the 2 July 2026 alignment
# This step uses the ND2 gene alignment to infer a species-level tree for comparison
# with the concatenated mitochondrial matrix.
###############################################################################
echo ""
echo "==================================================================="
echo "Running phylogenetic analysis on ND2 gene sequences..."
echo "==================================================================="

ALN_DIR=${WD}/data/alignments_2July2026
TREE_FILE=${ALN_DIR}/ND2.aln.fasta.contree
XLSX_FILE=${ALN_DIR}/samples_for_tit_leuc.xlsx
TIP_MAP=${ALN_DIR}/tip_labels.tsv
PLOT_DIR=${WD}/results/ND2_tree

mkdir -p "${PLOT_DIR}"

# Generate the old->new tip-label mapping if not already present
if [ ! -f "${TIP_MAP}" ]; then
    python3 ${WD}/scripts/generate_tip_labels.py \
        --tree "${TREE_FILE}" \
        --xlsx "${XLSX_FILE}" \
        --out  "${TIP_MAP}"
fi

Rscript ${WD}/scripts/plot_tree.r \
    ${TREE_FILE} \
    ${PLOT_DIR} \
    "Manacus manacus KF228541.1" \
    "nd2" \
    "${TIP_MAP}"

###############################################################################
# 18. Download complete mitogenomes from GenBank/RefSeq
# This step adds published complete mitochondrial genome references to provide
# a broader comparative backbone for the final mitogenome analysis.
###############################################################################
echo ""
echo "==================================================================="
echo "Downloading complete mitochondrial genomes from GenBank..."
echo "==================================================================="

bash ${WD}/shell/download_complete_mitogenomes.sh

echo ""
echo "==================================================================="
echo "Reorienting and circularizing Tityra mitochondrial genome..."
echo "==================================================================="

bash ${WD}/shell/reorient_circularize_tityra.sh

###############################################################################
# 19. Recreate the circular mitogenome map overview
# This step documents the manual OGDRAW workflow used to visualize the circular
# mitochondrial genome once the complete reference sequence is in place.
###############################################################################
# CircularPlot.pdf was generated using the OGDRAW web service:
#   https://chlorobox.mpimp-golm.mpg.de/OGdraw.html
#
# To recreate:
#   1. Go to https://chlorobox.mpimp-golm.mpg.de/OGDraw.html
#   2. Upload: ${WD}/data/refseq_genbank/Tityra_leucura.gb
#   3. Select output format: PDF
#   4. Click "Draw" and download the result to:
#      ${WD}/results/CircularPlot.pdf
echo ""
echo "==================================================================="
echo "NOTE: CircularPlot.pdf must be recreated manually via OGDRAW:"
echo "  https://chlorobox.mpimp-golm.mpg.de/OGdraw.html"
echo "  Input: ${WD}/data/refseq_genbank/Tityra_leucura.gb"
echo "  Output: ${WD}/results/CircularPlot.pdf"
echo "==================================================================="

###############################################################################
# 20. Align the complete mitochondrial genomes
# This step prepares a global alignment of complete mitogenomes so synteny and
# phylogenetic comparisons can be made across taxa.
###############################################################################
echo ""
echo "==================================================================="
echo "Aligning complete mitochondrial genomes..."
echo "==================================================================="

bash ${WD}/shell/align_complete_mitogenomes.sh

###############################################################################
# 21. Visualize the complete mitogenome alignments
# This step converts the aligned mitogenome data into plot-ready summaries for
# inspection of structural consistency and divergence.
###############################################################################
echo ""
echo "==================================================================="
echo "Visualizing complete mitogenome alignments..."
echo "==================================================================="

bash ${WD}/shell/visualize_mitogenome_alignment.sh

###############################################################################
# 22. Run phylogenetic analysis on the complete mitogenomes
# This step infers a final whole-mitogenome tree to compare the assembled Tityra
# genome with complete mitogenome references from related taxa.
###############################################################################
echo ""
echo "==================================================================="
echo "Phylogenetic analysis of complete mitogenomes..."
echo "==================================================================="

bash ${WD}/shell/phylogenetic_analysis_complete_mitogenome.sh

###############################################################################
# 23. Run mitochondrial genome synteny analysis
# This step compares genomic organization across mitogenomes to check for major
# structural rearrangements or concordance with the expected gene order.
###############################################################################
echo ""
echo "==================================================================="
echo "Mitochondrial genome synteny analysis..."
echo "==================================================================="

bash ${WD}/shell/synteny_analysis.sh

echo ""
echo "==================================================================="
echo "All analyses complete!"
echo "==================================================================="

###############################################################################
# 24. Plot the GBIF occurrence map
# This final step summarizes the geographic distribution of the focal species to
# place the molecular results in a broader biogeographic context.
###############################################################################
echo ""
echo "==================================================================="
echo "Plotting GBIF occurrence map for Tityra inquisitor and T. leucura..."
echo "==================================================================="

Rscript ${WD}/scripts/plot_GBIF_map.R ${WD}

echo ""
echo "==================================================================="
echo "GBIF map complete!"
echo "==================================================================="

