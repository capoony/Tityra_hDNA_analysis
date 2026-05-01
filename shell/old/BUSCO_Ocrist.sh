
WD=$1
## OK, last, but not least try mapping reads to Pachyramphus minor BUSCO genes
mkdir ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO

## concatenate single-copy BUSCO sequences and replace header with name of *.fna file
cd ${WD}/results/phylogeny/BUSCO/Oxyruncus_cristatus/run_aves_odb10/busco_sequences/single_copy_busco_sequences/
for file in *.fna; do
    name=$(basename "$file" .fna)
    #echo $name
    sed "s/^>/>${name} /g" "$file" | awk '{print $1}' >> ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Oxyruncus_cristatus.fna
done

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

pigz -f ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Oxyruncus_cristatus.fna

## map reads to Pachyramphus minor BUSCO genes
minimap2 -ax sr --secondary=no -t 200 \
    ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Oxyruncus_cristatus.fna.gz \
   ${WD}/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz |
    samtools view -bS -F 4 - |
    samtools sort -o ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.bam

samtools index ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.bam

pigz -d ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Oxyruncus_cristatus.fna.gz

## check mapping in IGViewer
#sh /opt/bioinformatics/IGV_Linux_2.10.2/igv.sh 

samtools coverage --reference ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Oxyruncus_cristatus.fna.gz \
    ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.bam \
    > ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.coverage.txt

## only retain genes with at least 2x average read depth and >= 60% coverage
awk '$7 >= 2 && $6 >= 60' ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.coverage.txt \
    > ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.coverage_filtered.txt

awk 'NR>1 {print $1}' ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.coverage_filtered.txt \
    > ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.genes.txt

## call variants

mkdir -p ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/consensus

samtools faidx ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Oxyruncus_cristatus.fna


while IFS=$'\n' read -r gene; do
    echo "Processing gene: $gene"

# Input arguments
BAM=${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.bam
REF=${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Oxyruncus_cristatus.fna
REGION=$(awk -v gene=${gene} '$1 == gene {print $1 ":1-" $2}' ${REF}.fai)
echo "Region: $REGION"
OUT=${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/consensus/Tityra_leucura_${gene}.fasta

# Temporary files
VCF_TMP=$(mktemp --suffix=.vcf)
BED_TMP=$(mktemp --suffix=.bed)
REF_TMP=$(mktemp --suffix=.fa)


# Step 1: Create VCF for the region
bcftools mpileup -Ou -f "$REF" -r "$REGION" "$BAM" | \
bcftools call -Ou -mv -Ov -o "$VCF_TMP"

bgzip -c "$VCF_TMP" > "${VCF_TMP}.gz"
tabix -p vcf "${VCF_TMP}.gz"

# Step 2: Generate BED of uncovered positions (coverage < 1)
samtools depth -a -r "$REGION" "$BAM" | \
awk '$3 < 1 {print $1"\t"($2-1)"\t"$2}' > "$BED_TMP"

# Step 3: Extract reference sequence for the region
samtools faidx "$REF" "$REGION" > "$REF_TMP"

# Step 4: Create consensus with variants and 'N' for no coverage
cat "$REF_TMP" | bcftools consensus --haplotype A -m "$BED_TMP" "${VCF_TMP}.gz" \
| sed "s/^>/>Tityra_leucura|/g" > "$OUT"

# Clean up
rm "$VCF_TMP" "$BED_TMP" "$REF_TMP"

echo "✅ Consensus saved to: $OUT"
done < ${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.genes.txt

###############################################################################
# 9. Make Phylogeny with BUSCO genes based on nucleotides and by aligning Tityra to the Pachyramphus minor BUSCO genes
###############################################################################

echo "Step 9: Making phylogeny with BUSCO genes..."

mkdir -p ${WD}/results/phylogeny_DNA_Ocrist/concatenated

cd ${WD}/results/phylogeny_DNA_Ocrist/BUSCO

## only use the gene names that are present in the Pachyramphus minor BUSCO genes and which have sufficient coverage in the Tityra reads

mkdir -p ${WD}/results/phylogeny_DNA_Ocrist/concatenated/busco_dna

while IFS=$"," read -r Name; do

    while IFS=$"," read -r gene; do

        cp ${WD}/results/phylogeny/BUSCO/${Name}/run_aves_odb10/busco_sequences/single_copy_busco_sequences/${gene}.fna \
            ${WD}/results/phylogeny_DNA_Ocrist/concatenated/busco_dna/${Name}_${gene}
        sed -i 's/^>/>'${Name}'|/g' ${WD}/results/phylogeny_DNA_Ocrist/concatenated/busco_dna/${Name}_${gene}
    done <${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.genes.txt

done <${WD}/results/phylogeny/data/genomes.names

### concatenate and reduce to shared genes, what a beautiful code!!!
mkdir -p ${WD}/results/phylogeny_DNA_Ocrist/prealigned

while read gene; do
    ## add remaining species from BUSCO results
    cat ${WD}/results/phylogeny_DNA_Ocrist/concatenated/busco_dna/*_${gene}* \
        >>${WD}/results/phylogeny_DNA_Ocrist/prealigned/${gene}_dna.fasta
    ## add consensus sequences from Tityra mapped to Pachyramphus minor BUSCO genes
    cat /media/inter/mkapun/projects/Tityra/results/phylogeny/Oxyruncus_cristatus_BUSCO/consensus/Tityra_leucura_${gene}.fasta \
        >>${WD}/results/phylogeny_DNA_Ocrist/prealigned/${gene}_dna.fasta
done <${WD}/results/phylogeny/Oxyruncus_cristatus_BUSCO/Tityra_leucura_Oxyruncus_cristatus_PE.genes.txt

mkdir ${WD}/results/phylogeny_DNA_Ocrist/mafft

conda activate mafft-7.487
for i in ${WD}/results/phylogeny_DNA_Ocrist/prealigned/*_dna.fasta; do

    TMP=${i##*/}
    ID=${TMP%_*}

    mafft \
        --thread 50 \
        --auto \
        ${i} \
        >${WD}/results/phylogeny_DNA_Ocrist/mafft/${ID}_aln.fasta

done

for i in ${WD}/results/phylogeny_DNA_Ocrist/prealigned/*_dna.fasta; do

    TMP=${i##*/}
    ID=${TMP%_*}

    python /media/inter/mkapun/projects/TorpedoDeNovo/scripts/fixIDAfterMafft.py \
        --Alignment ${WD}/results/phylogeny_DNA_Ocrist/mafft/${ID}_aln.fasta \
        --input ${i} \
        >${WD}/results/phylogeny_DNA_Ocrist/mafft/${ID}_aln_fixed.fasta

    rm -f ${WD}/results/phylogeny_DNA_Ocrist/mafft/${ID}_aln.fasta

done

## make phylogeny

mkdir ${WD}/results/phylogeny_DNA_Ocrist/phylogeny

python /media/inter/mkapun/projects/TorpedoDeNovo/scripts/proteins2genome.py \
    --input ${WD}/results/phylogeny_DNA_Ocrist/mafft \
    >${WD}/results/phylogeny_DNA_Ocrist/phylogeny/alignment.fa

module load Phylogeny/RAxML-2.8.10

## make new directory

cd ${WD}/results/phylogeny_DNA_Ocrist/phylogeny

## run ML tree reconstruction
raxmlHPC-PTHREADS-SSE3 \
    -m GTRGAMMA \
    -N 20 \
    -p 772374015 \
    -n Tityra \
    -s ${WD}/results/phylogeny_DNA_Ocrist/phylogeny/alignment.fa \
    -o Acanthisitta_chloris \
    -T 200

raxmlHPC-PTHREADS-SSE3 \
    -m GTRGAMMA \
    -N 25 \
    -p 772374015 \
    -b 444353738 \
    -n bootrep \
    -s ${WD}/results/phylogeny_DNA_Ocrist/phylogeny/alignment.fa \
    -o Acanthisitta_chloris \
    -T 200

# Now, reconcile the best ML tree w/ the bootreps:
raxmlHPC-SSE3 -f b \
    -m GTRGAMMA \
    -t RAxML_bestTree.Tityra \
    -z RAxML_bootstrap.bootrep \
    -n FINAL \
    -o Acanthisitta_chloris

Rscript -e """
# load necessary R libraries
library('ggtree')
library('gridExtra')
library('ggrepel')
library('ape')
library('ggplot2')
library('phangorn')
library('dplyr')
library(phytools) # to determine the maximum tree height and add midpoint root

## load tree file and root with outgroup taxa
tree<-read.tree('${WD}/results/phylogeny_DNA_Ocrist/phylogeny/RAxML_bipartitions.FINAL')


## plot tree

## caluculate tree height (on x-axis)
Xmax<-max(nodeHeights(tree))

tree<-root(tree,outgroup='Acanthisitta_chloris')

## plot tree
PLOT.tree<-ggtree(tree)+
  ggtitle('Tityra_BUSCO')+
  theme_tree2()+
  theme_bw()+
  ggplot2::xlim(0,
    Xmax+0.25)+
  xlab('av. subst./site') +
  geom_nodelab()+
  theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  geom_tiplab()

## export tree
ggsave(filename='${WD}/results/phylogeny_DNA_Ocrist/phylogeny/Tityra_BUSCO.pdf',
  PLOT.tree)
ggsave(filename='${WD}/results/phylogeny_DNA_Ocrist/phylogeny/Tityra_BUSCO.png',
  PLOT.tree)
"""
