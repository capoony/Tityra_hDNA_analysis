
###############################################################################
# 9. Combine mitochondrial genes by individual specimens
###############################################################################
echo ""
echo "==================================================================="
echo "Combining mitochondrial genes by individual specimens..."
echo "==================================================================="

bash ${WD}/shell/combine_individual_genes.sh

###############################################################################
# 10. Phylogenetic analysis of individual-level concatenated genes
###############################################################################
echo ""
echo "==================================================================="
echo "Running phylogenetic analysis (individual-level)..."
echo "==================================================================="

bash ${WD}/shell/individual_phylogenetic_analysis.sh

###############################################################################
# 8. Phylogenetic Analysis
###############################################################################
# Perform multiple alignment, phylogenetic inference, and tree plotting
echo ""
echo "==================================================================="
echo "Performing phylogenetic analysis..."
echo "==================================================================="

bash ${WD}/shell/phylogenetic_analysis.sh

###############################################################################
# 9. Run AutDeNovo Pipeline
# for documentation, see: https:/github.com/nhmvienna/AutDeNovo
###############################################################################
# Run AutDeNovo pipeline for de novo assembly and annotation
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
# 9. Graph-Based Alignment to Reduce Reference Bias
###############################################################################

echo "Step 9: Building variation graph from multiple references..."
echo "Note: This step requires Step 10 (BUSCO on references) to be run first"
echo "Running graph alignment pipeline..."

# Run the standalone graph alignment script
bash ${WD}/shell/graph_alignment_busco.sh ${WD} 50

echo "Graph-based alignment complete!"

###############################################################################
# 10. Make Phylogeny with BUSCO genes (Graph-Based + De Novo Assembly)
###############################################################################

echo "Step 10: Making phylogeny with BUSCO genes (including graph-based sequences)..."

mkdir -p ${WD}/results/phylogeny/data # Create directory for phylogeny data

# Clear genomes.names to avoid duplicates from multiple runs
> ${WD}/results/phylogeny/data/genomes.names

## download reference genomes (commented out, see genomes.txt for URLs)
while IFS=$"," read -r Name ID URL; do
    #curl -L -o "${WD}/results/phylogeny/data/${Name}.fna.gz" "${URL}"
    echo ${Name} >>${WD}/results/phylogeny/data/genomes.names
done <${WD}/results/phylogeny/data/genomes.txt

## copy Pachyrhampus genome to phylogeny data directory
cp ${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz \
    ${WD}/results/phylogeny/data/Pachyramphus_minor.fna.gz
echo "Pachyramphus_minor" >>${WD}/results/phylogeny/data/genomes.names

## copy Tityra genome (de novo assembly) to phylogeny data directory
cp ${WD}/results/denovo/output/Tityra_ILL.fa.gz \
    ${WD}/results/phylogeny/data/Tityra_leucura.fna.gz
echo "Tityra_leucura" >>${WD}/results/phylogeny/data/genomes.names

## copy Tityra graph-based BUSCO sequences (avoids reference bias)
mkdir -p ${WD}/results/phylogeny/data/Tityra_leucura_graph

# Check if graph-based BUSCO sequences exist
if [ -d "${WD}/results/phylogeny_graph/busco_aa" ] && [ $(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | wc -l) -gt 0 ]; then
    # Combine all graph-based BUSCO protein sequences into single file
    cat ${WD}/results/phylogeny_graph/busco_aa/*.faa > \
        ${WD}/results/phylogeny/data/Tityra_leucura_graph/busco_set.fasta
    echo "Tityra_leucura_graph" >>${WD}/results/phylogeny/data/genomes.names
    echo "Including graph-based Tityra sequences in phylogeny"
else
    echo "Warning: Graph-based BUSCO sequences not found - run Step 9 first"
    echo "Skipping graph-based sequences in phylogeny"
fi

mkdir -p ${WD}/results/phylogeny/BUSCO # Create BUSCO output directory

# Run BUSCO for each genome to identify single-copy orthologs
while IFS=$"," read -r Name; do

    # Skip graph-based sample - already has protein sequences from graph alignment
    if [ "${Name}" == "Tityra_leucura_graph" ]; then
        echo "Skipping BUSCO run for ${Name} - will copy graph-based sequences directly"
        continue
    fi

    # Check if BUSCO has already been run successfully for this genome
    if [ -d "${WD}/results/phylogeny/BUSCO/${Name}/run_aves_odb10" ] && \
       [ -f "${WD}/results/phylogeny/BUSCO/${Name}/run_aves_odb10/full_table.tsv" ]; then
        echo "BUSCO already completed for ${Name}, skipping..."
        continue
    fi

    echo """

  #!/bin/sh

  ## name of Job
  #PBS -N BUSCO_${Name}

  ## Redirect output stream to this file.
  #PBS -o ${WD}/results/phylogeny/BUSCO/${Name}_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select a maximum of 200 cores and 1000gb of RAM
  #PBS -l select=1:ncpus=50:mem=500g

  ######## load dependencies #######

  source /opt/anaconda3/etc/profile.d/conda.sh
  conda activate busco_5.4.3

  ######## run analyses #######

  cd ${WD}/results/phylogeny/BUSCO

  pigz -d ${WD}/results/phylogeny/data/${Name}.fna.gz

  busco -i ${WD}/results/phylogeny/data/${Name}.fna \
    -o ${Name} \
    -m genome \
    -c 200 \
    -l aves_odb10

  pigz -f ${WD}/results/phylogeny/data/${Name}.fna

  """ >${WD}/results/phylogeny/BUSCO/${Name}.sh

    sh ${WD}/results/phylogeny/BUSCO/${Name}.sh

done <${WD}/results/phylogeny/data/genomes.names

## concatenate BUSCO results
mkdir -p ${WD}/results/phylogeny/concatenated

cd ${WD}/results/phylogeny/BUSCO

## identify all "complete" BUSCO genes across all genomes and concatenate their IDs
for file in $(find -iname "full_table*.tsv"); do
    grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >>${WD}/results/phylogeny/concatenated/complete_busco_ids.txt
done

## Add graph-based BUSCO genes to the list (if they exist)
if [ -d "${WD}/results/phylogeny_graph/busco_aa" ]; then
    # Extract gene IDs from graph-based protein files
    for faa in ${WD}/results/phylogeny_graph/busco_aa/*.faa; do
        if [ -f "${faa}" ]; then
            basename ${faa} .faa >> ${WD}/results/phylogeny/concatenated/complete_busco_ids.txt
        fi
    done
    echo "Added $(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | wc -l) graph-based genes to BUSCO list"
fi

## sort all BUSCO IDs and count occurrences
sort ${WD}/results/phylogeny/concatenated/complete_busco_ids.txt |
    uniq -c \
        >${WD}/results/phylogeny/concatenated/complete_busco_ids_with_counts.txt

# Count number of genomes being compared
NUM_GENOMES=$(wc -l < ${WD}/results/phylogeny/data/genomes.names)

## filter for BUSCO genes that are present in all genomes
awk -v ng="$NUM_GENOMES" '$1 == ng {print $2}' ${WD}/results/phylogeny/concatenated/complete_busco_ids_with_counts.txt \
    >${WD}/results/phylogeny/concatenated/all_genomes_busco_ids.txt

echo "Found $(wc -l < ${WD}/results/phylogeny/concatenated/all_genomes_busco_ids.txt) BUSCO genes present in all ${NUM_GENOMES} genomes"

## Further filter to only include genes present in graph-based dataset
if [ -d "${WD}/results/phylogeny_graph/busco_aa" ]; then
    # Create list of graph-based genes
    ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | xargs -n1 basename | sed 's/.faa$//' > ${WD}/results/phylogeny/concatenated/graph_genes.txt
    
    # Keep only genes that are in both lists
    comm -12 \
        <(sort ${WD}/results/phylogeny/concatenated/all_genomes_busco_ids.txt) \
        <(sort ${WD}/results/phylogeny/concatenated/graph_genes.txt) \
        >${WD}/results/phylogeny/concatenated/final_busco_ids.txt
    
    echo "Filtered to $(wc -l < ${WD}/results/phylogeny/concatenated/final_busco_ids.txt) genes present in graph-based dataset"
else
    # No graph-based data, use all genes
    cp ${WD}/results/phylogeny/concatenated/all_genomes_busco_ids.txt ${WD}/results/phylogeny/concatenated/final_busco_ids.txt
    echo "Warning: No graph-based data found, using all $(wc -l < ${WD}/results/phylogeny/concatenated/final_busco_ids.txt) genes"
fi

mkdir -p ${WD}/results/phylogeny/concatenated/busco_aa

# Copy and rename BUSCO protein sequences for each genome and gene
while IFS=$"," read -r Name; do

    # For graph-based sequences, copy directly from phylogeny_graph/busco_aa
    if [ "${Name}" == "Tityra_leucura_graph" ]; then
        while IFS=$"," read -r gene; do
            # Check if this gene exists in graph-based results
            if [ -f "${WD}/results/phylogeny_graph/busco_aa/${gene}.faa" ]; then
                cp ${WD}/results/phylogeny_graph/busco_aa/${gene}.faa \
                    ${WD}/results/phylogeny/concatenated/busco_aa/${Name}_${gene}
                sed -i 's/^>/>'${Name}'|/g' ${WD}/results/phylogeny/concatenated/busco_aa/${Name}_${gene}
            else
                echo "Warning: Gene ${gene} not found in graph-based results, skipping..."
            fi
        done <${WD}/results/phylogeny/concatenated/final_busco_ids.txt
        continue
    fi

    while IFS=$"," read -r gene; do

        # Check if gene file exists before copying (graph-based may not have all genes)
        if [ ! -f "${WD}/results/phylogeny/BUSCO/${Name}/run_aves_odb10/busco_sequences/single_copy_busco_sequences/${gene}.faa" ]; then
            echo "Warning: Gene ${gene} not found for ${Name}, skipping..."
            continue
        fi

        cp ${WD}/results/phylogeny/BUSCO/${Name}/run_aves_odb10/busco_sequences/single_copy_busco_sequences/${gene}.faa \
            ${WD}/results/phylogeny/concatenated/busco_aa/${Name}_${gene}
        sed -i 's/^>/>'${Name}'|/g' ${WD}/results/phylogeny/concatenated/busco_aa/${Name}_${gene}
    done <${WD}/results/phylogeny/concatenated/final_busco_ids.txt

done <${WD}/results/phylogeny/data/genomes.names

### concatenate and reduce to shared genes

mkdir -p ${WD}/results/phylogeny/prealigned

# Concatenate protein sequences for each BUSCO gene
while read gene; do
    # Check if at least one genome has this gene
    if ls ${WD}/results/phylogeny/concatenated/busco_aa/*_${gene} 1> /dev/null 2>&1; then
        cat ${WD}/results/phylogeny/concatenated/busco_aa/*_${gene} \
            >${WD}/results/phylogeny/prealigned/${gene}_aa.fasta
    else
        echo "Warning: No sequences found for gene ${gene}, skipping..."
    fi
done <${WD}/results/phylogeny/concatenated/final_busco_ids.txt

mkdir -p ${WD}/results/phylogeny/mafft

conda activate mafft-7.487

# Align each BUSCO gene with MAFFT
for i in ${WD}/results/phylogeny/prealigned/*_aa.fasta; do
    # Check if file exists (pattern might not match any files)
    [ -f "$i" ] || continue

    TMP=${i##*/}
    ID=${TMP%_*}

    mafft \
        --thread 50 \
        --auto \
        ${i} \
        >${WD}/results/phylogeny/mafft/${ID}_aln.fasta

done

# Fix sequence IDs after MAFFT alignment
for i in ${WD}/results/phylogeny/prealigned/*_aa.fasta; do
    # Check if file exists (pattern might not match any files)
    [ -f "$i" ] || continue

    TMP=${i##*/}
    ID=${TMP%_*}

    python /media/inter/mkapun/projects/TorpedoDeNovo/scripts/fixIDAfterMafft.py \
        --Alignment ${WD}/results/phylogeny/mafft/${ID}_aln.fasta \
        --input ${i} \
        >${WD}/results/phylogeny/mafft/${ID}_aln_fixed.fasta

    rm -f ${WD}/results/phylogeny/mafft/${ID}_aln.fasta

done

## make phylogeny

mkdir ${WD}/results/phylogeny/phylogeny

# Concatenate all alignments for phylogenetic analysis
python /media/inter/mkapun/projects/TorpedoDeNovo/scripts/proteins2genome.py \
    --input ${WD}/results/phylogeny/mafft \
    >${WD}/results/phylogeny/phylogeny/alignment.fa

module load Phylogeny/RAxML-2.8.10

## make new directory

cd ${WD}/results/phylogeny/phylogeny

## run ML tree reconstruction with RAxML
raxmlHPC-PTHREADS-SSE3 \
    -m PROTGAMMAWAG \
    -N 20 \
    -p 772374015 \
    -n Tityra \
    -s ${WD}/results/phylogeny/phylogeny/alignment.fa \
    -o Acanthisitta_chloris \
    -T 200

# Run RAxML bootstrapping
raxmlHPC-PTHREADS-SSE3 \
    -m PROTGAMMAWAG \
    -N autoMRE \
    -p 772374015 \
    -b 444353738 \
    -n bootrep \
    -s ${WD}/results/phylogeny/phylogeny/alignment.fa \
    -o Acanthisitta_chloris \
    -T 200

# Reconcile best ML tree with bootstrap replicates
raxmlHPC-SSE3 -f b \
    -m GTRGAMMA \
    -t RAxML_bestTree.Tityra \
    -z RAxML_bootstrap.bootrep \
    -n FINAL \
    -o Acanthisitta_chloris

# Plot phylogenetic tree using R and ggtree
Rscript -e """
# load necessary R libraries
library('ggtree')
library('gridExtra')
library('ggrepel')
library('ape')
library('ggplot2')
library('dplyr')

## load tree file and root with outgroup taxa
tree<-read.tree('${WD}/results/phylogeny/phylogeny/RAxML_bipartitions.FINAL')

## root tree with outgroup
tree<-root(tree,outgroup='Acanthisitta_chloris')

## calculate tree height (on x-axis) using ape functions
tree_depth <- max(node.depth.edgelength(tree))

## plot tree
PLOT.tree<-ggtree(tree)+
  ggtitle('Tityra_BUSCO')+
  theme_tree2()+
  theme_bw()+
  ggplot2::xlim(0, tree_depth+0.1)+
  xlab('av. subst./site') +
  geom_nodelab()+
  theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  geom_tiplab()

## export tree
ggsave(filename='${WD}/results/phylogeny/phylogeny/Tityra_BUSCO.pdf',
  PLOT.tree)
ggsave(filename='${WD}/results/phylogeny/phylogeny/Tityra_BUSCO.png',
  PLOT.tree)
"""

## OK, now repeat with mapping reads to Pachyrhamphus minor BUSCO genes
sh ${WD}/shell/BUSCO_Pminor.sh ${WD}

## OK, now repeat with mapping reads to Pachyrhamphus minor BUSCO genes and use ALL genes
sh ${WD}/shell/BUSCO_Pminor_allBUSCO.sh ${WD}

## OK, now repeat with mapping Tityra reads to Ochrsit BUSCO genes
sh ${WD}/shell/BUSCO_Ocrist.sh ${WD}

## OK, now repeat with mapping Tityra reads to Ochrsit BUSCO genes
sh ${WD}/shell/BUSCO_Ocrist_allBUSCO.sh ${WD}

## OK, now repeat with mapping Tityra reads to Ochrsit BUSCO genes
sh ${WD}/shell/BUSCO_tsav.sh ${WD}

## Make Phylogenies for all genes based on the denovo assembly 
sh ${WD}/shell/BUSCO_phylogeny_denovo_indgenes.sh ${WD}

## Make Phylogenies for all genes based on the graph-based alignment
echo "Building phylogenies for graph-based BUSCO genes..."
if [ -d "${WD}/results/phylogeny_graph/busco_aa" ] && [ $(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | wc -l) -gt 0 ]; then
    sh ${WD}/shell/BUSCO_phylogeny_graph_indgenes.sh ${WD}
else
    echo "Warning: Graph-based BUSCO genes not found - skipping graph phylogeny"
    echo "Run Step 9 (graph_alignment_busco.sh) first to generate graph-based sequences"
fi
