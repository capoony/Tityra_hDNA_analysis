#!/bin/bash
###############################################################################
# Phylogeny Pipeline for Graph-Based BUSCO Genes
# Processes graph-based consensus sequences for phylogenetic analysis
###############################################################################

WD=$1

# Check if graph-based BUSCO directory exists
if [ ! -d "${WD}/results/phylogeny_graph/busco_aa" ]; then
    echo "Error: Graph-based BUSCO directory does not exist!"
    echo "Please run Step 9 (graph_alignment_busco.sh) first"
    exit 1
fi

# Count available graph-based BUSCO genes
NUM_GENES=$(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | wc -l)

if [ ${NUM_GENES} -eq 0 ]; then
    echo "Error: No graph-based BUSCO genes found!"
    echo "Please run Step 9 (graph_alignment_busco.sh) first"
    exit 1
fi

echo "Processing ${NUM_GENES} graph-based BUSCO genes for phylogeny..."

# Get list of genes from graph-based BUSCO
cd ${WD}/results/phylogeny_graph/busco_aa
ls *.faa | sed 's/.faa//g' > ${WD}/results/phylogeny_graph/gene_list.txt

# For each gene, combine graph-based Tityra with reference sequences and build tree
while IFS=$"," read -r gene; do
    echo "Processing gene: ${gene}"
    
    mkdir -p ${WD}/results/phylogeny_graph/phylogeny/${gene}

    # Check if this gene exists in the regular phylogeny pipeline
    if [ ! -f "${WD}/results/phylogeny/busco_aa/${gene}.faa" ]; then
        echo "Warning: Gene ${gene} not found in regular phylogeny pipeline, skipping..."
        continue
    fi

    # Combine graph-based Tityra sequence with reference sequences
    cd ${WD}/results/phylogeny_graph/phylogeny/${gene}
    
    # Copy reference sequences from regular phylogeny
    cat ${WD}/results/phylogeny/busco_aa/${gene}.faa > ${gene}_combined.fasta
    
    # Add graph-based Tityra sequence
    cat ${WD}/results/phylogeny_graph/busco_aa/${gene}.faa >> ${gene}_combined.fasta

    # Align with MAFFT
    conda activate /media/inter/mkapun/envs/mafft-7.487
    
    mafft --auto --thread 4 ${gene}_combined.fasta > ${gene}_aligned.fasta 2> /dev/null
    
    conda deactivate

    # Clean up sequence names (remove everything after first "|")
    sed -i 's/|.*//g' ${gene}_aligned.fasta

    # Run ML tree reconstruction with RAxML
    raxmlHPC-PTHREADS-SSE3 \
        -m PROTGAMMAWAG \
        -N 20 \
        -p 772374015 \
        -n Tityra_graph \
        -s ${gene}_aligned.fasta \
        -o Acanthisitta_chloris \
        -T 200

    # Bootstrap replicates
    raxmlHPC-PTHREADS-SSE3 \
        -m PROTGAMMAWAG \
        -N 25 \
        -p 772374015 \
        -b 444353738 \
        -n bootrep \
        -s ${gene}_aligned.fasta \
        -o Acanthisitta_chloris \
        -T 200

    # Reconcile best ML tree with bootstraps
    raxmlHPC-SSE3 -f b \
        -m GTRGAMMA \
        -t RAxML_bestTree.Tityra_graph \
        -z RAxML_bootstrap.bootrep \
        -n FINAL \
        -o Acanthisitta_chloris

    # Plot tree with R
    Rscript -e """
library('ggtree')
library('gridExtra')
library('ggrepel')
library('ape')
library('ggplot2')
library('phangorn')
library('dplyr')

## read tree
tree <- read.tree('RAxML_bipartitions.FINAL')

## caluculate tree height (on x-axis)
Xmax<-max(nodeHeights(tree))

tree<-root(tree,outgroup='Acanthisitta_chloris')

## plot tree
PLOT.tree<-ggtree(tree)+
  ggtitle('${gene} (Graph-based)')+
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
ggsave(filename='${gene}_graph_BUSCO.pdf',
  PLOT.tree)
ggsave(filename='${gene}_graph_BUSCO.png',
  PLOT.tree)
"""

    echo "Completed phylogeny for gene: ${gene}"

done < ${WD}/results/phylogeny_graph/gene_list.txt

echo "==================================================================="
echo "Graph-based phylogeny pipeline complete!"
echo "Results in: ${WD}/results/phylogeny_graph/phylogeny/"
echo "==================================================================="
