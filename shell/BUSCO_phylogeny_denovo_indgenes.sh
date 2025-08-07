WD=$1

while IFS=$"," read -r gene; do
    mkdir -p ${WD}/results/phylogeny/phylogeny/${gene}

    ## make new directory

    cd ${WD}/results/phylogeny/phylogeny/${gene}

    ## in the fasta files exclude everything after the first "|"
    sed -i 's/|.*//g' ${WD}/results/phylogeny/mafft/${gene}_aln_fixed.fasta

    ## run ML tree reconstruction
    raxmlHPC-PTHREADS-SSE3 \
        -m PROTGAMMAWAG \
        -N 20 \
        -p 772374015 \
        -n Tityra \
        -s ${WD}/results/phylogeny/mafft/${gene}_aln_fixed.fasta \
        -o Acanthisitta_chloris \
        -T 200
    ## only do 25 replicates for the first run, to get a quick idea of the tree
    ## and to avoid running out of memory
    raxmlHPC-PTHREADS-SSE3 \
        -m PROTGAMMAWAG \
        -N 25 \
        -p 772374015 \
        -b 444353738 \
        -n bootrep \
        -s ${WD}/results/phylogeny/mafft/${gene}_aln_fixed.fasta \
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
tree<-read.tree('${WD}/results/phylogeny/phylogeny/${gene}/RAxML_bipartitions.FINAL')


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
ggsave(filename='${WD}/results/phylogeny/phylogeny/Tityra_${gene}.pdf',
  PLOT.tree)
ggsave(filename='${WD}/results/phylogeny/phylogeny/Tityra_${gene}.png',
  PLOT.tree)
"""

done <${WD}/results/phylogeny/concatenated/final_busco_ids.txt

