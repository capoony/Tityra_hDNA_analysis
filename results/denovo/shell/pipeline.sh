# Copying data
# Mi Jul 16 16:34:58 CEST 2025
## FWD Illumina data copied
# Mi Jul 16 16:35:15 CEST 2025
########################

# Start raw QC
## ... of Illumina data
# Mi Jul 16 16:35:15 CEST 2025
bash FullPipeline_exp/fastqc.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra /media/inter/pipelines/AutDeNovo 2 200 no no
# ########################

# Estimation of genomesize
# Mi Jul 16 16:43:25 CEST 2025
bash FullPipeline_exp/genomesize.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra ILL no /media/inter/pipelines/AutDeNovo 150 200 1000 no no
# ########################

## Starting denovo assembly
# Mi Jul 16 16:51:31 CEST 2025
bash FullPipeline_exp/denovo.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra ILL no /media/inter/pipelines/AutDeNovo 8 1000 no no
# Assembly of ILL with Spades
# Mi Jul 16 16:51:31 CEST 2025
# ########################

## Starting polishing with Racon
# Do Jul 17 03:27:33 CEST 2025
bash FullPipeline_exp/racon.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra ILL /media/inter/pipelines/AutDeNovo 150 1000 4 no no no
# ########################

# Starting assembly QC with Quast
# Do Jul 17 05:42:22 CEST 2025
bash FullPipeline_exp/quast.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra ILL /media/inter/pipelines/AutDeNovo 150 200 no no
# ########################

# Starting assembly QC with BUSCO
# Do Jul 17 05:42:31 CEST 2025
bash FullPipeline_exp/busco.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra vertebrata_odb10 ILL /media/inter/pipelines/AutDeNovo 150 200 no no
# ########################

# Mapping reads against reference
# Do Jul 17 05:52:54 CEST 2025
bash FullPipeline_exp/mapping.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra ILL no /media/inter/pipelines/AutDeNovo 150 200 no no
# ########################

# BLASTing genome against the nt database
# Do Jul 17 06:14:06 CEST 2025
bash FullPipeline_exp/blast.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra ILL /media/inter/pipelines/AutDeNovo 150 200 /media/scratch/NCBI_nt_DB_210714/nt no no
# ########################

# Summarize with Blobtools
# Do Jul 17 06:33:49 CEST 2025
bash FullPipeline_exp/blobtools.sh /media/inter/mkapun/projects/Tityra/results/denovo Tityra vertebrata_odb10 ILL /media/inter/pipelines/AutDeNovo 150 200 no no  $
# ########################

# Anlayses done!!
# Now copying results to output folder and writing commands for HTML output
# check /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_HTML_outputs.sh for more details
# Do Jul 17 06:33:49 CEST 2025
########################

