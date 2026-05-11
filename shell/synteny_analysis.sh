#!/bin/bash
###############################################################################
# Mitochondrial Genome Synteny Analysis with Clinker
# Author: Generated for Tityra project
# Description: Visualize gene synteny across mitochondrial genomes
###############################################################################

WD="/media/inter/mkapun/projects/Tityra"

# Create output directories
mkdir -p "${WD}/results/synteny_analysis"
mkdir -p "${WD}/results/synteny_analysis/proteins"

echo "======================================================================="
echo "Mitochondrial Genome Synteny Analysis"
echo "======================================================================="

# Define environment path
ENV_PATH="${WD}/scripts/programs/synteny"

# Check if conda environment exists
if [ ! -d "${ENV_PATH}" ]; then
    echo "Creating synteny conda environment in scripts/programs/..."
    
    mamba create -p "${ENV_PATH}" \
        python=3.11 \
        biopython \
        mafft \
        muscle \
        fasttree \
        clinker-py \
        clustalo \
        blast \
        seqkit \
        emboss \
        -c bioconda -c conda-forge -y
    
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create synteny environment"
        exit 1
    fi
fi

echo ""
echo "Extracting protein sequences from GenBank files..."

# Create Python script to extract proteins
cat > "${WD}/scripts/extract_proteins.py" << 'PYEOF'
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import glob

wd = "/media/inter/mkapun/projects/Tityra"
gb_dir = f"{wd}/data/refseq_genbank"
out_dir = f"{wd}/results/synteny_analysis/proteins"

print("Processing GenBank files...")
for gbk_file in sorted(glob.glob(f"{gb_dir}/*.gb")):
    species = os.path.basename(gbk_file).replace(".gb", "")
    print(f"  {species}")
    
    out_records = []
    
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "translation" in feature.qualifiers:
                    protein = feature.qualifiers["translation"][0]
                    
                    gene = feature.qualifiers.get(
                        "gene",
                        ["unknown"]
                    )[0]
                    
                    locus = feature.qualifiers.get(
                        "locus_tag",
                        [f"{species}_{gene}"]
                    )[0]
                    
                    out_records.append(
                        SeqRecord(
                            Seq(protein),
                            id=locus,
                            description=gene
                        )
                    )
    
    if out_records:
        out_file = f"{out_dir}/{species}.faa"
        SeqIO.write(out_records, out_file, "fasta")
        print(f"    Wrote {len(out_records)} proteins to {os.path.basename(out_file)}")
    else:
        print(f"    Warning: No proteins found")

print("\nProtein extraction complete!")
PYEOF

# Activate environment and run extraction
source /opt/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate "${ENV_PATH}"

python "${WD}/scripts/extract_proteins.py"

if [ $? -ne 0 ]; then
    echo "Error: Protein extraction failed"
    conda deactivate
    exit 1
fi

echo ""
echo "Running clinker for synteny analysis..."
echo "Input GenBank files: ${WD}/data/refseq_genbank/*.gb"
echo ""

# Run clinker on all GenBank files
cd "${WD}/results/synteny_analysis"

clinker "${WD}/data/refseq_genbank"/*.gb \
    -p clinker_plot.html \
    -o clinker_output.txt \
    -na \
    --jobs 4

if [ $? -ne 0 ]; then
    echo "Error: Clinker analysis failed"
    conda deactivate
    exit 1
fi

conda deactivate

echo ""
echo "======================================================================="
echo "Synteny Analysis Complete!"
echo "======================================================================="
echo "Output files:"
echo "  - HTML plot: results/synteny_analysis/clinker_plot.html"
echo "  - Alignment: results/synteny_analysis/clinker_output.txt"
echo "  - Proteins:  results/synteny_analysis/proteins/"
echo ""
echo "Open clinker_plot.html in a web browser to view the interactive synteny plot"
