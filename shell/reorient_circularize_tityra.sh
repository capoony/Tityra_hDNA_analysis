#!/bin/bash
###############################################################################
# Reorient and Circularize Tityra Mitochondrial Genome
# Author: Generated for Tityra project
# Description: Fixes reverse-complement orientation and sets proper start position
###############################################################################

WD="/media/inter/mkapun/projects/Tityra"

# Create output directories
mkdir -p "${WD}/results/tityra_reorientation"
mkdir -p "${WD}/results/tityra_reorientation/logs"

echo "======================================================================="
echo "Tityra Mitochondrial Genome Reorientation"
echo "======================================================================="

# Define environment path
ENV_PATH="${WD}/scripts/programs/mito_tools"

# Check if conda environment exists, create if not
if [ ! -d "${ENV_PATH}" ]; then
    echo "Creating mito_tools conda environment in scripts/programs/..."
    
    mamba create -p "${ENV_PATH}" \
        python=3.11 \
        biopython \
        seqkit \
        minimap2 \
        samtools \
        blast \
        mummer4 \
        matplotlib \
        -c bioconda -c conda-forge -y
    
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create mito_tools environment"
        exit 1
    fi
fi

# Activate environment
source /opt/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate "${ENV_PATH}"

# Input files
TITYRA_ASSEMBLY="${WD}/results/mitofinder_assembly/Tityra_leucura/Tityra_leucura/Tityra_leucura_MitoFinder_mitfi_Final_Results/Tityra_leucura_mtDNA_contig_1.fasta"
TITYRA_GB_ORIG="${WD}/results/mitofinder_assembly/Tityra_leucura/Tityra_leucura/Tityra_leucura_MitoFinder_mitfi_Final_Results/Tityra_leucura_mtDNA_contig_1.gb"
TITYRA_GB="${WD}/data/refseq_genbank/Tityra_leucura.gb"
REFERENCE_GB="${WD}/data/refseq_genbank/Pachyramphus_minor.gb"

# Extract reference FASTA from GenBank if needed
REFERENCE="${WD}/results/tityra_reorientation/Pachyramphus_minor_ref.fasta"

# Use the start region of Pachyramphus as the seed for proper start position
GENE_SEED="${WD}/results/tityra_reorientation/start_region_seed.fasta"

# Check if input files exist
if [ ! -f "${TITYRA_ASSEMBLY}" ]; then
    echo "Error: Tityra assembly not found: ${TITYRA_ASSEMBLY}"
    echo "Please run mitofinder first (Section 6)"
    conda deactivate
    exit 1
fi

# Copy original Tityra GenBank file to refseq_genbank if not already there
if [ ! -f "${TITYRA_GB}" ] && [ -f "${TITYRA_GB_ORIG}" ]; then
    echo "Copying original Tityra GenBank file..."
    cp "${TITYRA_GB_ORIG}" "${TITYRA_GB}"
fi

if [ ! -f "${REFERENCE_GB}" ]; then
    echo "Error: Reference GenBank not found: ${REFERENCE_GB}"
    echo "Please run download_complete_mitogenomes.sh first (Section 11)"
    conda deactivate
    exit 1
fi

# Extract reference FASTA from GenBank
echo "Extracting reference FASTA from GenBank..."
python3 << PYEOF
from Bio import SeqIO
record = SeqIO.read("${REFERENCE_GB}", "genbank")
SeqIO.write(record, "${REFERENCE}", "fasta")
print(f"  Extracted: {len(record.seq)} bp")
PYEOF

echo ""
echo "Input files:"
echo "  Assembly:     ${TITYRA_ASSEMBLY}"
echo "  Reference GB: ${REFERENCE_GB}"
echo "  Reference FA: ${REFERENCE}"
echo "  GenBank:      ${TITYRA_GB}"

###############################################################################
# STEP 1: Extract start region from Pachyramphus as seed
###############################################################################
echo ""
echo "Step 1: Extracting start region from Pachyramphus_minor..."

python3 << 'PYEOF'
from Bio import SeqIO

# Extract first 1500 bp from Pachyramphus_minor as start seed
gb_file = "/media/inter/mkapun/projects/Tityra/data/refseq_genbank/Pachyramphus_minor.gb"

try:
    record = SeqIO.read(gb_file, "genbank")
    # Get first 1500 bp as start region seed
    start_seq = record.seq[:1500]
    
    with open("/media/inter/mkapun/projects/Tityra/results/tityra_reorientation/start_region_seed.fasta", "w") as f:
        f.write(f">start_region_seed\n")
        f.write(str(start_seq) + "\n")
    
    print(f"  Extracted start region ({len(start_seq)} bp) from Pachyramphus_minor")
except Exception as e:
    print(f"Error: {e}")
    import sys
    sys.exit(1)
PYEOF

if [ $? -ne 0 ]; then
    echo "Error: Failed to extract start region seed"
    conda deactivate
    exit 1
fi

###############################################################################
# STEP 2: Check orientation with minimap2
###############################################################################
echo ""
echo "Step 2: Checking orientation against reference..."

seqkit seq --line-width 0 "${TITYRA_ASSEMBLY}" \
    > "${WD}/results/tityra_reorientation/tityra_clean.fasta"

minimap2 -x asm5 -t 4 \
    "${REFERENCE}" \
    "${WD}/results/tityra_reorientation/tityra_clean.fasta" \
    > "${WD}/results/tityra_reorientation/minimap2_aln.paf" \
    2>"${WD}/results/tityra_reorientation/logs/minimap2.log"

# Detect dominant strand
STRAND=$(awk '{print $5}' "${WD}/results/tityra_reorientation/minimap2_aln.paf" | \
         sort | uniq -c | sort -rn | head -1 | awk '{print $2}')

echo "  Dominant strand: ${STRAND}"

if [[ "$STRAND" == "-" ]]; then
    echo "  Assembly is reverse complemented - fixing..."
    NEEDS_RC=1
else
    echo "  Orientation looks correct"
    NEEDS_RC=0
fi

###############################################################################
# STEP 3: Reverse complement if needed
###############################################################################
echo ""
echo "Step 3: Applying reverse complement if needed..."

if [[ "$NEEDS_RC" -eq 1 ]]; then
    seqkit seq -r -p \
        "${WD}/results/tityra_reorientation/tityra_clean.fasta" \
        > "${WD}/results/tityra_reorientation/tityra_rc.fasta"
    echo "  Reverse complemented"
    
    # Verify with minimap2
    minimap2 -x asm5 -t 4 \
        "${REFERENCE}" \
        "${WD}/results/tityra_reorientation/tityra_rc.fasta" \
        > "${WD}/results/tityra_reorientation/minimap2_aln_after_rc.paf" 2>/dev/null
    
    STRAND2=$(awk '{print $5}' "${WD}/results/tityra_reorientation/minimap2_aln_after_rc.paf" | \
              sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
    
    if [[ "$STRAND2" == "+" ]]; then
        echo "  Orientation confirmed correct after RC"
    fi
    
    ORIENTED_FA="${WD}/results/tityra_reorientation/tityra_rc.fasta"
else
    cp "${WD}/results/tityra_reorientation/tityra_clean.fasta" \
       "${WD}/results/tityra_reorientation/tityra_oriented.fasta"
    ORIENTED_FA="${WD}/results/tityra_reorientation/tityra_oriented.fasta"
fi

###############################################################################
# STEP 4: Find start position matching Pachyramphus
###############################################################################
echo ""
echo "Step 4: Finding start position to match Pachyramphus..."

# Build BLAST database
makeblastdb -in "${ORIENTED_FA}" -dbtype nucl \
    -out "${WD}/results/tityra_reorientation/blastdb" \
    -title "tityra_oriented" >/dev/null 2>&1

# BLAST start region seed
blastn -query "${GENE_SEED}" \
       -db "${WD}/results/tityra_reorientation/blastdb" \
       -outfmt "6 qseqid sseqid pident length sstart send bitscore" \
       -perc_identity 70 \
       -num_alignments 5 \
       -out "${WD}/results/tityra_reorientation/blast_start.tsv" 2>/dev/null

if [[ ! -s "${WD}/results/tityra_reorientation/blast_start.tsv" ]]; then
    echo "  Warning: No BLAST hits found with 70% identity, trying 60%..."
    blastn -query "${GENE_SEED}" \
           -db "${WD}/results/tityra_reorientation/blastdb" \
           -outfmt "6 qseqid sseqid pident length sstart send bitscore" \
           -perc_identity 60 \
           -num_alignments 5 \
           -out "${WD}/results/tityra_reorientation/blast_start.tsv" 2>/dev/null
fi

if [[ ! -s "${WD}/results/tityra_reorientation/blast_start.tsv" ]]; then
    echo "  Warning: Could not find matching start region - will not rotate genome"
    START_POS=1
else
    # Get start position from best hit
    START_POS=$(sort -k7 -rn "${WD}/results/tityra_reorientation/blast_start.tsv" | \
                head -1 | awk '{
                  s=$5; e=$6
                  if (s > e) { tmp=s; s=e; e=tmp }
                  print s
                }')
    
    echo "  Start region found at position: ${START_POS}"
    cat "${WD}/results/tityra_reorientation/blast_start.tsv"
fi

###############################################################################
# STEP 5: Rotate genome to match Pachyramphus start
###############################################################################
echo ""
echo "Step 5: Rotating genome to match Pachyramphus start position..."

if (( START_POS <= 1 )); then
    echo "  Start position is already at position 1 - no rotation needed"
    cp "${ORIENTED_FA}" "${WD}/results/tityra_reorientation/tityra_rotated.fasta"
else
    seqkit restart -i "${START_POS}" \
        "${ORIENTED_FA}" \
        > "${WD}/results/tityra_reorientation/tityra_rotated.fasta"
    echo "  Genome rotated to start at position ${START_POS}"
fi

# Update header
seqkit replace -p ".*" -r "Tityra_leucura_mtDNA_reoriented" \
    "${WD}/results/tityra_reorientation/tityra_rotated.fasta" \
    > "${WD}/results/tityra_reorientation/Tityra_leucura_final.fasta"

###############################################################################
# STEP 6: Update GenBank file with reoriented sequence
###############################################################################
echo ""
echo "Step 6: Creating updated GenBank file with corrected features..."

python3 << PYEOF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

wd = "/media/inter/mkapun/projects/Tityra"

# Load new sequence
new_seq_file = f"{wd}/results/tityra_reorientation/Tityra_leucura_final.fasta"
new_seq = SeqIO.read(new_seq_file, "fasta")

# Load original GenBank
gb_file_orig = f"{wd}/results/mitofinder_assembly/Tityra_leucura/Tityra_leucura/Tityra_leucura_MitoFinder_mitfi_Final_Results/Tityra_leucura_mtDNA_contig_1.gb"
gb_rec = SeqIO.read(gb_file_orig, "genbank")

# Check if we did reverse complement and rotation
needs_rc = ${NEEDS_RC}
start_pos = ${START_POS}
orig_length = len(gb_rec.seq)

print(f"Original length: {orig_length} bp")
print(f"Reverse complement: {needs_rc == 1}")
print(f"Rotation start: {start_pos}")

# Transform features based on operations performed
new_features = []

for feature in gb_rec.features:
    new_feature = feature
    
    # If we did reverse complement, update feature
    if needs_rc == 1:
        # For reverse complement:
        # - Swap strand (+ becomes -, - becomes +)
        # - Reverse position (pos becomes length - pos)
        
        old_start = int(feature.location.start)
        old_end = int(feature.location.end)
        old_strand = feature.location.strand
        
        # New positions (reverse complement)
        new_start = orig_length - old_end
        new_end = orig_length - old_start
        new_strand = -old_strand if old_strand else None
        
        # Create new feature location
        new_loc = FeatureLocation(new_start, new_end, strand=new_strand)
        new_feature = SeqFeature(
            location=new_loc,
            type=feature.type,
            qualifiers=feature.qualifiers.copy()
        )
    
    # If we rotated, adjust positions
    if start_pos > 1:
        old_start = int(new_feature.location.start)
        old_end = int(new_feature.location.end)
        strand = new_feature.location.strand
        
        # Rotate positions
        new_start = (old_start - start_pos + 1) % orig_length
        new_end = (old_end - start_pos + 1) % orig_length
        
        # Handle wrap-around features
        if new_start > new_end:
            # Feature wraps around origin - split it
            # For simplicity, we'll keep the larger portion
            if new_end < (orig_length - new_start):
                # Keep the end portion
                new_start = 0
            else:
                # Keep the start portion  
                new_end = orig_length
        
        new_loc = FeatureLocation(new_start, new_end, strand=strand)
        new_feature = SeqFeature(
            location=new_loc,
            type=new_feature.type,
            qualifiers=new_feature.qualifiers.copy()
        )
    
    new_features.append(new_feature)

# Create new record with updated sequence and transformed features
new_rec = SeqRecord(
    seq=new_seq.seq,
    id=gb_rec.id,
    name=gb_rec.name,
    description="Tityra leucura mitochondrion, complete genome (reoriented)",
    features=new_features,
    annotations=gb_rec.annotations.copy()
)

# Update topology to circular
new_rec.annotations["topology"] = "circular"

# Write updated GenBank - replace the original file
out_gb = f"{wd}/data/refseq_genbank/Tityra_leucura.gb"
out_gb_backup = f"{wd}/data/backup/Tityra_leucura_original.gb"

# Backup original if not already backed up
import os
os.makedirs(f"{wd}/data/backup", exist_ok=True)
if not os.path.exists(out_gb_backup) and os.path.exists(out_gb):
    import shutil
    shutil.copy2(out_gb, out_gb_backup)
    print(f"  Backed up original: {out_gb_backup}")

# Write reoriented version as the main file
SeqIO.write(new_rec, out_gb, "genbank")
print(f"  Updated GenBank: {out_gb}")
print(f"  Length: {len(new_rec.seq)} bp")
print(f"  Features: {len(new_rec.features)}")

# Also update the FASTA file in refseq_genbank
out_fa = f"{wd}/data/refseq_genbank/Tityra_leucura.fasta"
SeqIO.write(new_rec, out_fa, "fasta")
print(f"  Updated FASTA: {out_fa}")
PYEOF

###############################################################################
# STEP 7: Create diagnostic plots
###############################################################################
echo ""
echo "Step 7: Creating diagnostic visualizations..."

python3 << 'PYEOF'
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

wd = "/media/inter/mkapun/projects/Tityra"

# Read sequences
ref_file = f"{wd}/results/tityra_reorientation/Pachyramphus_minor_ref.fasta"
tityra_orig = f"{wd}/results/tityra_reorientation/tityra_clean.fasta"
tityra_final = f"{wd}/results/tityra_reorientation/Tityra_leucura_final.fasta"

ref_seq = SeqIO.read(ref_file, "fasta")
orig_seq = SeqIO.read(tityra_orig, "fasta")
final_seq = SeqIO.read(tityra_final, "fasta")

# Calculate GC content in sliding windows
def calc_gc_windows(seq, window=100):
    gc_content = []
    for i in range(0, len(seq) - window, window):
        subseq = seq[i:i+window]
        gc = (subseq.count('G') + subseq.count('C') + 
              subseq.count('g') + subseq.count('c')) / len(subseq)
        gc_content.append(gc * 100)
    return gc_content

print("  Calculating GC content profiles...")
ref_gc = calc_gc_windows(str(ref_seq.seq))
orig_gc = calc_gc_windows(str(orig_seq.seq))
final_gc = calc_gc_windows(str(final_seq.seq))

# Create comparison plot
fig, axes = plt.subplots(3, 1, figsize=(14, 10))

# Reference
axes[0].plot(ref_gc, color='blue', linewidth=0.5, alpha=0.7)
axes[0].set_title(f'Pachyramphus_minor (Reference) - {len(ref_seq.seq)} bp', fontsize=12, fontweight='bold')
axes[0].set_ylabel('GC%', fontsize=10)
axes[0].set_xlim(0, len(ref_gc))
axes[0].grid(True, alpha=0.3)

# Original Tityra
axes[1].plot(orig_gc, color='red', linewidth=0.5, alpha=0.7)
axes[1].set_title(f'Tityra_leucura (Original) - {len(orig_seq.seq)} bp', fontsize=12, fontweight='bold')
axes[1].set_ylabel('GC%', fontsize=10)
axes[1].set_xlim(0, len(orig_gc))
axes[1].grid(True, alpha=0.3)

# Reoriented Tityra
axes[2].plot(final_gc, color='green', linewidth=0.5, alpha=0.7)
axes[2].set_title(f'Tityra_leucura (Reoriented) - {len(final_seq.seq)} bp', fontsize=12, fontweight='bold')
axes[2].set_ylabel('GC%', fontsize=10)
axes[2].set_xlabel('Position (100 bp windows)', fontsize=10)
axes[2].set_xlim(0, len(final_gc))
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{wd}/results/tityra_reorientation/gc_content_comparison.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{wd}/results/tityra_reorientation/gc_content_comparison.png", dpi=300, bbox_inches='tight')
plt.close()

print("  Saved: gc_content_comparison.pdf/png")

# Create alignment comparison plot
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Read minimap2 PAF files
paf_before = f"{wd}/results/tityra_reorientation/minimap2_aln.paf"
paf_after_file = f"{wd}/results/tityra_reorientation/minimap2_aln_after_rc.paf"

def read_paf(paf_file):
    alignments = []
    try:
        with open(paf_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 11:
                    alignments.append({
                        'qstart': int(fields[2]),
                        'qend': int(fields[3]),
                        'strand': fields[4],
                        'tstart': int(fields[7]),
                        'tend': int(fields[8]),
                        'mapq': int(fields[11])
                    })
    except:
        pass
    return alignments

aln_before = read_paf(paf_before)
aln_after = read_paf(paf_after_file) if aln_after_file else []

# Plot alignments
if aln_before:
    for aln in aln_before[:10]:  # Plot top 10
        color = 'red' if aln['strand'] == '-' else 'blue'
        ax.plot([aln['qstart'], aln['qend']], 
                [aln['tstart'], aln['tend']], 
                color=color, alpha=0.5, linewidth=2, label='Original' if color == 'red' else '')

if aln_after:
    for aln in aln_after[:10]:
        color = 'green' if aln['strand'] == '+' else 'orange'
        ax.plot([aln['qstart'], aln['qend']], 
                [aln['tstart'], aln['tend']], 
                color=color, alpha=0.5, linewidth=1, linestyle='--')

ax.set_xlabel('Tityra position (bp)', fontsize=12)
ax.set_ylabel('Pachyramphus position (bp)', fontsize=12)
ax.set_title('Alignment: Tityra vs Pachyramphus\nRed=Original(RC), Green=Corrected', 
             fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{wd}/results/tityra_reorientation/alignment_comparison.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{wd}/results/tityra_reorientation/alignment_comparison.png", dpi=300, bbox_inches='tight')
plt.close()

print("  Saved: alignment_comparison.pdf/png")

print("\n  Diagnostic plots created successfully!")
PYEOF

echo ""
echo "Diagnostic plots saved in: results/tityra_reorientation/"

conda deactivate

echo ""
echo "======================================================================="
echo "Reorientation Complete!"
echo "======================================================================="
echo "Output files:"
echo "  - Final FASTA:  results/tityra_reorientation/Tityra_leucura_final.fasta"
echo "  - Updated GB:   data/refseq_genbank/Tityra_leucura_reoriented.gb"
echo "  - Updated FASTA: data/refseq_genbank/Tityra_leucura.fasta"
echo ""
echo "The genome has been:"
if [[ "$NEEDS_RC" -eq 1 ]]; then
    echo "  ✓ Reverse complemented"
else
    echo "  ✓ Already in correct orientation"
fi
if (( START_POS > 1 )); then
    echo "  ✓ Rotated to match Pachyramphus start (position ${START_POS})"
else
    echo "  ✓ Already starts at position 1"
fi
echo ""
echo "You can now re-run the alignment scripts to use the corrected genome."
