#!/bin/bash

###############################################################################
# Validation Script for Graph-Based Phylogeny Pipeline
# Checks directory structure and file integrity
###############################################################################

WD=/media/inter/mkapun/projects/Tityra

echo "==================================================================="
echo "Validating Graph-Based Phylogeny Pipeline"
echo "==================================================================="

# Check Step 9 outputs
echo ""
echo "Step 9: Graph Alignment Outputs"
echo "--------------------------------"

if [ -d "${WD}/results/graph_alignment" ]; then
    echo "✓ Graph alignment directory exists"
    
    NUM_GRAPHS=$(ls ${WD}/results/graph_alignment/graphs/*.xg 2>/dev/null | wc -l)
    echo "  - Variation graphs: ${NUM_GRAPHS}"
    
    NUM_CONSENSUS=$(ls ${WD}/results/graph_alignment/consensus/*_Tityra.fasta 2>/dev/null | wc -l)
    echo "  - Consensus sequences: ${NUM_CONSENSUS}"
    
    # Check if consensus sequences are non-empty
    if [ ${NUM_CONSENSUS} -gt 0 ]; then
        echo "  - Checking consensus sequence sizes..."
        EMPTY_COUNT=0
        for f in ${WD}/results/graph_alignment/consensus/*_Tityra.fasta; do
            SIZE=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f")
            if [ ${SIZE} -lt 50 ]; then
                EMPTY_COUNT=$((EMPTY_COUNT + 1))
            fi
        done
        if [ ${EMPTY_COUNT} -eq 0 ]; then
            echo "  ✓ All consensus sequences are non-empty"
        else
            echo "  ✗ Warning: ${EMPTY_COUNT} consensus sequences are suspiciously small"
        fi
    fi
else
    echo "✗ Graph alignment directory not found"
    echo "  Run Step 9 (graph_alignment_busco.sh) first"
fi

# Check graph-based BUSCO proteins
echo ""
echo "Graph-Based BUSCO Proteins"
echo "--------------------------"

if [ -d "${WD}/results/phylogeny_graph/busco_aa" ]; then
    echo "✓ Graph-based BUSCO directory exists"
    
    NUM_PROTEINS=$(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | wc -l)
    echo "  - Protein sequences: ${NUM_PROTEINS}"
    
    if [ ${NUM_PROTEINS} -gt 0 ]; then
        # Check for stop codons in proteins
        echo "  - Checking for premature stop codons..."
        STOP_COUNT=0
        for f in ${WD}/results/phylogeny_graph/busco_aa/*.faa; do
            # Count internal stop codons (not at the end)
            STOPS=$(grep -v "^>" "$f" | tr -d '\n' | sed 's/\*$//' | grep -o '\*' | wc -l)
            if [ ${STOPS} -gt 0 ]; then
                STOP_COUNT=$((STOP_COUNT + 1))
            fi
        done
        if [ ${STOP_COUNT} -eq 0 ]; then
            echo "  ✓ No premature stop codons detected"
        else
            echo "  ✗ Warning: ${STOP_COUNT} sequences have internal stop codons"
        fi
        
        # Show example sequence lengths
        echo "  - Example protein lengths:"
        for f in $(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa | head -5); do
            LEN=$(grep -v "^>" "$f" | tr -d '\n' | wc -c)
            GENE=$(basename "$f" .faa)
            echo "    ${GENE}: ${LEN} aa"
        done
    fi
else
    echo "✗ Graph-based BUSCO directory not found"
    echo "  This is created by graph_alignment_busco.sh Step 7"
fi

# Check Step 10 integration
echo ""
echo "Step 10: Phylogeny Integration"
echo "-------------------------------"

if [ -f "${WD}/results/phylogeny/data/genomes.names" ]; then
    echo "✓ Genome names file exists"
    
    if grep -q "Tityra_leucura_graph" "${WD}/results/phylogeny/data/genomes.names"; then
        echo "  ✓ Graph-based sequences are included"
    else
        echo "  ✗ Graph-based sequences not included"
        echo "    Run Step 10 after Step 9 to integrate"
    fi
    
    echo "  - Genomes in phylogeny:"
    cat "${WD}/results/phylogeny/data/genomes.names" | sed 's/^/    /'
else
    echo "✗ Phylogeny not yet built"
    echo "  Run Step 10 of main.sh"
fi

# Check if graph-based BUSCO structure exists for phylogeny
if [ -d "${WD}/results/phylogeny/BUSCO/Tityra_leucura_graph" ]; then
    echo "✓ Graph-based BUSCO structure created for phylogeny"
    
    NUM_GENES=$(ls ${WD}/results/phylogeny/BUSCO/Tityra_leucura_graph/run_aves_odb10/busco_sequences/single_copy_busco_sequences/*.faa 2>/dev/null | wc -l)
    echo "  - Individual gene files: ${NUM_GENES}"
else
    echo "✗ Graph-based BUSCO structure not created"
    echo "  This is created by main.sh Step 10 when graph sequences exist"
fi

# Check final phylogenetic tree
echo ""
echo "Final Phylogenetic Trees"
echo "------------------------"

if [ -f "${WD}/results/phylogeny/phylogeny/RAxML_bipartitions.FINAL" ]; then
    echo "✓ Combined phylogenetic tree exists"
    
    if grep -q "Tityra_leucura_graph" "${WD}/results/phylogeny/phylogeny/RAxML_bipartitions.FINAL"; then
        echo "  ✓ Graph-based sample included in tree"
    else
        echo "  ✗ Graph-based sample not in tree"
    fi
    
    echo "  - Taxa in tree:"
    grep "(" "${WD}/results/phylogeny/phylogeny/RAxML_bipartitions.FINAL" | \
        grep -o "[A-Za-z_]*:" | sed 's/:$//' | sort -u | sed 's/^/    /'
else
    echo "✗ Combined phylogenetic tree not found"
fi

# Check individual gene phylogenies
echo ""
echo "Individual Gene Phylogenies"
echo "---------------------------"

if [ -d "${WD}/results/phylogeny_graph/phylogeny" ]; then
    NUM_GENE_TREES=$(find ${WD}/results/phylogeny_graph/phylogeny -name "RAxML_bipartitions.FINAL" 2>/dev/null | wc -l)
    echo "✓ Graph-based gene phylogenies directory exists"
    echo "  - Individual gene trees: ${NUM_GENE_TREES}"
    
    if [ ${NUM_GENE_TREES} -gt 0 ]; then
        echo "  - Example genes with trees:"
        find ${WD}/results/phylogeny_graph/phylogeny -name "RAxML_bipartitions.FINAL" | head -5 | \
            xargs -I {} dirname {} | xargs -I {} basename {} | sed 's/^/    /'
    fi
else
    echo "✗ Graph-based gene phylogenies not found"
    echo "  Run BUSCO_phylogeny_graph_indgenes.sh to create"
fi

# Summary
echo ""
echo "==================================================================="
echo "Validation Summary"
echo "==================================================================="

ISSUES=0

# Critical checks
if [ ! -d "${WD}/results/phylogeny_graph/busco_aa" ]; then
    echo "❌ CRITICAL: Graph-based BUSCO proteins missing (Step 9)"
    ISSUES=$((ISSUES + 1))
elif [ $(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | wc -l) -eq 0 ]; then
    echo "❌ CRITICAL: No protein sequences found"
    ISSUES=$((ISSUES + 1))
fi

if [ -f "${WD}/results/phylogeny/data/genomes.names" ]; then
    if ! grep -q "Tityra_leucura_graph" "${WD}/results/phylogeny/data/genomes.names"; then
        echo "⚠️  WARNING: Graph sequences not integrated into phylogeny"
        echo "   Re-run Step 10 after Step 9 completes"
        ISSUES=$((ISSUES + 1))
    fi
fi

if [ ${ISSUES} -eq 0 ]; then
    echo "✅ All checks passed!"
    echo ""
    echo "Pipeline Status:"
    echo "  ✓ Graph alignment complete"
    echo "  ✓ Protein sequences generated"
    echo "  ✓ Integration successful"
    echo ""
    echo "You can now proceed with phylogenetic analysis."
else
    echo "⚠️  Found ${ISSUES} issue(s)"
    echo ""
    echo "Action Required:"
    echo "  1. Complete Step 9 (graph_alignment_busco.sh)"
    echo "  2. Run/Re-run Step 10 (main.sh phylogeny section)"
    echo "  3. Build individual gene trees (BUSCO_phylogeny_graph_indgenes.sh)"
fi

echo "==================================================================="
