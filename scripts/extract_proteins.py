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
