#!/bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=bamqc
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --error qualimap_%j.err
#SBATCH --output qualimap_%j.out

output_dir="/athena/angsd/scratch/jif4007/project/alignment"
sra_ids=("SRR27351934" "SRR27351933" "SRR27351932" "SRR27351931" "SRR27351918" "SRR27351917")

cd "$output_dir"

conda activate qualimap

for sample_id in "${sra_ids[@]}"; do
  qualimap bamqc -pe -bam "${sample_id}.star.Aligned.sortedByCoord.out.bam" \
                 -outdir "qualimap_${sample_id}" \
                 -outformat PDF
done