#!/bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=count
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --error featurecounts_%j.err
#SBATCH --output featurecounts_%j.out

output_dir="/athena/angsd/scratch/jif4007/project/gene_counts"
sra_ids=("SRR27351934" "SRR27351933" "SRR27351932" "SRR27351931" "SRR27351918" "SRR27351917")
bam_dir="/athena/angsd/scratch/jif4007/project/alignment"  # adjust this if your BAM files are elsewhere

mkdir -p "$output_dir"
cd "$output_dir"

conda activate angsd

# Build list of BAM files
bam_files=()
for sample_id in "${sra_ids[@]}"; do
  bam_files+=("${bam_dir}/${sample_id}.star.Aligned.sortedByCoord.out.bam")
done

# Run featureCounts once for all
featureCounts -T 16 -p --countReadPairs -t exon -g gene_id\
  -a /athena/cayuga_0019/scratch/jif4007/project/referenceGenomes/gencode.v45.basic.annotation.gtf \
  -o "${output_dir}/fc_output" "${bam_files[@]}"
