#!/bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16
#SBATCH --job-name=star_align
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --error star_align_%j.err
#SBATCH --output star_align_%j.out

# Activate the conda environment that contains STAR
conda activate angsd

output_dir="/athena/angsd/scratch/jif4007/project/alignment"
fastq_path="/athena/cayuga_0019/scratch/jif4007/project/fastq"
sra_ids=("SRR27351934" "SRR27351933" "SRR27351932" "SRR27351931" "SRR27351918" "SRR27351917")
STAR_index_path="/athena/angsd/scratch/jif4007/project/referenceGenomes/GRCh38_STARindex"

for sample_id in "${sra_ids[@]}"; do
  STAR --runMode alignReads \
       --readFilesCommand zcat \
       --runThreadN 16 \
       --genomeDir "$STAR_index_path" \
       --readFilesIn "${fastq_path}/${sample_id}_1_val_1.fq.gz" "${fastq_path}/${sample_id}_2_val_2.fq.gz" \
       --outFileNamePrefix "${output_dir}/${sample_id}.star." \
       --outSAMtype BAM SortedByCoordinate
  samtools index "${output_dir}/${sample_id}.star.Aligned.sortedByCoord.out.bam"
done
