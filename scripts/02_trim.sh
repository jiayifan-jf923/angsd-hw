#!/bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16
#SBATCH --job-name=trim
#SBATCH --time=08:00:00
#SBATCH --mem=64G

conda activate trim-galore

# Navigate to the fastq directory
fastq_path=/athena/cayuga_0019/scratch/jif4007/project/fastq
cd "$fastq_path"

# Define SRA IDs
sra_ids=("SRR27351934" "SRR27351933" "SRR27351932" "SRR27351931" "SRR27351918" "SRR27351917")

for id in "${sra_ids[@]}"; do
  trim_galore --illumina --fastqc --paired "$id"_1.fastq.gz "$id"_2.fastq.gz
done