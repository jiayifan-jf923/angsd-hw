#!/bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16
#SBATCH --job-name=download
#SBATCH --time=08:00:00
#SBATCH --mem=64G

# Activate the ANGSD environment
conda activate angsd

# Create and navigate to the fastq directory
fastq_path=/athena/cayuga_0019/scratch/jif4007/project/fastq
mkdir -p "$fastq_path"
cd "$fastq_path"

# Define SRA IDs
sra_ids=("SRR27351934" "SRR27351933" "SRR27351932" "SRR27351931" "SRR27351918" "SRR27351917")

# Loop through each SRA ID to download and convert to FASTQ
for id in "${sra_ids[@]}"; do
  prefetch "$id"
  fastq-dump --split-files "$id" # paired-end, so I specified --split-files here
  bgzip -@ 16 "$id"_1.fastq
  bgzip -@ 16 "$id"_2.fastq
  fastqc "$id"_1.fastq.gz --extract
  fastqc "$id"_2.fastq.gz --extract
done

conda activate multiqc
multiqc . --filename multiqc_trim_galore_report.html --outdir multiqc_report
