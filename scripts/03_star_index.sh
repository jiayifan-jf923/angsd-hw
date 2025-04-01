#!/bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=star_index
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --error star_index_%j.err
#SBATCH --output star_index_%j.out

# Activate the conda environment that contains STAR
conda activate angsd

# Run STAR genomeGenerate command
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir /athena/angsd/scratch/jif4007/project/referenceGenomes/GRCh38_STARindex \
     --genomeFastaFiles /athena/angsd/scratch/las2017/angsd/referenceGenomes/GRCh38.primary_assembly.v45.genome.fa \
     --sjdbGTFfile /athena/cayuga_0019/scratch/jif4007/project/referenceGenomes/gencode.v45.basic.annotation.gtf \
     --sjdbOverhang 150