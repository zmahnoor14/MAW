#!/bin/sh
#SBATCH --job-name=maw_r_run1
#SBATCH --mail-user=mahnoor.zulfiqar@uni-jena.de
#SBATCH --mail-type=ALL
#SBATCH --output=maw_r_run1.txt
#SBATCH --error=maw_r_run1.txt
#SBATCH --partition=s_standard,b_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=0

module load tools/singularity/3.7.0


srun singularity run --bind ./ docker://zmahnoor/maw-r@sha256:a47d3ecc4ac4489482229c00950733489e0340534e81932dc7121758a720fe6f 
