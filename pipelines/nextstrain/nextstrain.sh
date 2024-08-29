#!/bin/bash
#SBATCH --job-name=nextstrain
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=nextstrain.out
#SBATCH --error=nextstrain.err

module load miniconda
conda activate nextstrain

snakemake --cores all export --rerun-incomplete