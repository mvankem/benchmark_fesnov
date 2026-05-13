#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --job-name=fesnov
#SBATCH -o run.out
#SBATCH -e run.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate benchmark-fesnov

snakemake --cores 32 "$@"
