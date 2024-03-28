#!/bin/bash
#SBATCH --job-name=arrayJob_cvtobitnet_6jan
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --output=arrayJob_cvtobitnet_6jan_%A_%a.out
#SBATCH --error=arrayJob_cvtobitnet_6jan_%A_%a.err
#SBATCH --array=1-24

cd /home/mzeeb/data/tobitnet

module load mamba

source activate tobitnet

./tobitnet_main_array.sh $SLURM_ARRAY_TASK_ID

