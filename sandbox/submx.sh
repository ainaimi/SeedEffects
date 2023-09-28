#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=30-00:00:00
#SBATCH --job-name=seed_effects
#SBATCH --mem=190g
#SBATCH --partition=naimi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anaimi@emory.edu

module purge
module load R

Rscript --no-save --no-restore --verbose ./xgboost_test.R 128 $i > xg_run_$i.Rout 2>&1