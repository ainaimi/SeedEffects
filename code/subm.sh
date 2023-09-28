#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=30-00:00:00
#SBATCH --job-name=seed_effects
#SBATCH --mem=120g
#SBATCH --partition=naimi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anaimi@emory.edu

module purge
module load R

for i in {1..5}
do
Rscript --no-save --no-restore --verbose ./analysis_numom.R 32 $i > seed_run_$i.Rout 2>&1

Rscript --no-save --no-restore --verbose ./analysis_numom.R 64 $i > seed_run_$i.Rout 2>&1

Rscript --no-save --no-restore --verbose ./analysis_numom.R 96 $i > seed_run_$i.Rout 2>&1
 
Rscript --no-save --no-restore --verbose ./analysis_numom.R 128 $i > seed_run_$i.Rout 2>&1
done

Rscript --no-save --no-restore --verbose ./run_time.R > run_time.Rout 2>&1