#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --time=1:00:00 #Based on initial test
#SBATCH --nodes=1 #MulticoreParam needs to run on one node
#SBATCH --cpus-per-task=1 #CPUS for parellization
#SBATCH --mem=256G #To load seurat objects
#SBATCH --job-name=singleR_to_loupeR
#SBATCH --output=%x-%j.out #Log file for errors or output
#SBATCH --mail-user=ssaab@uoguelph.ca
#SBATCH --mail-type=ALL

#load modules
module load StdEnv/2023 udunits/2.2.28 r/4.4.0 gcc/12.3 hdf5/1.14

#Run my Rscript
Rscript /home/ssaab/links/scratch/R_scripts/4.making_loupe_objects.R

