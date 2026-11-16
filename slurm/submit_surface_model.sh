#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=log.surface_model_%j.out
module load apps/anaconda3/2024.06
python3 11_adsorption.py analyze > /dev/null
python3 12_surface_model.py
