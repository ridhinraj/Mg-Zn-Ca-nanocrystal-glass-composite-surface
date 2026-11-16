#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=log.oh_report_%j.out
module load apps/anaconda3/2024.06
python3 21_oh.py report
