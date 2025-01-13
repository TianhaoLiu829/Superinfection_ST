#!/bin/bash
#SBATCH --job-name rctd_lung_github
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --time=24:00:00
#SBATCH --mem=600GB

module load gcc/12.2.0 r/4.3.0

R CMD BATCH --no-restore /ix1/wchen/liutianhao/result/lung_ST/alcorn/script/original_analysis/RCTD.R out_RCTD.txt
