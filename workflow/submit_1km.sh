#!/bin/bash
#SBATCH --account=qiulab
#SBATCH --partition=common
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tvc11@duke.edu
# usage (array index = tile):  sbatch --array=1-30 --job-name=w1_5 submit_1km.sh <pred 1..26>
set -euo pipefail
cd /work/tvc11/PBGJAM/v4_gjam_workflow
module load R GDAL/3.2.1-rhel8 GEOS/3.9.1-rhel8
export W2W_NTILES=${W2W_NTILES:-30}
Rscript R/07_predict_1km_tiled.R "$1" "${SLURM_ARRAY_TASK_ID}"
