#!/bin/bash
#SBATCH --account=qiulab
#SBATCH --partition=common
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tvc11@duke.edu
# usage: sbatch --job-name=a1_5 submit_1km_assemble.sh <pred 1..26>
set -euo pipefail
cd /work/tvc11/PBGJAM/v4_gjam_workflow
module load R GDAL/3.2.1-rhel8 GEOS/3.9.1-rhel8
export W2W_NTILES=${W2W_NTILES:-30}
Rscript R/07c_assemble_1km.R "$1"
