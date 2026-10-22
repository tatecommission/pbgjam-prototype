#!/bin/bash
# usage: bash launch_1km.sh <first pred> <last pred>     e.g.  bash launch_1km.sh 1 8
# For each prediction: submit a 30-task array (one tile per task, 24 h limit), then an assemble job
# that starts only if every tile succeeded.  Set W2W_NSIM before running to override nsim.
set -euo pipefail
NT=${W2W_NTILES:-30}
LOG=/work/tvc11/PBGJAM/out/pred_logs_1km
mkdir -p "$LOG"
cd /work/tvc11/PBGJAM/v4_gjam_workflow
for P in $(seq "$1" "$2"); do
  JID=$(sbatch --parsable --job-name=w1km_$P --array=1-$NT \
        --output=$LOG/%x_%A_%a.out --error=$LOG/%x_%A_%a.err submit_1km.sh "$P")
  AID=$(sbatch --parsable --job-name=a1km_$P --dependency=afterok:$JID \
        --output=$LOG/%x_%j.out --error=$LOG/%x_%j.err submit_1km_assemble.sh "$P")
  echo "pred $P: tiles job $JID, assemble job $AID"
done
