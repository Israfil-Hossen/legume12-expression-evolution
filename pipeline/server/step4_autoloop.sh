#!/usr/bin/env bash
# Every 30 min: submit Step 4 quant jobs for samples whose downloads have finished.
P=/home/liuchuanfa/israfil/legume2
cd $P
while true; do
  out=$(./step4_submit_quant.sh)
  echo "=== $(date '+%F %H:%M') ==="; echo "$out"
  if echo "$out" | grep -q "waiting for index 0 | waiting for download 0"; then echo "all samples submitted - loop finished"; break; fi
  sleep 1800
done
