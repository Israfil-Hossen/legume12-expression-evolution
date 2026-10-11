#!/usr/bin/env bash
# Step 4: submit one quant job per biological replicate whose index and downloads are ready.
# An index counts as ready only when its PBS log exists and contains "done" (the log is written at job end).
# Re-runnable. Optional argument: maximum number of jobs to submit this time.
P=/home/liuchuanfa/israfil/legume2
cd $P; mkdir -p jobs/quant logs/quant
MAX=${1:-9999}; k=0; sub=0; held=0; wait_idx=0; wait_dl=0; done_n=0; queued=0
for br in $(awk -F'\t' 'NR>1{print $4}' metadata/step2_runs.tsv | sort -u); do
  k=$((k+1)); sp=${br%%_*}
  if [ -f quant/$br/.done ]; then done_n=$((done_n+1)); continue; fi
  if [ -f jobs/quant/$br.submitted ]; then queued=$((queued+1)); continue; fi
  if grep -qs '^done' logs/index_$sp.log; then :; else wait_idx=$((wait_idx+1)); continue; fi
  ready=1; for r in $(awk -F'\t' -v b="$br" 'NR>1 && $4==b{print $5}' metadata/step2_runs.tsv); do [ -f fastq/$r/.done ] || ready=0; done
  if [ $ready = 1 ]; then :; else wait_dl=$((wait_dl+1)); continue; fi
  if [ $sub -ge $MAX ]; then held=$((held+1)); continue; fi
  qsub -N q$(printf %03d $k)_$sp -v BR=$br,SP=$sp -o $P/logs/quant/$br.pbs.out quant_one.sh > /dev/null && touch jobs/quant/$br.submitted && sub=$((sub+1)) && echo "submitted $br"
done
echo "of 109 samples: done $done_n | in queue/running $queued | submitted now $sub | ready but held back $held | waiting for index $wait_idx | waiting for download $wait_dl"
