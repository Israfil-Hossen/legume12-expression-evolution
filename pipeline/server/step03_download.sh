#!/usr/bin/env bash
# Step 3: download FASTQ for every run in metadata/step2_runs.tsv (login node). Resumable; md5-checked.
P=/home/liuchuanfa/israfil/legume2
cd $P
dl_one() {
  run=$1; urls=$2; md5s=$3; d=$P/fastq/$run; mkdir -p $d
  [ -f $d/.done ] && { echo "[$run] already done"; return 0; }
  IFS=';' read -ra U <<< "$urls"; IFS=';' read -ra M <<< "$md5s"; ok=1
  for i in "${!U[@]}"; do
    f=$d/$(basename ${U[$i]})
    wget -q -c --tries=10 -O $f "https://${U[$i]}" || ok=0
    echo "${M[$i]}  $f" | md5sum -c --quiet 2>/dev/null || { ok=0; rm -f $f; }
  done
  if [ $ok = 1 ]; then touch $d/.done; echo "[$run] OK"; else echo "[$run] FAILED"; fi
}
export -f dl_one; export P
tail -n +2 metadata/step2_runs.tsv | awk -F'\t' '{print $5"\t"$11"\t"$12}' | parallel --will-cite -j 4 --colsep '\t' dl_one {1} {2} {3}
echo "=== summary: $(ls fastq/*/.done 2>/dev/null | wc -l) of $(($(wc -l < metadata/step2_runs.tsv)-1)) runs downloaded ==="
