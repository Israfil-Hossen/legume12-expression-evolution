#!/usr/bin/env bash
# Step 3b: robust FASTQ download with aria2c (multi-connection, resumable, md5-verified).
P=/home/liuchuanfa/israfil/legume2
cd $P
get_file() {
  url=$1; md5=$2; d=$3; b=$(basename $url); f=$d/$b
  for attempt in $(seq 1 30); do
    aria2c -q -c -x 8 -s 8 -k 1M --file-allocation=none --max-tries=5 --retry-wait=15 --timeout=60 --connect-timeout=30 -d $d -o $b "https://$url"
    if [ -f "$f" ]; then
      if [ -f "$f.aria2" ]; then :; else
        if [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then return 0; fi
        rm -f "$f"
      fi
    fi
    sleep 20
  done
  return 1
}
dl_one() {
  run=$1; urls=$2; md5s=$3; d=$P/fastq/$run; mkdir -p $d
  if [ -f $d/.done ]; then echo "[$run] already done"; return 0; fi
  IFS=';' read -ra U <<< "$urls"; IFS=';' read -ra M <<< "$md5s"; ok=1
  for i in $(seq 0 $((${#U[@]}-1))); do get_file "${U[$i]}" "${M[$i]}" $d || ok=0; done
  if [ $ok = 1 ]; then touch $d/.done; echo "[$run] OK $(date +%H:%M)"; else echo "[$run] GAVE UP"; fi
}
export -f get_file dl_one; export P
tail -n +2 metadata/step4_runs.tsv | awk -F'\t' '{print $5"\t"$11"\t"$12}' | parallel --will-cite -j 4 --colsep '\t' dl_one {1} {2} {3}
echo "=== summary: $(ls fastq/*/.done 2>/dev/null | wc -l) of 33 runs downloaded ==="
