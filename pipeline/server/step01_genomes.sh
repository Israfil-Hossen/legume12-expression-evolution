#!/usr/bin/env bash
# Step 1: download RefSeq genome, GFF, transcripts and proteins for all species; verify md5.
set -uo pipefail
cd /home/liuchuanfa/israfil/legume2/refs
declare -A ACC=(
  [Gmax]=GCF_000004515.6  [Pvul]=GCF_000499845.2  [Mtru]=GCF_003473485.1
  [Cari]=GCF_000331145.2  [Atha]=GCF_000001735.4
  [Vrad]=GCF_000741045.1  [Vung]=GCF_004118075.2  [Ccaj]=GCF_000230855.2
  [Adur]=GCF_000817695.3  [Lang]=GCF_001865875.1  [Psat]=GCF_024323335.1
  [Tpra]=GCF_020283565.1  [Ljap]=GCF_012489685.1 )
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all
for sp in "${!ACC[@]}"; do
  acc=${ACC[$sp]}; n=${acc#GCF_}; n=${n%.*}
  path="GCF/${n:0:3}/${n:3:3}/${n:6:3}"
  dir=$(curl -s --retry 5 "$BASE/$path/" | grep -o "${acc}_[^/\"]*" | head -1)
  if [ -z "$dir" ]; then echo "[$sp] ERROR: directory not found for $acc"; continue; fi
  mkdir -p "$sp"; cd "$sp"
  for suf in genomic.fna.gz genomic.gff.gz rna.fna.gz protein.faa.gz; do
    f="${dir}_${suf}"
    [ -s "$f" ] || wget -q --tries=5 "$BASE/$path/$dir/$f" || echo "[$sp] FAILED $f"
  done
  wget -q -O md5checksums.txt "$BASE/$path/$dir/md5checksums.txt"
  ok=$(grep -E "_(genomic.fna|genomic.gff|rna.fna|protein.faa).gz$" md5checksums.txt | sed 's#  \./#  #' | md5sum -c 2>/dev/null | grep -c ": OK")
  echo "[$sp] $acc  $dir  md5 OK: $ok/4"
  cd ..
done
echo "=== sizes ==="; du -sh */
