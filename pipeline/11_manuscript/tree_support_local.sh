#!/usr/bin/env bash
# Branch support, partitioned model and concordance factors for the species tree (requested in the pre-submission review).
# Same alignments as the registered tree (tree/concat.fa and tree/aln from the deposit archive); IQ-TREE 3.1.3 (Windows build
# from the official release, sha256 6e0e11c7...3aef). Run on the authors' laptop because the cluster nodes were unavailable.
# Per-orthogroup alignments name sequences <species>_GeneID_<id>; they are renamed to <species>, the rule concat_aln.py used.
# Usage: bash tree_support_local.sh <iqtree3.exe> <dir containing tree/> <threads>
set -euo pipefail
IQ=$1; D=$2; T=${3:-10}; SEED=20261010
cd "$D"; mkdir -p support tree/aln_sp
python - <<'EOF'
import glob, os
for f in sorted(glob.glob("tree/aln/*.aln")):
    with open(f) as fi, open(os.path.join("tree/aln_sp", os.path.basename(f)), "w") as fo:
        for line in fi:
            fo.write(">" + line[1:].split("_")[0].strip() + "\n" if line.startswith(">") else line)
print("renamed", len(glob.glob("tree/aln_sp/*.aln")), "alignments")
EOF
{
  "$IQ" --version | head -1
  echo "start $(date)  threads $T  alignments $(ls tree/aln_sp | wc -l)"
  if [ ! -s support/lg_boot.treefile ]; then
    "$IQ" -s tree/concat.fa -m LG+G4 -o Atha -B 1000 -alrt 1000 -T "$T" --seed $SEED --prefix support/lg_boot
  fi
  echo "step 1 done $(date)"
  "$IQ" -p tree/aln_sp -m MFP -o Atha -B 1000 -alrt 1000 -T "$T" --seed $SEED --prefix support/part -redo
  echo "step 2 done $(date)"
  "$IQ" -S tree/aln_sp -m MFP -T "$T" --seed $SEED --prefix support/loci
  echo "step 3 done $(date)"
  "$IQ" -t support/part.treefile --gcf support/loci.treefile -p tree/aln_sp --scfl 100 -T "$T" --seed $SEED --prefix support/concord
  echo "done $(date)"
} >> support/tree_support.log 2>&1
