#!/usr/bin/env bash
# Step 2: one PBS job per species -> genome/transcript prep, tx2gene, decoy-aware salmon index.
# Safe to re-run: skips species already indexed or already submitted.
P=/home/liuchuanfa/israfil/legume2
cd $P
for sp in Gmax Pvul Mtru Cari Vrad Vung Ccaj Adur Lang Psat Tpra Ljap; do
  if [ -f index/$sp/info.json ];     then echo "[$sp] index exists, skip"; continue; fi
  if [ -f jobs/index_$sp.submitted ]; then echo "[$sp] already submitted, skip"; continue; fi
  if ! grep -q "^\[$sp\].*md5 OK: 4/4" logs/step01_genomes.log 2>/dev/null; then
       echo "[$sp] genome download not verified yet, skip"; continue; fi
  wt=24:00:00; np=12; [ $sp = Psat ] && { wt=72:00:00; np=24; }
  cat > jobs/index_$sp.pbs <<JOB
#!/bin/bash
#PBS -N idx_$sp
#PBS -q cpu
#PBS -l nodes=1:ppn=$np
#PBS -l walltime=$wt
#PBS -j oe
#PBS -o $P/logs/index_$sp.log
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh && conda activate legume
cd $P/refs/$sp
echo "node \$(hostname)  start \$(date)  \$(salmon --version)"
pigz -dc *_genomic.fna.gz > genome.fa
pigz -dc *_rna.fna.gz     > transcripts.fa
python3 $P/make_maps.py *_genomic.gff.gz $P/tables/${sp}_tx2gene.tsv $P/tables/${sp}_prot2gene.tsv
nt=\$(grep -c '^>' transcripts.fa)
nm=\$(grep '^>' transcripts.fa | cut -d' ' -f1 | tr -d '>' | awk 'NR==FNR{a[\$1];next} !(\$1 in a)' $P/tables/${sp}_tx2gene.tsv - | wc -l)
echo "transcripts in fasta: \$nt   without a gene mapping: \$nm"
grep '^>' genome.fa | cut -d' ' -f1 | sed 's/^>//' > decoys.txt
cat transcripts.fa genome.fa > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -i $P/index/$sp -k 31 -p 12 --gencode
rm -f gentrome.fa genome.fa
echo "done \$(date)"
JOB
  qsub jobs/index_$sp.pbs && touch jobs/index_$sp.submitted && echo "[$sp] submitted"
done
