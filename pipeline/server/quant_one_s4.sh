#!/bin/bash
#PBS -q cpu
#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -j oe
# One biological replicate: merge lanes -> fastp -> salmon quant. Called with -v BR=<biorep>,SP=<species>.
set -uo pipefail
P=/home/liuchuanfa/israfil/legume2
source ~/miniconda3/etc/profile.d/conda.sh && conda activate legume
cd $P; mkdir -p logs/quant; exec > logs/quant/$BR.log 2>&1
runs=$(awk -F'\t' -v b="$BR" 'NR>1 && $4==b{print $5}' metadata/step4_runs.tsv)
W=$P/work/$BR; rm -rf $W; mkdir -p $W quant qc
echo "[$BR] species $SP  runs: $runs  node $(hostname)  start $(date)"
R1=""; R2=""; SE=""
for r in $runs; do
  if ls fastq/$r/*_1.fastq.gz > /dev/null 2>&1; then
    R1="$R1 $(ls fastq/$r/*_1.fastq.gz)"; R2="$R2 $(ls fastq/$r/*_2.fastq.gz)"
  else
    SE="$SE $(ls fastq/$r/*.fastq.gz | head -1)"
  fi
done
if [ -n "$R1" ]; then
  cat $R1 > $W/in_1.fq.gz; cat $R2 > $W/in_2.fq.gz
  fastp -i $W/in_1.fq.gz -I $W/in_2.fq.gz -o $W/t_1.fq.gz -O $W/t_2.fq.gz --detect_adapter_for_pe -q 20 -u 30 -l 36 --trim_poly_g -w 8 -j qc/$BR.fastp.json -h qc/$BR.fastp.html
  READS="-1 $W/t_1.fq.gz -2 $W/t_2.fq.gz"
else
  cat $SE > $W/in.fq.gz
  fastp -i $W/in.fq.gz -o $W/t.fq.gz -q 20 -u 30 -l 36 --trim_poly_g -w 8 -j qc/$BR.fastp.json -h qc/$BR.fastp.html
  READS="-r $W/t.fq.gz"
fi
salmon quant -i index/$SP -l A $READS -o quant/$BR -p 8 --seqBias --gcBias --numBootstraps 30
if [ -s quant/$BR/quant.sf ]; then touch quant/$BR/.done; echo "[$BR] DONE $(date)  $(grep -o '"percent_mapped": [0-9.]*' quant/$BR/aux_info/meta_info.json)"; else echo "[$BR] FAILED"; fi
rm -rf $W
