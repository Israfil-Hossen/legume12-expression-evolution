#!/usr/bin/env python3
# transcript->gene and protein->gene maps from a RefSeq GFF, keyed on NCBI GeneID (unique).
import gzip, re, sys
gff, out_tx, out_pr = sys.argv[1:4]
tx, pr = {}, {}
for line in gzip.open(gff, 'rt'):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9: continue
    m = re.search(r'GeneID:(\d+)', f[8])
    if not m: continue
    g = 'GeneID_' + m.group(1)
    t = re.search(r'transcript_id=([^;]+)', f[8])
    if t: tx[t.group(1)] = g
    p = re.search(r'protein_id=([^;]+)', f[8])
    if p and f[2] == 'CDS': pr[p.group(1)] = g
open(out_tx, 'w').write(''.join(f'{k}\t{v}\n' for k, v in tx.items()))
open(out_pr, 'w').write(''.join(f'{k}\t{v}\n' for k, v in pr.items()))
print(f'transcripts {len(tx)}  genes {len(set(tx.values()))}  proteins {len(pr)}')
