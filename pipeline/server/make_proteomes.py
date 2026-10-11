#!/usr/bin/env python3
# Longest protein per gene for OrthoFinder. Header = <species>_<GeneID>.
import gzip, glob, os
P = '/home/liuchuanfa/israfil/legume2'
os.makedirs(P + '/proteomes', exist_ok=True)
for sp in ['Gmax','Pvul','Mtru','Cari','Vrad','Vung','Ccaj','Adur','Lang','Psat','Tpra','Ljap','Atha']:
    p2g = dict(l.rstrip('\n').split('\t') for l in open(f'{P}/tables/{sp}_prot2gene.tsv'))
    best, name, seq = {}, None, []
    def flush():
        if name is None or name not in p2g: return
        s = ''.join(seq); g = p2g[name]
        if g not in best or len(s) > len(best[g]): best[g] = s
    for line in gzip.open(glob.glob(f'{P}/refs/{sp}/*_protein.faa.gz')[0], 'rt'):
        if line.startswith('>'):
            flush(); name = line[1:].split()[0]; seq = []
        else:
            seq.append(line.strip())
    flush()
    with open(f'{P}/proteomes/{sp}.fa', 'w') as out:
        for g, s in best.items(): out.write(f'>{sp}_{g}\n{s}\n')
    print(f'{sp}: {len(best)} genes')
