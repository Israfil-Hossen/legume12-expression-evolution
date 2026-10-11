#!/usr/bin/env python3
# Step 3A input preparation for one species (docs/step3_spec.md section A).
#   proteins.fa  longest protein per gene (same choice as the OrthoFinder input), header <sp>_GeneID_x
#   cds.fa       the CDS of exactly that protein (matched on protein_id), same header
#   <sp>.gff     MCScanX gene table: chromosome code, gene, start, end
# Genes whose CDS length does not match the protein (3 x aa, with or without stop) are dropped.
# Usage: prep_paralogs.py <species>
import glob, gzip, os, re, sys

P = os.environ.get('LEGUME2', '/home/liuchuanfa/israfil/legume2')
sp = sys.argv[1]
ref = f'{P}/refs/{sp}'
out = f'{P}/step3/{sp}'
os.makedirs(out, exist_ok=True)

p2g = dict(l.rstrip('\n').split('\t') for l in open(f'{P}/tables/{sp}_prot2gene.tsv'))

def fasta(path):
    name, seq = None, []
    for line in gzip.open(path, 'rt'):
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(seq)
            name, seq = line[1:].rstrip('\n'), []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)

best = {}                                   # gene -> (protein_id, protein sequence)
for h, s in fasta(glob.glob(f'{ref}/*_protein.faa.gz')[0]):
    pid = h.split()[0]
    g = p2g.get(pid)
    if g and (g not in best or len(s) > len(best[g][1])):
        best[g] = (pid, s)

cds = {}                                    # protein_id -> CDS
for h, s in fasta(f'{ref}/cds.fna.gz'):
    m = re.search(r'\[protein_id=([^\]]+)\]', h)
    if m:
        cds[m.group(1)] = s.upper()

kept = 0
with open(f'{out}/proteins.fa', 'w') as fp, open(f'{out}/cds.fa', 'w') as fc:
    for g, (pid, prot) in best.items():
        c = cds.get(pid)
        if not c or len(c) % 3 or len(c) not in (3 * len(prot), 3 * len(prot) + 3):
            continue
        name = f'{sp}_{g}'
        fp.write(f'>{name}\n{prot}\n')
        fc.write(f'>{name}\n{c}\n')
        kept += 1

genes, chrom = set(f'{sp}_{g}' for g in best), {}
with open(f'{out}/{sp}.gff', 'w') as fg:
    for line in gzip.open(glob.glob(f'{ref}/*_genomic.gff.gz')[0], 'rt'):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene':
            continue
        m = re.search(r'GeneID:(\d+)', f[8])
        if not m:
            continue
        name = f'{sp}_GeneID_{m.group(1)}'
        if name not in genes:
            continue
        code = chrom.setdefault(f[0], f'aa{len(chrom) + 1}')
        fg.write(f'{code}\t{name}\t{f[3]}\t{f[4]}\n')

print(f'{sp}: genes with longest protein {len(best)}; kept with matching CDS {kept}; sequences in gff {len(chrom)}')
